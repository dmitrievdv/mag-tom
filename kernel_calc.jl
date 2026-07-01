function calc_emission_kernel_matrix(star, geometry, orientation, Δv_z, n_freq, n_ζ, v_z_borders = (-3e7, 3e7); n_Rm = 10, n_vz = 10)
    kernel = zeros(n_freq, n_ζ)
    calc_emission_kernel_matrix!(kernel, star, geometry, orientation, Δv_z, v_z_borders; n_Rm = n_Rm, n_vz = n_vz)
    return kernel
end

function calc_emission_kernel_matrix!(kernel_matrix_raw, star, geometry, orientation, Δv_z, v_z_borders = (-3e7, 3e7); n_Rm = 10, n_vz = 10)
    n, n_ζ = size(kernel_matrix_raw)
    
    v_z_start, v_z_end = v_z_borders
    v_z_step = (v_z_end - v_z_start)/n

    n_threads = Threads.nthreads(:default)

    n_jobs = min(4n_threads, n)

    jobs = Channel{Int}(n_jobs)
    results = Channel{Tuple{Int, Int}}(n_jobs)

    println("Jobs creation")

    for i_job = 1:n_jobs
        put!(jobs, i_job)
    end

    n_threads = Threads.nthreads(:default)

    println("Work tasks, $n_threads")

    workers = []

    for i_thread = 1:n_threads
        w = Threads.@spawn for job in jobs
            # println(job)
            i_job = job
            i_vzs = [i_job:n_jobs:n;]
            # v_zs = v_z_start .+ v_z_step*i_vzs .- v_z_step/2
            # println(length(v_zs))
            # kernel = zeros(length(v_zs), n_ζ)
            for i_vz in i_job:n_jobs:n
                v_z = v_z_start + v_z_step*i_vz - v_z_step/2
                for i_ζ = 1:n_ζ
                    ζ = i_ζ/n_ζ - 1/2n_ζ
                    # println("$ζ $v_z Rm")
                    kernel_matrix_raw[i_vz,i_ζ] += calc_kernel_advanced(star, geometry, orientation, ζ, v_z, Δv_z, n_Rm = n_Rm)/π
                    # println("-$ζ $v_z Rm")
                    kernel_matrix_raw[i_vz,i_ζ] += calc_kernel_advanced(star, geometry, orientation, -ζ, v_z, Δv_z, n_Rm = n_Rm)/π
                    # println("$ζ $v_z vz")
                    kernel_matrix_raw[i_vz,i_ζ] += calc_velocity_kernel_advanced(star, geometry, orientation, ζ, v_z, Δv_z, n_vz = n_vz)/π
                    # println("-$ζ $v_z vz")
                    kernel_matrix_raw[i_vz,i_ζ] += calc_velocity_kernel_advanced(star, geometry, orientation, -ζ, v_z, Δv_z, n_vz = n_vz)/π
                end
            end
            put!(results, (Threads.threadid(), i_job))
        end
        push!(workers, w)
    end

    # println("Progress, $(Threads.nthreads(:interactive))")
    # println(istaskstarted.(workers))

    # io = open("kernel_calc.log", "w")

    print("Nothing is done. ")
    n_done = 0
    n_remain = n_jobs
    progress_worker = Threads.@spawn :interactive  while n_remain > 0
        print("Waiting...")
        i_thread, i_job = take!(results)
        i_vzs = [i_job:n_jobs:n;]
        # println(i_vzs)
        # v_z = v_z_start + v_z_step*i_v_z - v_z_step/2
        # kernel_matrix_raw[i_vzs, :] .= kernel
        # println(io, "$n_jobs $i_thread $i_v_z $v_z $kernel")
        n_remain -= 1
        n_done += 1
        print("\e[2K\e[1GTask $n_done is done (thread $i_thread). $n_remain remaining tasks. ")
    end
    wait(progress_worker)
    print("\n")

    # close(io)

    # println("All is done")
    close(jobs)
    close(results)
    # for i_v_z = 1:n
    #     v_z = v_z_start + v_z_step*i_v_z - v_z_step/2
    #     # println(v_z)
    #     for i_ζ = 1:n_ζ
    #         ζ = i_ζ/n_ζ - 1/2n_ζ
    #         kernel_matrix_raw[i_v_z, i_ζ] += TTauUtils.Profiles.calc_kernel_advanced(star, geometry, orientation, ζ, v_z, Δv_z, n_Rm = n_Rm)/π
    #         kernel_matrix_raw[i_v_z, i_ζ] += TTauUtils.Profiles.calc_kernel_advanced(star, geometry, orientation, -ζ, v_z, Δv_z, n_Rm = n_Rm)/π
    #         kernel_matrix_raw[i_v_z, i_ζ] += TTauUtils.Profiles.calc_velocity_kernel_advanced(star, geometry, orientation, ζ, v_z, Δv_z, n_vz = n_vz)/π
    #         kernel_matrix_raw[i_v_z, i_ζ] += TTauUtils.Profiles.calc_velocity_kernel_advanced(star, geometry, orientation, -ζ, v_z, Δv_z, n_vz = n_vz)/π
    #     end
    # end
end

function calc_simple_kernel_matrix(star, geometry, orientation, v_z_borders, Δv_z, n, n_ζ)
    kernel = zeros(n, n_ζ, 2)

    v_z_start, v_z_end = v_z_borders
    v_z_step = (v_z_end - v_z_start)/n

    for i_v_z = 1:n
        v_z = v_z_start + v_z_step*i_v_z - v_z_step/2
    #     # println(v_z)
        for i_ζ = 1:n_ζ
            ζ = i_ζ/n_ζ - 1/2n_ζ
            kernel[i_v_z, i_ζ, 1] += calc_kernel(star, geometry, orientation, Δv_z, ζ, v_z)/π
            kernel[i_v_z, i_ζ, 1] += calc_kernel(star, geometry, orientation, Δv_z, -ζ, v_z)/π
            kernel[i_v_z, i_ζ, 2] += calc_velocity_kernel(star, geometry, orientation, Δv_z, ζ, v_z)/π
            kernel[i_v_z, i_ζ, 2] += calc_velocity_kernel(star, geometry, orientation, Δv_z, -ζ, v_z)/π
        end
    end
    return kernel
end

function kernel_ζ_half(kernel_matrix)
    n, n_ζ = size(kernel_matrix)
    if n_ζ%2 != 0
        return kernel_matrix
    end
    half_kernel = zeros(n, n_ζ÷2)
    for i = 1:n, i_ζ = 1:n_ζ
        half_kernel[i,(i_ζ-1)÷2+1] += kernel_matrix[i,i_ζ]/2
    end
    return half_kernel
end

function kernel_vz_half(kernel_matrix)
    n, n_ζ = size(kernel_matrix)
    if n%2 != 0
        return kernel_matrix
    end
    half_kernel = zeros(n÷2, n_ζ)
    for i = 1:n, i_ζ = 1:n_ζ
        half_kernel[(i-1)÷2+1,i_ζ] += kernel_matrix[i,i_ζ]/2
    end
    return half_kernel
end

function kernel_ζ_n_half(kernel_matrix, n)
    _, n_ζ = size(kernel_matrix)
    if n_ζ%(2^n) != 0
        return kernel_matrix
    end
    half_kernel = kernel_matrix
    for i_half = 1:n
        half_kernel = kernel_ζ_half(half_kernel)
    end
    return half_kernel
end

function kernel_vz_n_half(kernel_matrix, n)
    n_freq, _ = size(kernel_matrix)
    if n_freq%(2^n) != 0
        return kernel_matrix
    end
    half_kernel = kernel_matrix
    for i_half = 1:n
        half_kernel = kernel_vz_half(half_kernel)
    end
    return half_kernel
end