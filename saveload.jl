function printf_array(io, arr)
    for a in arr
        @printf(io, "%10.3e ", a)
    end
end

function calc_and_save_kernels(star, R_ins, Ws, incs, Δv_zs, Hs, v_z_borders, n_freq, n_ζ; n_Rm = 10, n_vz = 10)
    # rm("kernels", recursive = true, force = true)
    mkpath("kernels")
    open("kernels/kernel_grid.dat", "w") do io
        @printf(io, "%8.2e %8.2e %5i %5i", v_z_borders[1], v_z_borders[2], n_freq, n_ζ)
    end

    open("kernels/models_grid.dat", "w") do io
        printf_array(io, incs); print(io, "\n")
        printf_array(io, R_ins); print(io, "\n")
        printf_array(io, Ws); print(io, "\n")
        printf_array(io, Δv_zs); print(io, "\n")
        printf_array(io, Hs); print(io, "\n")
    end

    v_z_start, v_z_end = v_z_borders
    v_z_step = (v_z_end - v_z_start)/n_freq
    v_zs = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n_freq]

    for i_inc in eachindex(incs), i_R_in in eachindex(R_ins), i_W in eachindex(Ws), i_Δv_z in eachindex(Δv_zs)
        W = Ws[i_W]
        R_in = R_ins[i_R_in]
        R_m = R_in + W/2
        inc = incs[i_inc]
        Δv_z = Δv_zs[i_Δv_z]
        
        geometry = DipoleGeometry(R_m - W/2, R_m + W/2)
        orientation = Orientation(inc,0,0)
        
        v_z_start, v_z_end = v_z_borders
        v_z_step = (v_z_end - v_z_start)/n_freq

        ker_string = "$(i_inc)_$(i_R_in)_$(i_W)_$(i_Δv_z)"

        print("$inc $R_in $W $Δv_z kernel calc... does file kernels/$(ker_string)_ker.dat exist? $(isfile("kernels/$(ker_string)_ker.dat"))")
        if !isfile("kernels/$(ker_string)_ker.dat")
            print("\n")
            kernel_matrix = calc_emission_kernel_matrix(star, geometry, orientation, Δv_z, n_freq, n_ζ, v_z_borders; n_Rm = n_Rm, n_vz = n_vz)
            print("\e[1A\e[2K\e[1G")
            print("\e[1A\e[2K\e[1G")
            print("\e[1A\e[2K\e[1G")
            open("kernels/$(ker_string)_ker.dat", "w") do io
                for i_freq = 1:n_freq
                    printf_array(io, kernel_matrix[i_freq,:]); print(io, "\n")
                end
            end
        end
        # print("\n")
        for i_H in eachindex(Hs)
            H = Hs[i_H]
            if !isfile("kernels/$(ker_string)_$(i_H)_abs.dat")
                print("$inc $R_in $W $Δv_z $H absorption calc... does file kernels/$(ker_string)_$(i_H)_abs.dat exist? $(isfile("kernels/$(ker_string)_$(i_H)_abs.dat"))")
                print("\n")
                absorption = calc_absorption_profile_parallel(star, geometry, orientation, v_zs, Δv_z, 0.01, 0.01, H)
                open("kernels/$(ker_string)_$(i_H)_abs.dat", "w") do io
                    printf_array(io, absorption); print(io, "\n")
                end
                print("\e[1A\e[2K\e[1G")
                print("\e[1A\e[2K\e[1G")
                print("\e[1A\e[2K\e[1G")
                print("\e[1A\e[2K\e[1G")
                print("\e[2K\e[1G")
            end
        end
        print("\e[1A\e[2K\e[1G")
    end
    println("Done!")
end

function load_kernel_models_grid(kernels_dir = "kernels")
    lines = readlines("$kernels_dir/models_grid.dat")
    incs = parse.(Float64, split(lines[1]))
    R_ms = parse.(Float64, split(lines[2]))
    Ws = parse.(Float64, split(lines[3]))
    Δv_zs = parse.(Float64, split(lines[4]))
    Hs = parse.(Float64, split(lines[5]))
    return incs, R_ms, Ws, Δv_zs, Hs
end

function load_kernel_grid(kernels_dir = "kernels")
    data = split(readline("$kernels_dir/kernel_grid.dat"))
    v_z_borders = (parse(Float64, data[1]), parse(Float64, data[2]))
    n_freq = parse(Int, data[3])
    n_ζ = parse(Int, data[4])
    return v_z_borders, n_freq, n_ζ
end

function load_kernel_from_index(i_inc, i_R_m, i_W, i_Δv_z, n_freq, n_ζ, kernels_dir = "kernels")
    kernel = zeros(n_freq, n_ζ)
    open("$kernels_dir/$(i_inc)_$(i_R_m)_$(i_W)_$(i_Δv_z)_ker.dat") do io
        for i_freq = 1:n_freq
            kernel[i_freq,:] = parse.(Float64, split(readline(io)))
        end
    end
    return kernel
end

function load_absorption_from_index(i_inc, i_R_m, i_W, i_Δv_z, i_H, n_freq, kernels_dir = "kernels")
    absorption = zeros(n_freq)
    absorption = parse.(Float64, split(readline("$kernels_dir/$(i_inc)_$(i_R_m)_$(i_W)_$(i_Δv_z)_$(i_H)_abs.dat")))
end

function load_Δv_z(kernels_dir = "kernels")
    data = readlines("$kernels_dir/models_grid.dat")
    return parse(Float64, data[end])
end