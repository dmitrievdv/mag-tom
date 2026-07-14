function printf_array(io, arr)
    for a in arr
        @printf(io, "%10.3e ", a)
    end
end

function write_array(io, arr)
    size_arr = size(arr)
    n_dim = Int64(length(size_arr))
    write(io,n_dim)
    for i_dim = 1:n_dim
        write(io, size_arr[i_dim])
    end
    for a in arr
        write(io, a)
    end
end

function read_array(io, type :: DataType)
    n_dim = read(io, Int64)
    size_arr = zeros(Int64, n_dim)
    for i_dim = 1:n_dim
        size_arr[i_dim] = read(io, Int64)
    end
    arr = zeros(size_arr...)
    for i in eachindex(arr)
        arr[i]=read(io, type)
    end
    return arr
end

function calc_and_save_kernels(star, R_ins, Ws, incs, Δv_zs, Hs, v_z_borders, n_freq, n_ζ; n_Rm = 10, n_vz = 10)
    # rm("kernels", recursive = true, force = true)
    mkpath("kernels")
    open("kernels/kernel_grid.bin", "w") do io
        write(io, v_z_borders[1], v_z_borders[2], n_freq, n_ζ)
        # @printf(io, "%8.2e %8.2e %5i %5i", v_z_borders[1], v_z_borders[2], n_freq, n_ζ)
    end

    open("kernels/models_grid.bin", "w") do io
        write_array(io, incs)
        write_array(io, R_ins)
        write_array(io, Ws)
        write_array(io, Δv_zs)
        write_array(io, Hs)
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

        print("$inc $R_in $W $Δv_z kernel calc... does file kernels/$(ker_string)_ker.bin exist? $(isfile("kernels/$(ker_string)_ker.bin"))")
        if !isfile("kernels/$(ker_string)_ker.bin")
            print("\n")
            kernel_matrix = calc_emission_kernel_matrix(star, geometry, orientation, Δv_z, n_freq, n_ζ, v_z_borders; n_Rm = n_Rm, n_vz = n_vz)
            print("\e[1A\e[2K\e[1G")
            print("\e[1A\e[2K\e[1G")
            print("\e[1A\e[2K\e[1G")
            open("kernels/$(ker_string)_ker.bin", "w") do io
                write_array(io, kernel_matrix)
            end
        end
        # print("\n")
        for i_H in eachindex(Hs)
            H = Hs[i_H]
            print("$inc $R_in $W $Δv_z $H absorption calc... does file kernels/$(ker_string)_$(i_H)_abs.bin exist? $(isfile("kernels/$(ker_string)_$(i_H)_abs.bin"))")
            if !isfile("kernels/$(ker_string)_$(i_H)_abs.bin")
                print("\n")
                absorption = calc_absorption_profile_parallel(star, geometry, orientation, v_zs, Δv_z, 0.01, 0.01, H)
                open("kernels/$(ker_string)_$(i_H)_abs.bin", "w") do io
                    write_array(io, absorption)
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

function load_kernel_models_grid(;kernels_dir = "kernels")
    io = open("$kernels_dir/models_grid.bin", "r")
    incs = read_array(io, Float64)
    R_ms = read_array(io, Float64)
    Ws = read_array(io, Float64)
    Δv_zs = read_array(io, Float64)
    Hs = read_array(io, Float64)
    close(io)
    return incs, R_ms, Ws, Δv_zs, Hs
end

function load_kernel_grid(;kernels_dir = "kernels")
    io = open("$kernels_dir/kernel_grid.bin", "r")
    v_z_borders = (read(io, Float64), read(io, Float64))
    n_freq = read(io, Int64)
    n_ζ = read(io, Int64)
    close(io)
    return v_z_borders, n_freq, n_ζ
end

function load_kernel_from_index(i_inc, i_R_m, i_W, i_Δv_z; kernels_dir = "kernels")
    io = open("$kernels_dir/$(i_inc)_$(i_R_m)_$(i_W)_$(i_Δv_z)_ker.bin", "r")
    kernel = read_array(io, Float64)
    close(io)
    return kernel
end

function load_absorption_from_index(i_inc, i_R_m, i_W, i_Δv_z, i_H; kernels_dir = "kernels")
    io = open("$kernels_dir/$(i_inc)_$(i_R_m)_$(i_W)_$(i_Δv_z)_$(i_H)_abs.bin", "r")
    absorption = read_array(io, Float64)
    close(io)
    return absorption
end

function load_kernel_models_grid_txt(;kernels_dir = "kernels")
    lines = readlines("$kernels_dir/models_grid.dat")
    incs = parse.(Float64, split(lines[1]))
    R_ms = parse.(Float64, split(lines[2]))
    Ws = parse.(Float64, split(lines[3]))
    Δv_zs = parse.(Float64, split(lines[4]))
    Hs = parse.(Float64, split(lines[5]))
    return incs, R_ms, Ws, Δv_zs, Hs
end

function load_kernel_grid_txt(;kernels_dir = "kernels")
    data = split(readline("$kernels_dir/kernel_grid.dat"))
    v_z_borders = (parse(Float64, data[1]), parse(Float64, data[2]))
    n_freq = parse(Int, data[3])
    n_ζ = parse(Int, data[4])
    return v_z_borders, n_freq, n_ζ
end

function load_kernel_from_index_txt(i_inc, i_R_m, i_W, i_Δv_z, n_freq, n_ζ; kernels_dir = "kernels")
    kernel = zeros(n_freq, n_ζ)
    open("$kernels_dir/$(i_inc)_$(i_R_m)_$(i_W)_$(i_Δv_z)_ker.dat") do io
        for i_freq = 1:n_freq
            kernel[i_freq,:] = parse.(Float64, split(readline(io)))
        end
    end
    return kernel
end

function load_absorption_from_index_txt(i_inc, i_R_m, i_W, i_Δv_z, i_H, n_freq; kernels_dir = "kernels")
    absorption = zeros(n_freq)
    absorption = parse.(Float64, split(readline("$kernels_dir/$(i_inc)_$(i_R_m)_$(i_W)_$(i_Δv_z)_$(i_H)_abs.dat")))
end

function correct_absorption_files(star;kernels_dir = "kernels")
    incs, R_ins, Ws, Δv_zs, Hs = load_kernel_models_grid(;kernels_dir)
    for i_inc in eachindex(incs), i_R_in in eachindex(R_ins), i_W in eachindex(Ws), i_Δv_z in eachindex(Δv_zs), i_H in eachindex(Hs)
        absorption = load_absorption_from_index(i_inc, i_R_in, i_W, i_Δv_z, i_H;kernels_dir)
        R_in  = R_ins[i_R_in]; W = Ws[i_W]
        geometry = DipoleGeometry(R_in, R_in+W)
        inc = incs[i_inc]
        orientation = Orientation(inc)
        H = Hs[i_H]
        absorption_correct = correct_absorption(star, geometry, orientation, 0.05, absorption, H)
        ker_string = "$(i_inc)_$(i_R_in)_$(i_W)_$(i_Δv_z)"
        open("$kernels_dir/$(ker_string)_$(i_H)_abs.bin", "w") do io
            write_array(io, absorption_correct)
        end
        print("\e[2K\e[1G$(ker_string)_$(i_H)_abs.bin $H $(findmax(absorption)) $(findmax(absorption_correct))")
    end
end