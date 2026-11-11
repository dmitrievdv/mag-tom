include("include.jl")


function main()
star = Star("PDS70", 1.26, 0.76, 4040, 20)

incs = [40:2:60.0;]
R_ins = [2.0:0.3:5;]
Ws = [0.01; 0.03; 0.05; 0.1; 0.3:0.2:1.1;]
Δv_zs = [2e6:0.4e6:6e6;]#[2e6:0.2e6:4.0e6;]
Hs = 10 .^ [0.0:0.1:1;]

v_z_borders = (-3.0e7, 6e7)
n_freq = 256
n_ζ = 512

@time calc_and_save_kernels(star, R_ins, Ws, incs, Δv_zs, Hs, v_z_borders, n_freq, n_ζ; n_Rm = 20, n_vz = 20)
# Ws = [0.1,0.5,1.0]
# @time redo_models_grid(star, incs, R_ins, Ws, Δv_zs, Hs; n_Rm = 20, n_vz = 20)

v_zs = Float64[]
profile = Float64[]
obs_file = "pds70obs-archive_2-korg"
open("$obs_file.dat", "r") do io
    for line in readlines(io)
        v_z, r = parse.(Float64, split(line))
        push!(v_zs, v_z)
        push!(profile, r)
    end
end

@time residuals = fit_known(reshape(v_zs*1e5, (1,length(v_zs))), reshape(profile, (1,length(v_zs))))

open("residuals/residuals_$obs_file.bin", "w") do io
    write_array(io, residuals[1,:,:,:,:,:])
end

v_zs = Float64[]
profile = Float64[]
obs_file = "pds70obs-archive_2"
open("$obs_file.dat", "r") do io
    for line in readlines(io)
        v_z, r = parse.(Float64, split(line))
        push!(v_zs, v_z)
        push!(profile, r)
    end
end

@time residuals = fit_known(reshape(v_zs*1e5, (1,length(v_zs))), reshape(profile, (1,length(v_zs))))

open("residuals/residuals_$obs_file.bin", "w") do io
    write_array(io, residuals[1,:,:,:,:,:])
end

v_zs = Float64[]
profile = Float64[]
obs_file = "pds70obs-archive_2"
open("$obs_file.dat", "r") do io
    for line in readlines(io)
        v_z, r = parse.(Float64, split(line))
        r = r - 0.6*exp(-(v_z^2/(2*20^2)))
        push!(v_zs, v_z)
        push!(profile, r)
    end
end

@time residuals = fit_known(reshape(v_zs*1e5, (1,length(v_zs))), reshape(profile, (1,length(v_zs))))

open("residuals/residuals_$obs_file-gauss.bin", "w") do io
    write_array(io, residuals[1,:,:,:,:,:])
end

# Ws = [0.1,0.5,1.0]
# @time redo_models_grid(star, incs, R_ins, Ws, Δv_zs, Hs; n_Rm = 20, n_vz = 20)

# geometry_true = DipoleGeometry(4,5)
# orientation_true = Orientation(60)

# source_known = open("hart_to_fit/mag_3-4_60_53_S.bin", "r") do io
#     read_array(io, Float64)
# end

# files = ["3-4_15_32", "3-4_35_32", "3-4_55_32", "3-4_75_32",
#          "3-4_15_42", "3-4_35_42", "3-4_55_42", "3-4_75_42",
#          "3-4_15_52", "3-4_35_52", "3-4_55_52", "3-4_75_52",
#          "3-4_15_43", "3-4_35_43", "3-4_55_43", "3-4_75_43",
#          "3-4_15_53", "3-4_35_53", "3-4_55_53", "3-4_75_53",
#          "3-4_15_74", "3-4_35_74", "3-4_55_74", "3-4_75_74"]

# mkpath("residuals_H")

# io = open("hart_to_fit/mag_$(files[1])_P.bin", "r") 
# v_zs = reverse(read_array(io, Float64))'
# profiles = reverse(read_array(io, Float64))'
# close(io)

# for file in files[2:end]
#     io = open("hart_to_fit/mag_$(file)_P.bin", "r")
#     v_zs = vcat(v_zs, reverse(read_array(io, Float64))')
#     profiles = vcat(profiles, reverse(read_array(io, Float64))')
#     close(io)
# end
    
# Hs = [3.705, 5.788, 7.133, 1.749, 2.070, 1.671]

# recalc_and_save_absorption(star, Hs; kernels_dir = "kernels_bin")

# residuals = fit_known_H(v_zs, profiles, 4; absorption_dir = "kernels_abs", kernels_dir = "kernels_bin")

# for i_file in eachindex(files)
#     file = files[i_file]
#     open("residuals_H/residuals_$file.bin", "w") do io
#         write_array(io, residuals[i_file, :,:,:,:,:])
#     end
# end
end

main()