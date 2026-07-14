include("include.jl")


function main()
star = Star("CTTS", 2, 0.8, 4000, 15)

# incs = [10:10:80.0;]
# R_ins = [3.0,4,5]
# Ws = [0.5,1,1.5]
# Δv_zs = [1.5e6]
# Hs = [1:3.0:10.0;]

# v_z_borders = (-3.5e7, 3.5e7)
# n_freq = 256
# n_ζ = 512

# calc_and_save_kernels(star, R_ins, Ws, incs, Δv_zs, Hs, v_z_borders, n_freq, n_ζ; n_Rm = 20, n_vz = 20)

# geometry_true = DipoleGeometry(4,5)
# orientation_true = Orientation(60)

source_known = open("hart_to_fit/mag_3-4_60_53_S.bin", "r") do io
    read_array(io, Float64)
end

files = ["3-4_15_43", "3-4_35_43", "3-4_55_43", "3-4_75_43",
         "3-4_15_52", "3-4_35_52", "3-4_55_52", "3-4_75_52"]

mkpath("residuals")

io = open("hart_to_fit/mag_$(files[1])_P.bin", "r") 
v_zs = reverse(read_array(io, Float64))'
profiles = reverse(read_array(io, Float64))'
close(io)

for file in files[2:end]
    io = open("hart_to_fit/mag_$(file)_P.bin", "r")
    v_zs = vcat(v_zs, reverse(read_array(io, Float64))')
    profiles = vcat(profiles, reverse(read_array(io, Float64))')
    close(io)
end
    
    

residuals = fit_known(v_zs, profiles; kernels_dir = "kernels_bin")

for i_file in eachindex(files)
    file = files[i_file]
    open("residuals/residuals_$file.bin", "w") do io
        write_array(io, residuals[i_file, :,:,:,:,:])
    end
end
end

main()