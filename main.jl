using Polynomials
using Roots
import PolynomialRoots
using LinearAlgebra
using StaticArrays
using Printf
using LeastSquaresOptim
using Statistics
using FFTW
# using GLMakie

include("const.jl")
include("star.jl")
include("orientation.jl")
include("geometry.jl")
include("grids.jl")
include("tom.jl")
include("kernel_calc.jl")
include("saveload.jl")


star = Star("test", 2, 0.8, 4000, 10)
geometry = DipoleGeometry(4,5)
orientation = Orientation(60)

n = 512
v_z_borders = (-3.5e7, 3.5e7) 

kernel = zeros(n, 2n)

# @time kernel_simp = calc_simple_kernel_matrix(star, geometry, orientation, v_z_borders, 1e6, n, 2n)
R_ins = [3]
Ws = [1]
incs = [30,45]
Δv_zs = [1e6]
Hs = [1.0,3.0]

precompile(calc_and_save_kernels, tuple(typeof.([R_ins, Ws, incs, Δv_zs, Hs, v_z_borders, n, 2n])...))
@time calc_and_save_kernels(star, R_ins, Ws, incs, Δv_zs, Hs, v_z_borders, n, 2n)

# kernel = zeros(8, 16)
# @time calc_emission_kernel_matrix!(kernel, star, geometry, orientation, 1e6, v_z_borders; n_Rm = 20, n_vz = 20)
# kernel = zeros(n, 2n)
# # @time calc_emission_kernel_matrix!(kernel, star, geometry, orientation, 1e6, v_z_borders; n_Rm = 20, n_vz = 20)

# n_freq, n_ζ = size(kernel)
# ζs = [2i_ζ/n_ζ - 1/n_ζ for i_ζ = 1:(n_ζ÷2)]
# v_z_start, v_z_end = v_z_borders
# v_z_step = (v_z_end - v_z_start)/n_freq
# v_zs = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n_freq]
# @time calc_absorption_profile_parallel(star, geometry, orientation, [0.0], 1e6, 0.5, 0.5, 1.0)
# @time calc_absorption_profile_parallel(star, geometry, orientation, v_zs, 1e6, 0.01, 0.01, 1.0)


# kernel = kernel_ζ_n_half(kernel, 1)
# kernel = kernel_vz_n_half(kernel, 2)

# 

# fig = Figure()
# ax = Axis(fig[1,1])
# # ax_simp = Axis(fig[1,2])

# # linkaxes!(ax, ax_simp)

# # heatmap!(ax_simp, log10.(kernel_simp[:,:,1]*0.2 + kernel_simp[:,:,2]*1e5), colorrange = [-3,3])
# heatmap!(ax, v_zs, ζs, log10.(kernel), colorrange = [-3,3])

# screen = display(fig)
# wait(screen) 