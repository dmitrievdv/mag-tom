using Polynomials
using Roots
import PolynomialRoots
using LinearAlgebra
using StaticArrays
using Printf
using LeastSquaresOptim
using Statistics
using FFTW
using GLMakie

include("const.jl")
include("star.jl")
include("orientation.jl")
include("geometry.jl")
include("tom.jl")
include("kernel_calc.jl")


star = Star("test", 2, 0.8, 4000, 10)
geometry = DipoleGeometry(4,5)
orientation = Orientation(75)

n = 512

kernel = zeros(n, 2n)

# calc_kernel_advanced(star, geometry, orientation, -0.693359375, -1.17578125e7, 1e6; n_Rm = 40)
# @time kernel_simp = calc_simple_kernel_matrix(star, geometry, orientation, (-3.5e7, 3.5e7), 1e6, n, 2n)
@time calc_emission_kernel_matrix!(kernel, star, geometry, orientation, 1e6, (-3.5e7, 3.5e7); n_Rm = 20, n_vz = 20)
# kernel = kernel_ζ_n_half(kernel, 3)
# kernel = kernel_vz_n_half(kernel, 2)

fig = Figure()
ax = Axis(fig[1,1])
# ax_simp = Axis(fig[1,2])

# linkaxes!(ax, ax_simp)

# heatmap!(ax_simp, log10.(kernel_simp[:,:,1]*0.2 + kernel_simp[:,:,2]*1e5), colorrange = [-3,3])
heatmap!(ax, log10.(kernel), colorrange = [-3,3])

fig