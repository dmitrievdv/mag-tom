abstract type AbstractGeometry end

spheretogrid(geometry :: AbstractGeometry, r, θ) = (r, θ)
gridtosphere(geometry :: AbstractGeometry, grid_1, grid_2) = (grid_1, grid_2)

spheretogridproper(geometry :: AbstractGeometry, r, θ) = spheretogrid(geometry, r, θ)
gridtosphereproper(geometry :: AbstractGeometry, grid_x, grid_y) = gridtosphere(geometry, grid_1, grid_2)

struct DipoleGeometry <: AbstractGeometry
    r_mi :: Float64
    r_mo :: Float64
end

function truncatebystar(R, borders)
    n = length(borders)
    borders_truncated = borders
    if R <= 1 
        truncate = √(1 - R^2)
        borders_truncated = borders[borders .>  truncate]
        n_truncated = length(borders_truncated)
        if (n - n_truncated)%2 == 1
            return [truncate; borders_truncated]
        else
            return borders_truncated
        end
    end
    return borders
end

function calcborders(xs :: Array{Float64, 1}, ys :: Array{Float64, 1}, 
                     geometry :: AbstractGeometry, orientation :: Orientation)
    borders = fill([0.0], length(xs))
    for i ∈ 1:length(xs)
        x, y = xs[i], ys[i]
        borders[i] = calcborders(x, y, geometry, orientation)
    end
    return borders
end

"""
    calcborders(x, y, geometry :: AbstractGeometry, orientation :: Orientation)

Returns array of z coordinates (z axis is along the line of sight) of points where the line
of sight with coordinates `x` and `y` in the picture plane crosses the borders of `geometry`,
which orientation is specified in `orientation`. `x` and `y` can be arrays of coordinates
(with the same size), then an array of arrays of z coordinates is returned.
"""
function calcborders(x :: Real, y :: Real, magnetosphere :: DipoleGeometry, orientation :: Orientation)
    ψ = orientation.ψ
    R = √(x^2 + y^2)
    R_d = orientation.dipole_axis[1]*x + orientation.dipole_axis[2]*y
    ρ = R^2 - R_d^2
    P0 = Polynomial([R^6, 0, 3R^4, 0, 3R^2, 0, 1])
    z_d = orientation.dipole_axis[3]
    Pmag = Polynomial([ρ^2, -4ρ*R_d*z_d, 4R_d^2*z_d^2 + 2ρ - 2ρ*z_d^2, 
    4R_d*z_d^3 - 4R_d*z_d, 1 + z_d^4 - 2z_d^2])
    P_in = P0 - magnetosphere.r_mi^2*Pmag
    P_out = P0 - magnetosphere.r_mo^2*Pmag
    all_roots = [roots(P_in); roots(P_out)]
    borders = truncatebystar(R, sort(@. real(all_roots[imag(all_roots) ≈ 0])))
    return borders
end

"""
    gridtosphere(geometry :: AbstractGeometry, r_m, t)

Returns spherical coordinates `(r, θ)`
"""
function gridtosphere(geometry :: DipoleGeometry, r_m, t)
    θ_star = asin(√(1/r_m))
    θ = θ_star + (π/2 - θ_star)*t 
    r = r_m*sin(θ)^2
    return r, θ
end

function spheretogrid(geometry :: DipoleGeometry, r, θ)
    r_m = r/sin(θ)^2
    θ_star = asin(√(1/r_m))
    t = (π/2 - abs(π/2 - θ) - θ_star)/(π/2 - θ_star) # if θ ≤ π/2: (θ - θ_star)/(π/2 - θ_star)
                                                     # if θ > π/2: (π - θ - θ_star)/(π/2 - θ_star)
    return r_m, t
end

function spheretogridproper(geometry :: DipoleGeometry, r, θ)
    r_m = r/sin(θ)^2
    θ_star = asin(√(1/r_m))
    ζ = (θ - π/2)/(θ_star - π/2)
    return r_m, ζ
end

function gridtosphereproper(geometry :: DipoleGeometry, r_m, ζ)
    θ = π/2 - ζ*acos(√(1/r_m))
    r = r_m*sin(θ)^2
    return r, θ
end

function dipolegridinversejacobianmatrix(R_m, ζ, ϕ)
    dθ_dζ = -acos(√(1/R_m))
    θ = π/2 + ζ*dθ_dζ
    sinθ = sin(θ); sinϕ = sin(ϕ)
    cosθ = cos(θ); cosϕ = cos(ϕ)

    ξ_2 = 3*sinθ^2-2
    ξ_1 = √(1/(R_m - 1))

	return SA[(ξ_2*sinϕ)/sinθ^3                                     -ξ_2*cosϕ/sinθ^3                                        3*cosθ/sinθ^2;
	          (ξ_1*ζ*ξ_2+2*cosθ*sinθ)*sinϕ/(2*R_m*sinθ^3*(θ-π/2)/ζ) -(ξ_1*ζ*ξ_2+2*cosθ*sinθ)*cosϕ/(2*R_m*sinθ^3*(θ-π/2)/ζ) (3*ξ_1*ζ*cosθ/sinθ-2)/(2*R_m*sinθ*(θ-π/2)/ζ);
	          cosϕ/(R_m*sinθ^3)                                       sinϕ/(R_m*sinθ^3)                                        0]
end

function dipolegridjacobianmatrix(R_m, ζ, ϕ)
    dθ_dζ = -acos(√(1/R_m))
    dθ_dRm = -ζ/2R_m*√(1/(R_m-1))
    θ = π/2 + ζ*dθ_dζ
    sinθ = sin(θ); sinϕ = sin(ϕ)
    cosθ = cos(θ); cosϕ = cos(ϕ)
    
    
    dx_dϕ = R_m*sinθ^3*cosϕ
    dy_dϕ = R_m*sinθ^3*sinϕ
    dz_dϕ = 0.0

    dx_dζ = 3R_m*dθ_dζ*sinθ^2*cosθ*sinϕ
    dy_dζ =  -3R_m*dθ_dζ*sinθ^2*cosθ*cosϕ
    dz_dζ =  R_m*dθ_dζ*sinθ*(2cosθ^2 - sinθ^2)

    dx_dRm = (sinθ + 3R_m*dθ_dRm*cosθ)*sinθ^2*sinϕ
    dy_dRm = -(sinθ + 3R_m*dθ_dRm*cosθ)*sinθ^2*cosϕ
    dz_dRm = sinθ^2*cosθ + R_m*dθ_dRm*sinθ*(2cosθ^2 - sinθ^2)

    return SMatrix{3,3,Float64}(dx_dRm, dy_dRm, dz_dRm, dx_dζ, dy_dζ, dz_dζ, dx_dϕ, dy_dϕ, dz_dϕ)
end

function gridinversejacobianmatrix(geometry :: DipoleGeometry, r_m, ζ, ϕ)
    dipolegridinversejacobianmatrix(r_m, ζ, ϕ)
end

function gridjacobianmatrix(geometry :: DipoleGeometry, r_m, ζ, ϕ)
    dipolegridjacobianmatrix(r_m, ζ, ϕ)
end

function workgridjacobianmatrix(geometry :: DipoleGeometry, r_m, ζ, ϕ, orientation :: Orientation)
    J = gridjacobianmatrix(geometry, r_m, ζ, ϕ)
    return orientation.w_from_d*J
end

function workgridinversejacobianmatrix(geometry :: DipoleGeometry, r_m, ζ, ϕ, orientation :: Orientation)
    inv_J = gridinversejacobianmatrix(geometry, r_m, ζ, ϕ)
    return inv_J*orientation.d_from_w
end

function gridstreamlinevelocity(geometry :: DipoleGeometry, r_m, t)
    θ_star = asin(√(1/r_m))
    Δθ = (π/2 - θ_star)
    θ = θ_star + Δθ*t 
    r = r_m*sin(θ)^2
    ξ = sign(cos(θ))*√(4 - 3*sin(θ)^2)
    return sin(θ)/(ξ*Δθ*r)
end

innerradii(geometry :: DipoleGeometry) = geometry.r_mi
outerradii(geometry :: DipoleGeometry) = geometry.r_mo

function isinside(geometry :: DipoleGeometry, decart_pos)
    r = norm(decart_pos)
    if r < 1
        return false
    end
    cosθ = decart_pos[3]/r
    sin²θ = 1 - cosθ^2
    r_m = r/sin²θ
    if r_m < geometry.r_mi || r_m > geometry.r_mo
        return false
    else
        return true
    end
end

function geometrysize(geometry :: DipoleGeometry)
    return geometry.r_mo
end