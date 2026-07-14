function calc_poloidal_velocity_proper_component(v_esc, R_m, ζ)
    dθ_dζ = -acos(√(1/R_m))
    θ = π/2 + dθ_dζ*ζ
    return -v_esc/√(R_m) * cos(θ) /(R_m*sin(θ)^2*√(4-3*sin(θ)^2)*dθ_dζ)
end

function calc_poloidal_radial_velocity(v_esc, R_m, ζ, dz_dζ)
    v_ζ = calc_poloidal_velocity_proper_component(v_esc, R_m, ζ)
    return -v_ζ * dz_dζ
end

function calc_rotational_velocity_proper_component(v_eq)
    return v_eq
end

function calc_rotational_radial_velocity(v_eq, dz_dϕ)
    v_ϕ = calc_rotational_velocity_proper_component(v_eq)
    return -v_ϕ * dz_dϕ
end

function calc_radial_velocity(v_esc, v_eq, R_m, ζ, dz_dζ, dz_dϕ)
    return calc_poloidal_radial_velocity(v_esc, R_m, ζ, dz_dζ) + calc_rotational_radial_velocity(v_eq, dz_dϕ)
end

function calc_radial_velocity(geometry, orientation, v_esc, v_eq, R_m, ζ, ϕ)
    J = workgridjacobianmatrix(geometry, R_m, ζ, ϕ, orientation)
    dz_dζ = J[3,2]; dz_dϕ = J[3,3]
    return calc_radial_velocity(v_esc, v_eq, R_m, ζ, dz_dζ, dz_dϕ)
end

function calc_radial_velocity(star :: AbstractStar, geometry, orientation, R_m, ζ, ϕ)
    J = workgridjacobianmatrix(geometry, R_m, ζ, ϕ, orientation)
    v_esc = starvescincms(star); v_eq = starveqincms(star)
    dz_dζ = J[3,2]; dz_dϕ = J[3,3]
    return calc_radial_velocity(v_esc, v_eq, R_m, ζ, dz_dζ, dz_dϕ)
end

function calc_radial_velocity_gradient(star :: AbstractStar, geometry, orientation, R_m, ζ, ϕ)
    J_d = gridjacobianmatrix(geometry, R_m, ζ, ϕ)
    M = orientation.w_from_d

    v_esc = starvescincms(star)
    v_eq = starveqincms(star)

    calc_radial_velocity_gradient(v_esc, v_eq, R_m, ζ, ϕ, J_d, M)
end

function calc_radial_velocity_gradient(v_esc, v_eq, R_m, ζ, ϕ, J_d, M)
    z_transform = M[3,:]

    dz_dζ  = z_transform ⋅ J_d[:,2]
    dz_dϕ = z_transform ⋅ J_d[:,3]

    v_ζ = calc_poloidal_velocity_proper_component(v_esc, R_m, ζ)
    v_ϕ = calc_rotational_velocity_proper_component(v_eq)
    
    dθ_dζ = -acos(√(1/R_m))
    dθ_dRm = -ζ/2R_m*√(1/(R_m - 1))
    θ = π/2 + dθ_dζ*ζ

    cosθ = cos(θ); sinϕ = sin(ϕ)
    sinθ = sin(θ); cosϕ = cos(ϕ)
    cotθ = cot(θ)
    tanθ = tan(θ)

    dvζ_dζ =  -v_ζ*(2*cotθ + tanθ/(4 - 3*sinθ^2))*dθ_dζ
    dvζ_dRm = -v_ζ*(3/2R_m + dθ_dRm*(1/ζ/dθ_dζ + 2*cotθ + tanθ/(4 - 3*sinθ^2)))

    d²xd_dζ² =  3R_m*dθ_dζ^2 * (2*cosθ^2 - sinθ^2)*sinθ*sinϕ
    d²yd_dζ² = -3R_m*dθ_dζ^2 * (2*cosθ^2 - sinθ^2)*sinθ*cosϕ
    d²zd_dζ² = R_m*dθ_dζ^2*cosθ * (2*cosθ^2 - 7*sinθ^2)

    d²xd_dζdRm =  3R_m*((dθ_dζ/R_m + dθ_dRm/ζ)sinθ*cosθ + dθ_dζ*dθ_dRm*(2*cosθ^2 - sinθ^2))*sinθ*sinϕ
    d²yd_dζdRm = -3R_m*((dθ_dζ/R_m + dθ_dRm/ζ)sinθ*cosθ + dθ_dζ*dθ_dRm*(2*cosθ^2 - sinθ^2))*sinθ*cosϕ
    d²zd_dζdRm = R_m*((dθ_dζ/R_m + dθ_dRm/ζ)sinθ*(2*cosθ^2 - sinθ^2) + dθ_dζ*dθ_dRm*cosθ*(2*cosθ^2 - 7*sinθ^2))

    d²_dζ² = SA[d²xd_dζ², d²yd_dζ², d²zd_dζ²]
    d²_dζdRm = SA[d²xd_dζdRm, d²yd_dζdRm, d²zd_dζdRm]

    d²_dϕ² = SA[-J_d[2,3], J_d[1,3], 0.0] #   [-dyd_dϕ,  dxd_dϕ,  0.0]
    d²_dϕdζ = SA[-J_d[2,2], J_d[1,2], 0.0] #  [-dyd_dζ,  dxd_dζ,  0.0]
    d²_dϕdRm = SA[-J_d[2,1], J_d[1,1], 0.0] # [-dyd_dRm, dxd_dRm, 0.0]

    d²z_dϕ² = z_transform ⋅ d²_dϕ²
    d²z_dϕdRm = z_transform ⋅ d²_dϕdRm
    d²z_dϕdζ = z_transform ⋅ d²_dϕdζ
    
    d²z_dζdϕ = d²z_dϕdζ
    d²z_dζ² = z_transform ⋅ d²_dζ²
    d²z_dζdRm = z_transform ⋅ d²_dζdRm

    dvz_dRm = -dvζ_dRm*dz_dζ - v_ζ*d²z_dζdRm - v_ϕ*d²z_dϕdRm
    dvz_dζ = -dvζ_dζ*dz_dζ - v_ζ*d²z_dζ² - v_ϕ*d²z_dϕdζ
    dvz_dϕ = - v_ζ*d²z_dζdϕ - v_ϕ*d²z_dϕ²
    return SA[dvz_dRm, dvz_dζ, dvz_dϕ]
end

function calc_dipole_magnetosphere_radius_los_derivative(R_m, ζ, ϕ, cartesian_dipole_from_work_transform :: AbstractMatrix)
    dθ_dζ = -acos(√(1/R_m))
    θ = π/2 + dθ_dζ*ζ
    r = R_m*sin(θ)^2
    
    x_d = r*sin(θ)*sin(ϕ)
    y_d = -r*sin(θ)*cos(ϕ)
    z_d = r*cos(θ)

    dipole_cartesian_gradient = SA[x_d/r/sin(θ)^2*(3 - 2/sin(θ)), y_d/r/sin(θ)^2*(3 - 2/sin(θ)), 3z_d/r/sin(θ)^2]
    dRm_dz = dipole_cartesian_gradient ⋅ cartesian_dipole_from_work_transform[:,3]
    return dRm_dz
end

function calc_dipole_magnetosphere_radius_los_derivative(R_m, ζ, ϕ, orientation :: Orientation)
    calc_dipole_magnetosphere_radius_los_derivative(R_m, ζ, ϕ, orientation.d_from_w)
end

function calc_dipole_radial_velocity_los_derivative(R_m, ζ, ϕ, v_esc, v_eq, inv_J, J_d, w_from_d)
    dθ_dζ = -acos(√(1/R_m))
    θ = π/2 + dθ_dζ*ζ
    r = R_m*sin(θ)^2

    dvz_d = calc_radial_velocity_gradient(v_esc, v_eq, R_m, ζ, ϕ, J_d, w_from_d)

    return dvz_d ⋅ inv_J[:,3]
end

function calc_dipole_radial_velocity_los_derivative(inv_J, dvz_d)
    return dvz_d ⋅ inv_J[:,3]
end

function calc_dipole_radial_velocity_los_derivative(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, R_m, ζ, ϕ)
    inv_J = workgridinversejacobianmatrix(geometry, R_m, ζ, ϕ, orientation)
    J_d = gridjacobianmatrix(geometry, R_m, ζ, ϕ)
    calc_dipole_radial_velocity_los_derivative(R_m, ζ, ϕ, starvescincms(star), starveqincms(star), inv_J, J_d, orientation.w_from_d)
end

function calc_picture_plane_jacobian(J_d, M, dvz_d)
    J = M*J_d
    jacobian_correction = SA[J[1,3]*dvz_d[1]/dvz_d[3] J[1,3]*dvz_d[2]/dvz_d[3]; J[2,3]*dvz_d[1]/dvz_d[3] J[2,3]*dvz_d[2]/dvz_d[3]]

    return J[1:2,1:2] - jacobian_correction
end

function calc_picture_plane_jacobian(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation,
                                R_m, ζ, ϕ)
    # J = workgridjacobianmatrix(geometry, R_m, ζ, ϕ, orientation)
    J_d = gridjacobianmatrix(geometry, R_m, ζ, ϕ)
    M = orientation.w_from_d

    v_esc = starvescincms(star)
    v_eq = starveqincms(star)

    dvz_d = calc_radial_velocity_gradient(v_esc, v_eq, R_m, ζ, ϕ, J_d, M)

    calc_picture_plane_jacobian(J_d, M, dvz_d)
end

function calc_picture_plane_velocity_jacobian(J_d, M, dvz_d)
    J = M*J_d

    return SA[J[1,2]-J[1,3]*dvz_d[2]/dvz_d[3] J[1,3]/dvz_d[3]; J[2,2]-J[2,3]*dvz_d[2]/dvz_d[3] J[2,3]/dvz_d[3]]
end

function calc_picture_plane_velocity_jacobian(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation,
                                R_m, ζ, ϕ)
    # J = workgridjacobianmatrix(geometry, R_m, ζ, ϕ, orientation)
    J_d = gridjacobianmatrix(geometry, R_m, ζ, ϕ)
    M = orientation.w_from_d
    J = M*J_d

    v_esc = starvescincms(star)
    v_eq = starveqincms(star)

    dvz_d = calc_radial_velocity_gradient(v_esc, v_eq, R_m, ζ, ϕ, J_d, M)

    return SA[J[1,2]-J[1,3]*dvz_d[2]/dvz_d[3] J[1,3]/dvz_d[3]; J[2,2]-J[2,3]*dvz_d[2]/dvz_d[3] J[2,3]/dvz_d[3]]
end

function calc_picture_plane_jacobian_ζϕ(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation,
                                R_m, ζ, ϕ)
    # J = workgridjacobianmatrix(geometry, R_m, ζ, ϕ, orientation)
    J_d = gridjacobianmatrix(geometry, R_m, ζ, ϕ)
    M = orientation.w_from_d
    J = M*J_d

    v_esc = starvescincms(star)
    v_eq = starveqincms(star)

    dvz_d = calc_radial_velocity_gradient(v_esc, v_eq, R_m, ζ, ϕ, J_d, M)

    jacobian_correction = SA[J[1,1]*dvz_d[2]/dvz_d[1] J[1,1]*dvz_d[3]/dvz_d[1]; J[2,1]*dvz_d[2]/dvz_d[1] J[2,1]*dvz_d[3]/dvz_d[1]]

    return J[1:2,2:3] - jacobian_correction
end

function find_v_z(star :: AbstractStar, geometry, orientation, R_m, ζ, v_z)
    v_esc = starvescincms(star); v_eq = starveqincms(star)
    v_ζ = calc_poloidal_velocity_proper_component(v_esc, R_m, ζ)
    v_ϕ = calc_rotational_velocity_proper_component(v_eq)

    # calculating the jacobian without azimuthal component
    J_d_1 = gridjacobianmatrix(geometry, R_m, ζ, π/4) # cos(π/4) = sin(π/4) = 1/√2, see next line
    sqrt2 = √2
    ϕ_remove =  SMatrix{3,3,Float64}(sqrt2, sqrt2, 1.0, sqrt2, sqrt2, 1.0, sqrt2, sqrt2, 1.0)# SA[√2 √2 √2; √2 √2 √2; 1 1 1]
    J_d = map(*, J_d_1, ϕ_remove)
    M = orientation.w_from_d

    a = v_ζ*M[3,1]*J_d[1,2] + v_ϕ*M[3,2]*J_d[2,3]
    b = v_ζ*M[3,2]*J_d[2,2] + v_ϕ*M[3,1]*J_d[1,3]
    c = v_z + v_ζ*M[3,3]*J_d[3,2] + v_ϕ*M[3,3]*J_d[3,3]
    D = (b^2*c^2 + (a^2 - c^2)*(a^2 + b^2))
    if D < 0
        return -1e10, -1e10
    else
        cosϕ₀ = -b*c/(a^2 + b^2)
        cosϕ₊ = cosϕ₀ + √(D)/(a^2 + b^2)
        cosϕ₋ = cosϕ₀ - √(D)/(a^2 + b^2)
        sinϕ₊ = -(cosϕ₊*b + c)/a
        sinϕ₋ = -(cosϕ₋*b + c)/a

        acosϕ₊ = if abs(cosϕ₊) ≤ 1
            acos(cosϕ₊)
        else
            acos(cosϕ₊/abs(cosϕ₊))
        end

        acosϕ₋ = if abs(cosϕ₋) ≤ 1
            acos(cosϕ₋)
        else
            acos(cosϕ₋/abs(cosϕ₋))
        end

        ϕ₊ = if sinϕ₊ ≥ 0
            acosϕ₊
        else
            2π - acosϕ₊
        end

        ϕ₋ = if sinϕ₋ ≥ 0
            acosϕ₋
        else
            2π - acosϕ₋
        end

        return ϕ₋, ϕ₊
    end
end

function dipole_decart_coordinates(R_m, ζ, ϕ)
    dθ_dζ = -acos(√(1/R_m))
    θ = π/2 + dθ_dζ*ζ
    x_d = R_m*sin(θ)^3*sin(ϕ)
    y_d = -R_m*sin(θ)^3*cos(ϕ)
    z_d = R_m*sin(θ)^2*cos(θ)
    return SA[x_d, y_d, z_d]
end


# calc_radial_velocity(geometry, orientation, v_esc, v_eq, R_m, ζ, ϕ)

function calc_v_z(x,y,z,orientation,geometry :: DipoleGeometry,star)
    xyz = SA[x,y,z]
    xyz_d = orientation.d_from_w * xyz
    r = √(x^2 + y^2 + z^2)
    θ = acos(xyz_d[3]/r)
    ϕ = atan(xyz_d[1],-xyz_d[2])
    R_m = r/sin(θ)^2
    θ_star = asin(√(1/R_m))
    ζ = (θ - π/2)/(θ_star - π/2)
    -calc_radial_velocity(star, geometry, orientation, R_m, ζ, ϕ)
end

function calcborders_prealloc(x :: Real, y :: Real, 
                                magnetosphere :: DipoleGeometry, orientation :: Orientation)
    ψ = orientation.ψ
    R = √(x^2 + y^2)
    R_d = orientation.dipole_axis[1]*x + orientation.dipole_axis[2]*y
    ρ = R^2 - R_d^2
    P0 = [R^6, 0, 3R^4, 0, 3R^2, 0, 1]
    z_d = orientation.dipole_axis[3]
    Pmag = [ρ^2, -4ρ*R_d*z_d, 4R_d^2*z_d^2 + 2ρ - 2ρ*z_d^2, 
    4R_d*z_d^3 - 4R_d*z_d, 1 + z_d^4 - 2z_d^2, 0, 0]
    P_in = P0 - magnetosphere.r_mi^2*Pmag
    P_out = P0 - magnetosphere.r_mo^2*Pmag
    all_roots = [PolynomialRoots.roots(P_in); PolynomialRoots.roots(P_out)]
    borders = GeometryAndOrientations.truncatebystar(R, sort(@. real(all_roots[abs(imag(all_roots)) < 1e-8])))
end

function calc_dipole_line_intersections!(zs, roots_arr :: Vector{C}, poly :: Vector{C}, x :: Real, y :: Real, 
                                R_m, orientation :: Orientation) where C <: Complex
    for i = 1:length(zs)
        zs[i] = 0.0
    end

    R = √(x^2 + y^2)
    R_d = orientation.dipole_axis[1]*x + orientation.dipole_axis[2]*y
    ρ = R^2 - R_d^2
    # P0 = [R^6, 0, 3R^4, 0, 3R^2, 0, 1]
    z_d = orientation.dipole_axis[3]
    # Pmag = [ρ^2, -4ρ*R_d*z_d, 4R_d^2*z_d^2 + 2ρ - 2ρ*z_d^2, 
    # 4R_d*z_d^3 - 4R_d*z_d, 1 + z_d^4 - 2z_d^2, 0, 0]
    # P_in = P0 - magnetosphere.r_mi^2*Pmag
    poly[1] = R^6 - R_m^2*ρ^2
    poly[2] = 4ρ*R_d*z_d * R_m^2
    poly[3] = 3R^4 - (4R_d^2*z_d^2 + 2ρ - 2ρ*z_d^2) * R_m^2
    poly[4] = -(4R_d*z_d^3 - 4R_d*z_d) * R_m^2
    poly[5] = 3R^2 - (1 + z_d^4 - 2z_d^2) * R_m^2
    poly[6] = 0; poly[7] = 1
    PolynomialRoots.roots!(roots_arr, poly, 1e-10, 6, false)
    n_z = 0
    for i = 1:6
        if abs(imag(roots_arr[i])) < 1e-8
            z = real(roots_arr[i])
            if (x^2 + y^2 ≥ 1.0) | ((z > 0.0) & (z^2 ≥ 1 - x^2 - y^2)) 
                zs[n_z+1] = z
                n_z += 1
            end
        end
    end

    sort!(view(zs, 1:n_z))
    if 0 < n_z < length(zs)
        zs[n_z + 1] = zs[n_z] - 1.0
    end
end

function is_self_absorbed_single_line!(zs, roots_arr, poly, star, geometry, orientation, x, y, z, v_z, Δv_z)
    R_m = (geometry.r_mi + geometry.r_mo)/2
    calc_dipole_line_intersections!(zs, roots_arr, poly, x, y, R_m, orientation)

    if sum(abs.(zs)) < 1e-6
        return false
    end

    n_z = 1
    while zs[n_z + 1] ≥ zs[n_z]
        n_z += 1
    end

    for i_z = 1:n_z
        if zs[i_z] > z + 1e-6
            v_z_cur = -calc_v_z(x, y, zs[i_z], orientation, geometry, star)
            if abs(v_z - v_z_cur) ≤ Δv_z
                return true
            end
        end
    end
    return false
end

function calcborders_prealloc!(borders :: Vector{Float64}, roots_arr :: Vector{C}, poly :: Vector{C}, x :: Real, y :: Real, 
                                magnetosphere :: DipoleGeometry, orientation :: Orientation) where C <: Complex
    for i = 1:length(borders)
        borders[i] = 0.0
    end
    ψ = orientation.ψ
    R = √(x^2 + y^2)
    R_d = orientation.dipole_axis[1]*x + orientation.dipole_axis[2]*y
    ρ = R^2 - R_d^2
    r_mi = magnetosphere.r_mi
    r_mo = magnetosphere.r_mo
    # P0 = [R^6, 0, 3R^4, 0, 3R^2, 0, 1]
    z_d = orientation.dipole_axis[3]
    # Pmag = [ρ^2, -4ρ*R_d*z_d, 4R_d^2*z_d^2 + 2ρ - 2ρ*z_d^2, 
    # 4R_d*z_d^3 - 4R_d*z_d, 1 + z_d^4 - 2z_d^2, 0, 0]
    # P_in = P0 - magnetosphere.r_mi^2*Pmag
    poly[1] = R^6 - r_mi^2*ρ^2
    poly[2] = 4ρ*R_d*z_d * r_mi^2
    poly[3] = 3R^4 - (4R_d^2*z_d^2 + 2ρ - 2ρ*z_d^2) * r_mi^2
    poly[4] = -(4R_d*z_d^3 - 4R_d*z_d) * r_mi^2
    poly[5] = 3R^2 - (1 + z_d^4 - 2z_d^2) * r_mi^2
    poly[6] = 0; poly[7] = 1
    PolynomialRoots.roots!(roots_arr, poly, 1e-10, 6, false)
    n_b = 0
    for i = 1:6
        if abs(imag(roots_arr[i])) < 1e-8
            z = real(roots_arr[i])
            if (x^2 + y^2 ≥ 1.0) | ((z > 0.0) & (z^2 ≥ 1 - x^2 - y^2)) 
                borders[n_b+1] = z
                n_b += 1
            end
        end
    end

    # P_out = P0 - magnetosphere.r_mo^2*Pmag
    poly[1] = R^6 - r_mo^2*ρ^2
    poly[2] = 4ρ*R_d*z_d * r_mo^2
    poly[3] = 3R^4 - (4R_d^2*z_d^2 + 2ρ - 2ρ*z_d^2) * r_mo^2
    poly[4] = -(4R_d*z_d^3 - 4R_d*z_d) * r_mo^2
    poly[5] = 3R^2 - (1 + z_d^4 - 2z_d^2) * r_mo^2
    poly[6] = 0; poly[7] = 1
    PolynomialRoots.roots!(roots_arr, poly, 1e-10, 6, false)
    for i = 1:6
        if abs(imag(roots_arr[i])) < 1e-8 
            z = real(roots_arr'[i])
            if (x^2 + y^2 ≥ 1.0) | ((z > 0.0) & (z^2 ≥ 1 - x^2 - y^2)) 
                borders[n_b+1] = z
                n_b += 1
            end
        end
    end

    if n_b % 2 == 1
        if R > 1
            borders .= 0.0
        else
            borders[n_b + 1] = √(1.0 + 1e-10 - x^2 - y^2)
            n_b += 1
        end
    end

    sort!(view(borders, 1:n_b))
end

function is_self_absorbed(borders :: Vector{Float64}, roots :: Vector{C}, poly :: Vector{C}, 
                                star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, x, y, z, v_z_0, Δv_z, z_step) where C <: Complex
    # kin = Models.SolidDipoleKin(star, 0.0)

    z_borders = calcborders_prealloc!(borders, roots, poly, x, y, geometry, orientation)
    n_skip = 0

    n_borders = 0

    for i_b = 1:2:length(borders)
        if abs(borders[i_b] - borders[i_b+1]) > 1e-10
            n_borders = i_b + 1
        end
    end
    # println(z_borders, " ", z)

    if n_borders == 0
        return false
    end

    for i_z = 1:1:n_borders
        if (z - z_borders[i_z]) > -z_step/10
            n_skip += 1
        else
            break
        end
    end

    if n_skip ≥ n_borders
        return false
    end

    if n_skip % 2 == 0
        z = z_borders[n_skip + 1]
        n_skip += 2
    else
        n_skip += 1

        z_in = z
        z_out = z_borders[n_skip]

        n_z = ceil(Int, (z_out - z_in)/z_step)
        z_step_rounded = (z_out - z_in)/n_z

        z = z + z_step_rounded/2
    end

    for i_border = n_skip:2:n_borders
        z_in = z
        z_out = z_borders[i_border]
        # println(z_in, " ", z_out)
        n_z = ceil(Int, (z_out - z_in)/z_step)
        z_step_rounded = (z_out - z_in)/n_z
        v_z_2 = -calc_v_z(x, y, z_out, orientation, geometry, star)
        v_z_1 = v_z_2

        for i_z = n_z-1:-1:1
            z = z_in + i_z*z_step_rounded
            v_z = -calc_v_z(x, y, z, orientation, geometry, star)
            # println(z_in, " ", z_out, " ", z, " ", v_z_0, " ", v_z)
            if abs(v_z - v_z_0) < Δv_z
                return true
            end
        end
        if i_border < n_borders
            z = z_borders[i_border + 1]
        end
    end

    return false
end

function is_self_absorbed(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, x, y, z, v_z_0, Δv_z, z_step)
    borders = zeros(12)
    roots = zeros(ComplexF64, 6)
    poly = zeros(ComplexF64, 7)
    is_self_absorbed(borders, roots, poly, star, geometry, orientation, x, y, z, v_z_0, Δv_z, z_step)
end

function calc_kernel(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, Δv_z, ζ, v_z)
    R_m = (geometry.r_mi + geometry.r_mo)/2
    ϕs = find_v_z(star, geometry, orientation, R_m, ζ, v_z)
    kernel = 0.0

    zs = zeros(12)
    roots_arr = zeros(ComplexF64, 6)
    poly = zeros(ComplexF64, 7)

    for ϕ in ϕs
        if abs(ϕ) > 3π
            continue
        end
        decart_coordinates = orientation.w_from_d * dipole_decart_coordinates(R_m, ζ, ϕ)
        x, y, z = decart_coordinates
        if (z > 0) | ((x^2 + y^2) ≥ 1)
            self_abs = is_self_absorbed_single_line!(zs, roots_arr, poly, star, geometry, orientation, x, y, z, v_z, Δv_z)
            if !self_abs
                kernel += abs(det(calc_picture_plane_jacobian(star, geometry, orientation, R_m, ζ, ϕ)))
            end
        end
    end
    return kernel
end

function calc_velocity_kernel(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, Δv_z, ζ, v_z)
    R_m = (geometry.r_mi + geometry.r_mo)/2
    ϕs = find_v_z(star, geometry, orientation, R_m, ζ, v_z)
    kernel = 0.0

    zs = zeros(12)
    roots_arr = zeros(ComplexF64, 6)
    poly = zeros(ComplexF64, 7)

    for ϕ in ϕs
        if abs(ϕ) > 3π
            continue
        end
        # borders .= 0.0
        decart_coordinates = orientation.w_from_d * dipole_decart_coordinates(R_m, ζ, ϕ)
        x, y, z = decart_coordinates
        if (z > 0) | ((x^2 + y^2) ≥ 1)
            self_abs = is_self_absorbed_single_line!(zs, roots_arr, poly, star, geometry, orientation, x, y, z, v_z, Δv_z)
            if !self_abs
                kernel += abs(det(calc_picture_plane_velocity_jacobian(star, geometry, orientation, R_m, ζ, ϕ)))
            end
        end
    end
    return kernel
end

function calc_kernel_advanced(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation,
                                            ζ, v_z0, Δv_z; n_Rm = 10, n_div = 100)
    kernel = 0.0
    dRm = (geometry.r_mo - geometry.r_mi)/n_Rm

    borders = zeros(12)
    roots = zeros(ComplexF64, 6)
    poly = zeros(ComplexF64, 7)

    M = orientation.w_from_d

    v_esc = starvescincms(star)
    v_eq = starveqincms(star)

    for i_Rm=1:n_Rm
        R_m = geometry.r_mi + dRm*i_Rm - dRm/2
        for surf_sign in -1:2:1
            v_z = v_z0 + surf_sign*Δv_z
            ϕ_1,ϕ_2 = find_v_z(star, geometry, orientation, R_m, ζ, v_z)
            if abs(ϕ_1) > 2π
                continue
            end
            ϕs = if abs(ϕ_1 - ϕ_2) < 1e-3
                (ϕ_1)
            else
                (ϕ_1, ϕ_2)
            end
            for ϕ in ϕs
                J_d = gridjacobianmatrix(geometry, R_m, ζ, ϕ)
                inv_J = workgridinversejacobianmatrix(geometry, R_m, ζ, ϕ, orientation)
                x, y, z = orientation.w_from_d * dipole_decart_coordinates(R_m, ζ, ϕ)
                if (z > 0) | ((x^2 + y^2) ≥ 1)
                    dvz_d = calc_radial_velocity_gradient(v_esc, v_eq, R_m, ζ, ϕ, J_d, M)
                    dvz_dz = calc_dipole_radial_velocity_los_derivative(inv_J, dvz_d)
                    if dvz_dz*surf_sign > 0.0
                        self_abs = !is_self_absorbed(borders, roots, poly, star, geometry, orientation, x, y, z, v_z0, Δv_z, 0.01)
                        if self_abs
                            kernel += abs(det(calc_picture_plane_jacobian(J_d, M, dvz_d)))*dRm
                        end
                    end
                end
            end
            # print("; ")
        end
    end
    return kernel
end

function calc_velocity_kernel_advanced(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation,
                        ζ, v_z0, Δv_z; n_vz = 10, n_div = 100)
    kernel = 0.0
    dvz = 2Δv_z/n_vz
    R_in = geometry.r_mi
    R_out = geometry.r_mo
    R_mid = (R_in + R_out)/2
    W = R_out - R_in

    borders = zeros(12)
    roots = zeros(ComplexF64, 6)
    poly = zeros(ComplexF64, 7)

    M = orientation.w_from_d

    v_esc = starvescincms(star)
    v_eq = starveqincms(star)

    for i_vz in 1:n_vz
        v_z = v_z0 - Δv_z + dvz*i_vz - dvz/2
        for surf_sign in -1:2:1
            R_m = R_mid + surf_sign*W/2
            ϕ_1,ϕ_2 = find_v_z(star, geometry, orientation, R_m, ζ, v_z)
            if abs(ϕ_1) > 2π
                continue
            end
            ϕs = if abs(ϕ_1 - ϕ_2) < 1e-3
                (ϕ_1)
            else
                (ϕ_1, ϕ_2)
            end
            for ϕ in ϕs
                J_d = gridjacobianmatrix(geometry, R_m, ζ, ϕ)
                x, y, z = orientation.w_from_d * dipole_decart_coordinates(R_m, ζ, ϕ)
                if (z > 0) | ((x^2 + y^2) ≥ 1)
                    dvz_d = calc_radial_velocity_gradient(v_esc, v_eq, R_m, ζ, ϕ, J_d, M)     
                    dRm_dz = calc_dipole_magnetosphere_radius_los_derivative(R_m, ζ, ϕ, orientation)
                    if dRm_dz*surf_sign > 0.0
                        self_abs = !is_self_absorbed(borders, roots, poly, star, geometry, orientation, x, y, z, v_z0, Δv_z, 0.01)
                        if self_abs
                            kernel += abs(det(calc_picture_plane_velocity_jacobian(J_d, M, dvz_d)))*dvz
                        end
                    end
                end
            end
        end
    end 
    return kernel
end

function calc_absorption_profile(x_arr, y_arr, dS_arr, star :: AbstractStar, geometry :: DipoleGeometry, 
                                            orientation :: Orientation, v_z_arr, Δv_z, grid_step, z_step, hotspot_val)
    n_grid = length(x_arr)
    absorption = 0.0

    n_v_z = length(v_z_arr)
    absorption_profile = zeros(n_v_z)
    absorption_profile_ray = zeros(n_v_z)

    θ_1 = asin(√(1/geometry.r_mo))
    θ_2 = asin(√(1/geometry.r_mi))

    hotspot_star = MagnetosphereSpotStar(star, 1e4, θ_1, θ_2)

    borders = zeros(12)
    roots = zeros(ComplexF64, 6)
    poly = zeros(ComplexF64, 7)

    for i_grid = 1:n_grid
        # println("$x, $y, $n_borders, $dS, $(√(x^2 + y^2))")
        x = x_arr[i_grid]; y = y_arr[i_grid]; dS = dS_arr[i_grid] 

        calcborders_prealloc!(borders, roots, poly, x, y, geometry, orientation)
        
        n_borders = 0
        for i_b = 1:2:length(borders)
            if abs(borders[i_b] - borders[i_b+1]) > 1e-10
                n_borders = i_b + 1
            end
        end

        if n_borders == 0; continue; end
        absorption_profile_ray .= 0.0
        # println("$i_grid $n_grid $x, $y, $n_borders, $dS, $(√(x^2 + y^2))")
        for i_out = n_borders:-2:2
            z_in = borders[i_out-1]; z_out = borders[i_out]
            # println("\t$z_in, $z_out")
            n_z = ceil(Int, (z_out - z_in)/z_step)
            z_step_rounded = (z_out - z_in)/n_z
            v_z_2 = -calc_v_z(x, y, z_out, orientation, geometry, star)
            v_z_1 = v_z_2
            for i_z = n_z-1:-1:1
                z = z_in + i_z*z_step_rounded
                v_z = -calc_v_z(x, y, z, orientation, geometry, star)
                if v_z ≤ v_z_1
                    v_z_1 = v_z
                elseif v_z ≥ v_z_2
                    v_z_2 = v_z
                end
            end
            v_z_1 = v_z_1 - Δv_z
            v_z_2 = v_z_2 + Δv_z
            # println("\t$v_z_1, $v_z_2")
            for i_v_z = 1:n_v_z
                v_z = v_z_arr[i_v_z]
                if (v_z_1 ≤ v_z ≤ v_z_2) & (absorption_profile_ray[i_v_z] < dS)
                    absorption_profile_ray[i_v_z] = dS                
                end
            end
        end
        absorption_profile += absorption_profile_ray
    end
    absorption_profile
end

function calc_absorption_profile(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, v_z_arr, Δv_z, grid_step, z_step, hotspot_val = 1.0)
    x_arr, y_arr, dS_arr = polargrid(1.0, grid_step)

    θ_1 = asin(√(1/geometry.r_mo))
    θ_2 = asin(√(1/geometry.r_mi))

    hotspot_star = MagnetosphereSpotStar(star, 1e4, θ_1, θ_2)
    n_grid = length(x_arr)
    for i_grid = 1:n_grid
        # println("$x, $y, $n_borders, $dS, $(√(x^2 + y^2))")
        x = x_arr[i_grid]; y = y_arr[i_grid]; dS = dS_arr[i_grid] 
        if abs(hotspot_val - 1) > 1e-8
            if isthispointonspot(hotspot_star, x, y, orientation.star_axis)
                dS_arr[i_grid] *= hotspot_val
            end
        end
    end

    calc_absorption_profile(x_arr, y_arr, dS_arr, star, geometry, orientation, v_z_arr, Δv_z, grid_step, z_step, hotspot_val) / sum(dS_arr)
end

function calc_absorption_profile_parallel(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, v_z_arr, Δv_z, grid_step, z_step, hotspot_val = 1.0)
    x_arr, y_arr, dS_arr = polargrid(1.0, grid_step)

    θ_1 = asin(√(1/geometry.r_mo))
    θ_2 = asin(√(1/geometry.r_mi))

    hotspot_star = MagnetosphereSpotStar(star, 1e4, θ_1, θ_2)
    n_grid = length(x_arr)
    for i_grid = 1:n_grid
        # println("$x, $y, $n_borders, $dS, $(√(x^2 + y^2))")
        x = x_arr[i_grid]; y = y_arr[i_grid]; dS = dS_arr[i_grid] 
        if abs(hotspot_val - 1) > 1e-8
            if isthispointonspot(hotspot_star, x, y, orientation.star_axis)
                dS_arr[i_grid] *= hotspot_val
            end
        end
    end

    abs_prof = zeros(length(v_z_arr))

    n_threads = Threads.nthreads(:default)

    n_jobs = 2n_threads

    jobs = Channel{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}(n_jobs)
    results = Channel{Tuple{Int, Vector{Float64}}}(n_jobs)

    println("Jobs creation")

    for i_job = 1:n_jobs
        put!(jobs, (x_arr[i_job:n_jobs:end], y_arr[i_job:n_jobs:end], dS_arr[i_job:n_jobs:end]))
    end

    println("Work tasks, $n_threads")

    workers = []

    for i_thread = 1:n_threads
        w = Threads.@spawn for job in jobs
            x_arr_job, y_arr_job, dS_arr_job = job
            # println(length(x_arr_job))
            abs_prof_job = calc_absorption_profile(x_arr_job, y_arr_job, dS_arr_job, star, geometry, orientation, v_z_arr, 
                                                            Δv_z, grid_step, z_step, hotspot_val)
            # println(threadid())
            put!(results, (Threads.threadid(), abs_prof_job))
            # println("put")
        end
        push!(workers, w)
    end

    print("Nothing is done. ")
    n_done = 0
    progress_worker = Threads.@spawn :interactive  while n_jobs > 0
        print("Waiting...")
        i_thread, abs_prof_job = take!(results)
        n_done += 1
        abs_prof += abs_prof_job
        n_jobs -= 1
        print("\e[2K\e[1GTask $n_done is done (thread $i_thread). $n_jobs remaining tasks. ")
    end
    wait(progress_worker)
    print("\n")

    close(jobs)
    close(results)

    abs_prof / sum(dS_arr)
end

function calc_absorption_profile(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, n_freq :: Int, 
                                                v_z_borders, Δv_z, grid_step, z_step)
    v_z_start, v_z_end = v_z_borders
    v_z_step = (v_z_end - v_z_start)/n_freq                                            
    v_z_arr = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n_freq]
    calc_absorption_profile(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, v_z_arr, Δv_z, grid_step, z_step)
end

function correct_absorption(star :: AbstractStar, geometry :: DipoleGeometry, orientation :: Orientation, grid_step, absorption, hotspot_val)
    θ_1 = asin(√(1/geometry.r_mo))
    θ_2 = asin(√(1/geometry.r_mi))

    hotspot_star = MagnetosphereSpotStar(star, 1e4, θ_1, θ_2)

    x_arr, y_arr, dS_arr = polargrid(1.0, grid_step)

    sumS = sum(dS_arr)

    n_grid = length(x_arr)
    for i_grid = 1:n_grid
        # println("$x, $y, $n_borders, $dS, $(√(x^2 + y^2))")
        x = x_arr[i_grid]; y = y_arr[i_grid]; dS = dS_arr[i_grid] 
        if abs(hotspot_val - 1) > 1e-8
            if isthispointonspot(hotspot_star, x, y, orientation.star_axis)
                dS_arr[i_grid]  *= hotspot_val
            end
        end
    end

    absorption_correct =  absorption * sumS/sum(dS_arr)
    return absorption_correct
end