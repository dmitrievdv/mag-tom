struct Orientation
    i :: Float64
    α :: Float64
    ϕ :: Float64
    ζ :: Float64
    ψ :: Float64
    # w -- Work; s -- Star; d -- Dipole 
    d_from_w :: SArray{Tuple{3,3},Float64,2,9}
    w_from_d :: SArray{Tuple{3,3},Float64,2,9}
    s_from_w :: SArray{Tuple{3,3},Float64,2,9}
    w_from_s :: SArray{Tuple{3,3},Float64,2,9}
    d_from_s :: SArray{Tuple{3,3},Float64,2,9}
    s_from_d :: SArray{Tuple{3,3},Float64,2,9}
    dipole_axis :: SArray{Tuple{3},Float64,1,3}
    star_axis :: SArray{Tuple{3},Float64,1,3}

    function Orientation(i_deg, α_deg, ϕ_deg)
        i = i_deg/180*π
        α = α_deg/180*π
        ϕ = ϕ_deg/180*π
        R_α = SA[1    0      0;
                    0  cos(α) sin(α); 
                    0 -sin(α) cos(α)]
        R_ϕ = SA[cos(ϕ) sin(ϕ) 0;
                -sin(ϕ) cos(ϕ) 0; 
                  0       0    1]
        R_i =    SA[1    0      0;
                    0  cos(i) sin(i); 
                    0 -sin(i) cos(i)]
        
        w_from_d = R_i*R_ϕ'*R_α' # R_α is inverted, this way at ϕ = 0 dipole axis is at the closest to line of sight
        d_from_w = R_α*R_ϕ*R_i'
        w_from_s = R_i*R_ϕ'
        s_from_w = R_ϕ*R_i'
        d_from_s = R_α
        s_from_d = R_α'

        dipole_axis = w_from_d*SA[0, 0, 1]
        star_axis = R_i*[0,0,1]
        ψ = acos(dipole_axis[3])
        ζ = 0
        return new(i, α, ϕ, ζ, ψ, d_from_w, w_from_d, s_from_w, w_from_s, d_from_s, s_from_d, dipole_axis, star_axis)
    end

    function Orientation(i_deg)
        Orientation(i_deg, 0.0, 0.0)
    end
end