function decartmatrixgrid(square_bound, m)
    decartmatrixgrid((-square_bound, square_bound), (-square_bound, square_bound), m)
end

function decartmatrixgrid(x_bounds, y_bounds, m)
    x_min, x_max = x_bounds
    y_min, y_max = y_bounds
    xs = zeros((m, m))
    ys = zeros((m, m))
    dS = zeros((m, m))
    x_step = (x_max - x_min)/m
    y_step = (y_max - y_min)/m

    for i=1:m
        x = x_min + i*x_step - x_step/2
        for j=1:m
            y = y_min + j*y_step - y_step/2
            xs[i,j] = x
            ys[i,j] = y
            dS[i,j] = x_step*y_step
        end
    end
    return xs, ys, dS
end

function decartgrid(xbounds :: Tuple{Real, Real}, ybounds :: Tuple{Real, Real}, m :: Int)
    x_min, x_max = xbounds
    y_min, y_max = ybounds
    xs = zeros(m*m)
    ys = zeros(m*m)
    dS = zeros(m*m)
    x_step = (x_max - x_min)/m
    y_step = (y_max - y_min)/m
    dS .= x_step*y_step
    for i=1:m
        for j=1:m
            k =  (i-1)*m + j
            xs[k] = x_min + i*x_step - x_step/2
            ys[k] = y_min + j*y_step - y_step/2
        end
    end
    return xs, ys, dS
end

function decartgrid(square_bound :: Real, m :: Int)
    decartgrid((-square_bound, square_bound), (-square_bound, square_bound), m)
end

function polargrid(rbound :: Real, r_step :: Float64)
    m = Int(ceil(rbound/r_step))
    polargrid(rbound, m)
end

function polargrid(rbounds :: Tuple{Real, Real}, r_step :: Float64)
    Δr = rbounds[2] - rbounds[1]
    m = Int(ceil(Δr/r_step))
    polargrid(rbounds, m)
end

function polargrid(rbounds :: Tuple{Real, Real}, m :: Int)
    r_min, r_max = rbounds
    r_step = (r_max - r_min)/m
    xs = Float64[]; ys = Float64[]; dS = Float64[]
    for i=1:m  
        r = r_min + r_step*i - r_step/2
        m_ϕ = Int(floor(2π*r/r_step))
        ϕ_step = 2π/m_ϕ
        for j=1:m_ϕ
            ϕ = ϕ_step*j - ϕ_step/(1 + i%2)
            push!(xs, r*cos(ϕ))
            push!(ys, r*sin(ϕ)) 
            push!(dS, r*r_step*ϕ_step)
        end
    end
    return xs, ys, dS
end

function polargrid(rbound :: Real, m :: Int)
    r_step = rbound/m
    xs = Float64[]; ys = Float64[]; dS = Float64[]
    # xs[1] = 0.0; ys[1] = 0.0; dS[1] = π*(r_step*1.5)^2
    for i=1:m  
        r = r_step*i - r_step/2
        m_ϕ = Int(floor(2π*r/r_step))
        ϕ_step = 2π/m_ϕ
        for j=1:m_ϕ
            ϕ = ϕ_step*j - ϕ_step/(1 + i%2)
            push!(xs, r*cos(ϕ))
            push!(ys, r*sin(ϕ)) 
            push!(dS, r*r_step*ϕ_step)
        end
    end
    return xs, ys, dS
end

function simplepolargrid(rbounds :: Tuple{Real, Real}, r_step)
    r_min, r_max = rbounds
    Δr = r_max - r_min
    m_r = Int(ceil(Δr/r_step))
    m_ϕ = Int(ceil(2π*r_min/r_step))
    ϕ_step = 2π/m_ϕ
    xs = Float64[]; ys = Float64[]; dS = Float64[]

    for i=1:m_r
        r = r_min + r_step*i - r_step/2
        for j=1:m_ϕ
            ϕ = ϕ_step*j - ϕ_step/(1 + i%2)
            push!(xs, r*cos(ϕ))
            push!(ys, r*sin(ϕ)) 
            push!(dS, r*r_step*ϕ_step)
        end
    end
    return xs, ys, dS
end

function appendgrid!(grid, grid2)
    append!.(grid, grid2)
end