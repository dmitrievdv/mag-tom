using CairoMakie
using LaTeXStrings

include("../include.jl")

series = Dict(1 => "\\mathrm{La}", 2 => "\\mathrm{H}", 3 => "\\mathrm{Pa}", 4 => "\\mathrm{Br}")
greek = Dict(1 => "\\alpha", 2 => "\\beta", 3 => "\\gamma")

# function line_name(u.)

function plot_hart(lines)
    # lines = [(3,2), (4,2), (5,2), (4,3), (5,3), (7,4)] 
    # lines = [(3,2), (4,2), (5,2)] 
    # lines = [(4,3), (5,3), (7,4)]
    incs = [15:20:75;]
    profile_files = ["mag_3-4_$(inc)_$u$(l)_P.bin" for inc in incs, (u,l) in lines]

    profiles_v_zs = Vector{Float64}[]
    profiles = Vector{Float64}[]

    for profile_file in profile_files
        open("../hart_to_fit/$profile_file", "r") do io
            v_zs = reverse(read_array(io, Float64))
            profile = reverse(read_array(io, Float64))
            push!(profiles_v_zs, v_zs)
            push!(profiles, profile)
        end
    end

    fig = Figure()
    axes = []
    for i_inc = eachindex(incs)
        ax = Axis(fig[(i_inc-1)÷2 + 1, (i_inc-1)%2 + 1], title=latexstring("i = ", incs[i_inc]))
        push!(axes, ax)
    end

    for (i_profile, profile_file) in enumerate(profile_files)
        inc = parse(Float64, profile_file[9:10])
        u = parse(Int64, profile_file[12])
        l = parse(Int64, profile_file[13])

        line_latex = latexstring(series[l], greek[u-l])

        i_inc = findfirst(x -> x==inc, incs)
        lines!(axes[i_inc], profiles_v_zs[i_profile]/1e5, profiles[i_profile], label = line_latex)

    end
    fig[3, 1:2] = Legend(fig, axes[1], orientation = :horizontal, tellheight = true)
    return fig
end

param_id = Dict(:inc => 1, :R_in => 2, :W => 3, :Δv_z => 4, :H => 5)
param_latex = Dict(:inc => "i", :R_in => "R_\\mathrm{in}", :W => "W", :Δv_z => "\\Delta v_z", :H => "H")

function plot_residuals(line, params)
    parameter_grid = load_kernel_models_grid(kernels_dir = "../kernels_bin")
    x_param_grid = parameter_grid[param_id[params[1]]]
    y_param_grid = parameter_grid[param_id[params[2]]]
    incs = [15:20:75;]
    fig = Figure()
    axes = []

    for i_inc in eachindex(incs)
        ax = Axis(fig[(i_inc-1)÷2 + 1, (i_inc-1)%2 + 1], title=latexstring("i = ", incs[i_inc]))
        push!(axes, ax)
    end

    u, l = line

    for i_inc in eachindex(incs)
        inc = incs[i_inc]
        residuals = open("../residuals/residuals_3-4_$(inc)_$u$l.bin", "r") do io
            read_array(io, Float64)
        end
        min_id = findmin(residuals)[2]
        println(min_id)
        println([arr[ind] for (arr, ind) in zip(parameter_grid, Tuple(min_id))])
        hm_data = zeros(length(x_param_grid), length(y_param_grid))
        for i_x in eachindex(x_param_grid), i_y in eachindex(y_param_grid)
            if param_id[params[2]] > param_id[params[1]]
                data = selectdim(selectdim(residuals, param_id[params[2]], i_y), param_id[params[1]], i_x)
                hm_data[i_x, i_y] = minimum(data)
            elseif param_id[params[1]] > param_id[params[2]]
                data = selectdim(selectdim(residuals, param_id[params[1]], i_x), param_id[params[2]], i_y)
                hm_data[i_x, i_y] = minimum(data)
            else
                return Figure()
            end
        end
        heatmap!(axes[i_inc], x_param_grid, y_param_grid, log10.(hm_data))

        # min_id = findmin(hm_data)[2]
        min_x = x_param_grid[min_id[param_id[params[1]]]]
        min_y = y_param_grid[min_id[param_id[params[2]]]]
        scatter!(axes[i_inc], [min_x], [min_y])
    end
    return fig
end