using CairoMakie
using LaTeXStrings

include("../include.jl")

CairoMakie.activate!()

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

param_symb = [:inc, :R_in, :W, :Δv_z, :H]
param_id = Dict(:inc => 1, :R_in => 2, :W => 3, :Δv_z => 4, :H => 5)
param_latex = Dict(:inc => "i,\\ ^\\circ", :R_in => "R_\\mathrm{in},\\ \\mathrm{R_\\star}", :W => "W,\\ \\mathrm{R_\\star}", :Δv_z => "\\Delta v_z,\\ \\text{cm/s}", :H => "H")

function plot_residuals_pds70(residuals_file; kernels_dir = "../kernels", residuals_dir = "../residuals")
    parameter_grid = load_kernel_models_grid(; kernels_dir)
    

    fig = Figure()
    axes = []

    residuals = open("$residuals_dir/$residuals_file", "r") do io
        read_array(io, Float64)
    end

    cb_range = (-1.5,-1.2)#(log10.(minimum(residuals)), log10.(√2*minimum(residuals)))

    for (x_param_id, y_param_id, ax_id) in tuple.([1,3,3,5], [2,2,4,4], 1:4) 

        ax = Axis(fig[(ax_id-1) ÷ 2 + 1, (ax_id-1) % 2 + 1][1,1],
                    xlabel = latexstring(param_latex[param_symb[x_param_id]]), ylabel = latexstring(param_latex[param_symb[y_param_id]]))

        x_param_grid = parameter_grid[x_param_id]
        y_param_grid = parameter_grid[y_param_id]

        hm_data = zeros(length(x_param_grid), length(y_param_grid))
        for i_x in eachindex(x_param_grid), i_y in eachindex(y_param_grid)
            if y_param_id > x_param_id
                data = selectdim(selectdim(residuals, y_param_id, i_y), x_param_id, i_x)
                hm_data[i_x, i_y] = minimum(data)
            elseif x_param_id > y_param_id
                data = selectdim(selectdim(residuals, x_param_id, i_x), y_param_id, i_y)
                hm_data[i_x, i_y] = minimum(data)
            else
                return Figure()
            end
        end

        hm = heatmap!(ax, x_param_grid, y_param_grid, log10.(hm_data),  colormap = :grayC)#, colorrange = cb_range)
        ct = contour!(ax, x_param_grid, y_param_grid, log10.(hm_data), levels = [-2:0.05:-1;], color = :black)
        Colorbar(fig[(ax_id-1) ÷ 2 + 1, (ax_id-1) % 2 + 1][1,2], hm, label = L"\lg \delta")
    end

    fig
end


function plot_residuals(line, params)
    parameter_grid = load_kernel_models_grid(kernels_dir = "../kernels_bin")
    x_param_grid = parameter_grid[param_id[params[1]]]
    y_param_grid = parameter_grid[param_id[params[2]]]
    incs = [15:20:75;]
    fig = Figure()
    axes = []

    for i_inc in eachindex(incs)
        ax = Axis(fig[(i_inc-1)÷2 + 1, (i_inc-1)%2 + 1][1,1], title=latexstring("i_\\mathrm{mod} = ", incs[i_inc], "^\\circ"),
                    xlabel = latexstring(param_latex[params[1]]), ylabel = latexstring(param_latex[params[2]]))
        push!(axes, ax)
    end

    u, l = line
    hm = 1
    for i_inc in eachindex(incs)
        inc = incs[i_inc]
        param_hart = Dict(:inc => inc, :R_in => 3.0, :W => 1.0, :Δv_z => NaN)
        residuals = open("../residuals_H/residuals_3-4_$(inc)_$u$l.bin", "r") do io
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
        hm = heatmap!(axes[i_inc], x_param_grid, y_param_grid, log10.(hm_data), colormap = :grayC)
        Colorbar(fig[(i_inc-1)÷2 + 1, (i_inc-1)%2 + 1][1,2], hm, label = L"\lg \delta")

        

        # min_id = findmin(hm_data)[2]
        min_x = x_param_grid[min_id[param_id[params[1]]]]
        min_y = y_param_grid[min_id[param_id[params[2]]]]        
        scatter!(axes[i_inc], [param_hart[params[1]]], [param_hart[params[2]]], marker = :circle, markersize = 12, color = :white, strokewidth = 1; strokecolor = :black)
        scatter!(axes[i_inc], [min_x], [min_y], marker = :xcross, color = :white, strokewidth = 1; strokecolor = :black)
    end
    return fig
end

function plot_source(line, i_H)
    incs = [15:20:75;]
    n_inc = length(incs)

    fig = Figure()

    source_axes = []
    prof_axes = []
    u,l = line

    # grid_incs, R_ins, Ws, Δv_zs, Hs = load_kernel_models_grid(; kernels_dir = "../kernels_abs")
    v_z_borders, n_freq, n_ζ = load_kernel_grid(; kernels_dir = "../kernels_bin")

    for (i_inc, inc) in enumerate(incs)
        name = "3-4_$(inc)_$u$l"
        residuals = open("../residuals_H/residuals_3-4_$(inc)_$u$l.bin", "r") do io
            read_array(io, Float64)
        end

        io = open("../hart_to_fit/mag_$(name)_P.bin", "r") 
        profile_v_zs = reverse(read_array(io, Float64))
        profile = reverse(read_array(io, Float64))
        profile[1] = 1.0; profile[end] = 1.0
        close(io)

        v_z_start, v_z_end = v_z_borders
        v_z_step = (v_z_end - v_z_start)/n_freq
        v_zs = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n_freq]

        profile_model_grid = fill(1.0, n_freq)

        profile_int = LinearInterpolation(profile, profile_v_zs)
        for i_vz in 1:n_freq
            v_z = v_zs[i_vz]
            if minimum(profile_v_zs) < v_z < maximum(profile_v_zs)
                profile_model_grid[i_vz] = profile_int(v_z)
            end
        end

        mod_i_inc, i_Rin, i_W, i_Δvz, _ = Tuple(findmin(residuals)[2])
        kernel_matrix = kernel_ζ_n_half(load_kernel_from_index(mod_i_inc, i_Rin, i_W, i_Δvz; kernels_dir = "../kernels_bin"), 3)
        absorption = load_absorption_from_index(mod_i_inc, i_Rin, i_W, i_Δvz, i_H; kernels_dir = "../kernels_abs")
        emission_profile = profile_model_grid + absorption
        lg_restored_source = fit_emission_profile(emission_profile, kernel_matrix, v_z_borders, 2.0, 1.0, 0.0)
        restored_emission = (kernel_matrix / (n_ζ/8)) * (10 .^ lg_restored_source)
        restored_profile = restored_emission - absorption .+ 1.0

        source_known = open("../hart_to_fit/mag_$(name)_S.bin", "r") do io
            read_array(io, Float64)
        end

        source_ax = Axis(fig[i_inc, 1], xlabel = L"\zeta", ylabel = L"\lg(S/<I_\mathrm{cont}>)", xtickalign = 1, xticksmirrored = true)
        source_ax.xticks = [0:0.2:1;]

        R_ms = [3,3.5,4]
        n_S = size(source_known)[2]
        for i_Rm = 1:3
            R_m = R_ms[i_Rm]
            rs = collect(range(1,R_m,n_S))
            θs = asin.(sqrt.(rs/R_m))
            ζs = zeros(n_S)
            for i_ζ = 1:n_S
                r_m, ζ = spheretogridproper(DipoleGeometry(3,4), rs[i_ζ], θs[i_ζ])
                ζs[i_ζ] = ζ
            end
            lines!(source_ax, ζs, log10.(source_known[i_Rm,:]), linestyle = :dot, color = :black, label = "Model")
        end

        ζs = [i_ζ/n_ζ*8 - 1/2n_ζ*8 for i_ζ = 1:(n_ζ/8)]
        lines!(source_ax, ζs, lg_restored_source, color = :black, label = "Restored")

        profile_ax = Axis(fig[i_inc, 2], xlabel = L"v_z,\ \text{cm/s}", ylabel = L"F_\lambda/F_\star", xtickalign = 1, xticksmirrored = true)

        

        lines!(profile_ax, profile_v_zs, profile, linestyle = :dot, color = :black, label = "Model")
        lines!(profile_ax, v_zs, restored_profile, color = :black, label = "Optimal")
        lines!(profile_ax, v_zs, 1.0 .+ restored_emission, linestyle = :dash, color = :black, label = L"E_\lambda,\ A_\lambda")
        lines!(profile_ax, v_zs, 1.0 .- absorption, linestyle = :dash, color = :black)
        push!(source_axes, source_ax)
        push!(prof_axes, profile_ax)
        text!(source_ax, latexstring("i_\\mathrm{mod} = $inc^\\circ,\\ $(series[l])$(greek[u-l])"), space = :relative, position = Point2f(0.1, 0.9), align = (:left, :top))
    end
    
    for i_inc = 1:(n_inc-1)
        hidexdecorations!(prof_axes[i_inc], ticks = false, grid = false)
        hidexdecorations!(source_axes[i_inc], ticks = false, grid = false)
    end

    axislegend(prof_axes[1])


    rowgap!(fig.layout, 0)

    # global_title = Label(fig[0, :], latexstring(series[l], greek[u-l]), 
    #                  fontsize = 20, 
    #                  font = :bold)

    fig
end

function plot_pabr_source(i_inc, i_H_43, i_H_53, i_H_74)
    incs = [15:20:75;]

    inc = incs[i_inc]

    i_H_arr = [i_H_43, i_H_53, i_H_74]
    line_arr = [(4,3), (5,3), (7,4)]

    fig = Figure()
    source_axes = []
    prof_axes = []

    v_z_borders, n_freq, n_ζ = load_kernel_grid(; kernels_dir = "../kernels_bin")
    n_line = 3

    for (i_line, line) in enumerate(line_arr)
        u,l = line
        i_H = i_H_arr[i_line]

        name = "3-4_$(inc)_$u$l"
        residuals = open("../residuals_H/residuals_3-4_$(inc)_$u$l.bin", "r") do io
            read_array(io, Float64)
        end

        io = open("../hart_to_fit/mag_$(name)_P.bin", "r") 
        profile_v_zs = reverse(read_array(io, Float64))
        profile = reverse(read_array(io, Float64))
        profile[1] = 1.0; profile[end] = 1.0
        close(io)

        v_z_start, v_z_end = v_z_borders
        v_z_step = (v_z_end - v_z_start)/n_freq
        v_zs = [v_z_start + v_z_step*i_v_z - v_z_step/2 for i_v_z = 1:n_freq]

        profile_model_grid = fill(1.0, n_freq)

        profile_int = LinearInterpolation(profile, profile_v_zs)
        for i_vz in 1:n_freq
            v_z = v_zs[i_vz]
            if minimum(profile_v_zs) < v_z < maximum(profile_v_zs)
                profile_model_grid[i_vz] = profile_int(v_z)
            end
        end

        mod_i_inc, i_Rin, i_W, i_Δvz, _ = Tuple(findmin(residuals)[2])
        kernel_matrix = kernel_ζ_n_half(load_kernel_from_index(mod_i_inc, i_Rin, i_W, i_Δvz; kernels_dir = "../kernels_bin"), 3)
        absorption = load_absorption_from_index(mod_i_inc, i_Rin, i_W, i_Δvz, i_H; kernels_dir = "../kernels_abs")
        emission_profile = profile_model_grid + absorption
        lg_restored_source = fit_emission_profile(emission_profile, kernel_matrix, v_z_borders, 2.0, 1.0, 0.0)
        restored_emission = (kernel_matrix / (n_ζ/8)) * (10 .^ lg_restored_source)
        restored_profile = restored_emission - absorption .+ 1.0

        source_known = open("../hart_to_fit/mag_$(name)_S.bin", "r") do io
            read_array(io, Float64)
        end

        source_ax = Axis(fig[i_line, 1], xlabel = L"\zeta", ylabel = L"\lg(S/<I_\mathrm{cont}>)", xtickalign = 1, xticksmirrored = true)
        source_ax.xticks = [0:0.2:1;]

        R_ms = [3,3.5,4]
        n_S = size(source_known)[2]
        for i_Rm = 1:3
            R_m = R_ms[i_Rm]
            rs = collect(range(1,R_m,n_S))
            θs = asin.(sqrt.(rs/R_m))
            ζs = zeros(n_S)
            for i_ζ = 1:n_S
                r_m, ζ = spheretogridproper(DipoleGeometry(3,4), rs[i_ζ], θs[i_ζ])
                ζs[i_ζ] = ζ
            end
            lines!(source_ax, ζs, log10.(source_known[i_Rm,:]), linestyle = :dot, color = :black, label = "Model")
        end

        ζs = [i_ζ/n_ζ*8 - 1/2n_ζ*8 for i_ζ = 1:(n_ζ/8)]
        lines!(source_ax, ζs, lg_restored_source, color = :black, label = "Restored")

        profile_ax = Axis(fig[i_line, 2], xlabel = L"v_z,\ \text{cm/s}", ylabel = L"F_\lambda/F_\star", xtickalign = 1, xticksmirrored = true)

        

        lines!(profile_ax, profile_v_zs, profile, linestyle = :dot, color = :black, label = "Model")
        lines!(profile_ax, v_zs, restored_profile, color = :black, label = "Optimal")
        lines!(profile_ax, v_zs, 1.0 .+ restored_emission, linestyle = :dash, color = :black, label = L"E_\lambda,\ A_\lambda")
        lines!(profile_ax, v_zs, 1.0 .- absorption, linestyle = :dash, color = :black)
        push!(source_axes, source_ax)
        push!(prof_axes, profile_ax)
        text!(source_ax, latexstring("i_\\mathrm{mod} = $inc^\\circ,\\ $(series[l])$(greek[u-l])"), space = :relative, position = Point2f(0.1, 0.9), align = (:left, :top))
    end

    for i_line = 1:(n_line-1)
        hidexdecorations!(prof_axes[i_line], ticks = false, grid = false)
        hidexdecorations!(source_axes[i_line], ticks = false, grid = false)
    end

    axislegend(prof_axes[1])


    rowgap!(fig.layout, 0)

    # global_title = Label(fig[0, :], latexstring(series[l], greek[u-l]), 
    #                  fontsize = 20, 
    #                  font = :bold)

    fig
end