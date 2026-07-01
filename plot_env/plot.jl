using GLMakie

include("../include.jl")

ker_bin = kernel_vz_n_half(kernel_ζ_n_half(load_kernel_from_index(2,1,1,1, kernels_dir = "../kernels"),2),2)
ker_txt = kernel_vz_n_half(kernel_ζ_n_half(load_kernel_from_index_txt(2,1,1,1,512,1024, kernels_dir = "../kernels_txt"),4),1)

fig = Figure()
ax_bin = Axis(fig[1,1])
ax_txt = Axis(fig[1,2])

color_range = Observable((0,1e3))

color_sl_max = Slider(fig[2,1], range = -5:0.01:4, value = 3)
color_sl_lab = Label(fig[2,2], tellwidth = false, text = @lift(sprint(show, $color_range)))


on(color_sl_max.value) do val
    color_range[] = (0, 10^val)    
    print("\e[2K\e[1G", color_range)
end

linkaxes!(ax_bin, ax_txt)

heatmap!(ax_bin, ker_bin, colorrange = color_range)
heatmap!(ax_txt, ker_txt, colorrange = color_range)

scene = display(fig)
wait(scene)
print("\n")