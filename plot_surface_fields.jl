using JLD2
using CairoMakie
using Oceananigans
using Oceananigans.Operators

function check_ranges(folder, ranks; H = 7, fields = false, time = 0, iteration = 0)
    Nx = Vector(undef, ranks)
    iranges = Vector(undef, ranks)
    for rank in 0:ranks - 1
        var = fields ? jldopen(folder * "RealisticOcean_fields_$(rank).jld2")["timeseries/u/"*string(time)][H+1:end-H, H+1:end-H, 1] :
                       jldopen(folder * "RealisticOcean_checkpoint_$(rank)_iteration$(iteration).jld2")["u/data"][H+1:end-H, H+1:end-H, H+1:end-H] 
        Nx[rank+1] = size(var, 1)
    end
    iranges[1] = UnitRange(1, Nx[1])
    for rank in 2:ranks
        iranges[rank] = UnitRange(iranges[rank-1][end]+1,iranges[rank-1][end]+Nx[rank])
    end

    return iranges
end

grid = LatitudeLongitudeGrid(size = (4320, 1800, 1), latitude = (-75, 75), longitude = (-180, 180), z = (0, 1))

u = XFaceField(grid)
v = YFaceField(grid)
T = CenterField(grid)

ζ_op = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)
ζ    = Field(ζ_op)

iterations = parse.(Int, keys(jldopen("RealisticOcean_fields_0.jld2")["timeseries/u"]))
iter = Observable(iterations[1])

iranges = check_ranges("./", 8; fields = true)

ζl = @lift begin
    ud = zeros(size(u))
    vd = zeros(size(v))

    for rank in 0:7
        irange = iranges[rank+1]
        ud[irange, :, :] .= jldopen("RealisticOcean_fields_$(rank).jld2")["timeseries/u/" * string($iter)][8:end-7, 8:end-7, :]
        vd[irange, :, :] .= jldopen("RealisticOcean_fields_$(rank).jld2")["timeseries/v/" * string($iter)][8:end-7, 8:end-7, :]
    end

    set!(u, ud)
    set!(v, vd)

    compute!(ζ)

    out = interior(ζ, :, :, 1)
    out[out .== 0] .= NaN

    out
end

Tl = @lift begin
    Td = zeros(size(T))

    for rank in 0:7
        irange = iranges[rank+1]
        Td[irange, :, :] .= jldopen("RealisticOcean_fields_$(rank).jld2")["timeseries/T/" * string($iter)][8:end-7, 8:end-7, :]
    end

    set!(T, Td)

    out = interior(T, :, :, 1)
    out[out .== 0] .= NaN

    out
end

fig = Figure(resolution = (5000, 1200))
ax  = Axis(fig[1, 1])
heatmap!(ax, ζl, colorrange = (-5e-5, 5e-5)) 
ax  = Axis(fig[1, 2])
heatmap!(ax, Tl, colorrange = (0, 30), colormap = :magma)

CairoMakie.record(fig, "sim-2years.mp4", iterations, framerate = 10) do i
    iter[] = i; @info "step $i"
end