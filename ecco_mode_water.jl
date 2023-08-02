using JLD2
using Oceananigans
using Oceananigans.Utils
using Oceananigans.Grids: architecture, node
using KernelAbstractions: @kernel, @index
using Oceananigans.AbstractOperations: GridMetricOperation
using NetCDF
using CUDA
using CairoMakie

MetricField(loc, grid, metric) = compute!(Field(GridMetricOperation(loc, metric, grid)))
VolumeField(grid, loc = (Center, Center, Center)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume)

@kernel function _flag_volumes(flag, T, grid)
    i, j, k = @index(Global, NTuple)
    λ, φ, _ = node(i, j, k, grid, Center(), Center(), Center())

    if (λW < λ < λE) && (φS < φ < φN) && (17 < T[i, j, k] < 19)
        flag[i, j, k] = 1
    end
end

function calc_mode_water_volume(T, grid, volumes)
    flag = zeros(size(T)...)
    launch!(architecture(grid), grid, :xyz, _flag_volumes, flag, T, grid)
    return sum(flag .* volumes)
end

function load_and_setup_data(file)
    data = ncread(file, "THETA")[:, :, :, 1]
    data[data .< -1e2] .= 0
    data = reverse(data, dims = 3)
    nx, ny, nz = size(data)
    tmp  = zeros(nx ÷2, ny, nz)
    transpose_flux!(data, tmp)

    return data
end

function transpose_flux!(var, tmp)
    nx = size(tmp, 1)
    tmp .= var[nx+1:end, :, :]
    var[nx+1:end, :, :] .= var[1:nx, :, :]
    var[1:nx, :, :]     .= tmp

    return nothing
end

const λW, λE = (-81.5, -40)
const φS, φN = ( 26.0, 43.0)
const year   = 1995

z_faces = ECCO_z_faces()
ecco_grid = LatitudeLongitudeGrid(size = (1440, 720, 50),
                                 longitude = (-180, 180),
                                 latitude = (-90, 90),
                                 z = z_faces)

volumes = Array(interior(VolumeField(ecco_grid)))

monthly_days() = [1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]

iters = []
files = []

for month in 1:12
    for day in monthly_days()[month]
        push!(iters, string(year) * string(month, pad=2) * string(day, pad=2))
        push!(files, "THETA.1440x720x50." * iters[end] * ".nc")
    end
end

mode_water = []
days       = [] 

for (file, day) in zip(files, days)
    T = load_and_setup_data(file)
    V = calc_mode_water(T, ecco_grid, volumes)
    push!(mode_water, V)
    push!(days, day)
end

jldsave("ecco_mode_water.jl", days = days, mode_water = mode_water)