using JLD2
using Oceananigans
using Oceananigans.Utils
using KernelAbstractions: @kernel, @index
using NetCDF
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
    flag = CenterField(grid)
    launch!(architecture(grid), grid, _flag_volumes, flag, T, grid)

    return sum(flag * volumes)
end

function load_and_setup_data(file, grid)
    data = ncread(file, "THETA")[:, :, :, 1]
    data[data .< -1e2] .= 0
    data = reverse(data, dims = 3)
    nx, ny, nz = size(data)
    tmp  = zeros(nx ÷2, ny, nz)
    transpose_flux!(data, tmp)

    return set!(CenterField(grid), data)
end

function transpose_flux!(var, tmp)
    nx = size(tmp, 1)
    tmp .= var[nx+1:end, :, :]
    var[nx+1:end, :, :] .= var[1:nx, :, :]
    var[1:nx, :, :]     .= tmp
end

const λW, λE = ( 26,   43)
const φS, φN = (-81.5, -5.6)

z_faces = ECCO_z_faces()
ecco_grid = LatitudeLongitudeGrid(size = (1440, 720, 50),
                                 longitude = (-180, 180),
                                 latitude = (-90, 90),
                                 z = z_faces)

volumes = VolumeField(ecco_grid)

V_all = []

for file in all_files
    T = load_and_setup_data(file, ecco_grid)
    V = calc_mode_water(T, ecco_grid, volumes)
    push!(V_all, V)
end
