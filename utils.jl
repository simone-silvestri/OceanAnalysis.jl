using JLD2
using Oceananigans
using Oceananigans.Utils
using Oceananigans.Architectures: arch_array
using Oceananigans.Grids: architecture, node
using KernelAbstractions: @kernel, @index
using Oceananigans.AbstractOperations: GridMetricOperation
using NetCDF
using CUDA

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
    flag = arch_array(architecture(grid), zeros(size(T)...))
    launch!(architecture(grid), grid, :xyz, _flag_volumes, flag, T, grid)
    return sum(flag .* volumes)
end

function load_and_setup_data(filepath, grid)
    file = jldopen(filepath)
    data = file["T"]
    return arch_array(architecture(grid), data)
end

function load_and_setup_ecco_data(filepath, grid)
    data = ncread(filepath, "THETA")[:, :, :, 1]
    data[data .< -1e2] .= 0
    data = reverse(data, dims = 3)
    nx, ny, nz = size(data)
    tmp  = zeros(nx ÷2, ny, nz)
    transpose_flux!(data, tmp)

    return arch_array(architecture(grid), data)
end

function transpose_flux!(var, tmp)
    nx = size(tmp, 1)
    tmp .= var[nx+1:end, :, :]
    var[nx+1:end, :, :] .= var[1:nx, :, :]
    var[1:nx, :, :]     .= tmp

    return nothing
end
