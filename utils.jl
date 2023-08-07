using JLD2
using Oceananigans
using Oceananigans.Utils
using Oceananigans.Architectures: arch_array
using Oceananigans.Grids: architecture, node
using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.Loopinfo: @unroll
using Oceananigans.AbstractOperations: GridMetricOperation
using NetCDF
using CUDA

MetricField(loc, grid, metric) = compute!(Field(GridMetricOperation(loc, metric, grid)))
VolumeField(grid, loc = (Center, Center, Center)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume)

@kernel function _flag_volumes!(flag, T, grid)
    i, j, k = @index(Global, NTuple)
    λ, φ, _ = node(i, j, k, grid, Center(), Center(), Center())

    if (λW < λ < λE) && (φS < φ < φN) && (17 < T[i, j, k] < 19)
        flag[i, j, k] = 1
    end
end

function calc_mode_water_volume(T, grid, volumes)
    flag = arch_array(architecture(grid), zeros(size(T)...))
    launch!(architecture(grid), grid, :xyz, _flag_volumes!, flag, T, grid)
    return sum(flag .* volumes)
end

@kernel function _compute_mixed_layer!(h, b, grid, t)
    i, j = @index(Global, NTuple)
    k    = grid.Nz
    @inbounds begin    
        b★ = 0.5 * (b[i, j, k] + b[i, j, k-1]) 
        k  = k-2
        h[i, j] = zero(grid)
        while abs(b[i, j, k] - b★) < t && !(b[i, j, k] == NaN || b★ == NaN) && k > 1
            h[i, j] = znode(i, j, k, grid, Center(), Center(), Center())
            k -= 1
        end
    end
end

@kernel function _fill_NaNs!(b)
    i, j, k = @index(Global, NTuple)
    @inbounds begin
        if b[i, j, k] == 0
            b[i, j, k] = NaN
        end
    end
end

function compute_buoyancy_mixed_layer(b, grid; threshold = 0.0003)
    h = zeros(size(grid, 1), size(grid, 2)) 
    launch!(architecture(grid), grid, :xyz, _fill_NaNs!, b)
    launch!(architecture(grid), grid, :xy, _compute_mixed_layer!, h, grid, b, threshold)
    return h
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
