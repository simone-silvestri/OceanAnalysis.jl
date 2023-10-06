using DataDeps

include("utils.jl")
include("oceanscalingtest_faces.jl")

arch = CPU()

@info "building grid"
z_faces = exponential_z_faces(100)
twelth_grid = LatitudeLongitudeGrid(arch; size = (4320, 1800, 100),
                                 longitude = (-180, 180),
                                 latitude = (-75, 75),
                                 z = z_faces)

bathymetry = jldopen("bathymetry12.jld2")["bathymetry"]
grid = ImmersedBoundaryGrid(twelth_grid, GridFittedBottom(bathymetry))

T = CenterField(grid)
η = Field((Center, Center, Nothing), grid)
vol = VolumeField(grid)

tm = mean(T, dims = 1)
fill!(tm, 0.0)

folder = "/nobackup/users/ssilvest/perlmutter-test/OceanScalingTests.jl/solution-to-transfer/"

files, iters = get_iters(folder)

for (idx, iter) in enumerate(iters)
    file = folder * "compressed_iteration_$(iter).jld2"
    @info "doing file $file and iter $iter"
    set!(T, jldopen(file)["T"])
    tm .+= sum(T * vol, dims = 1) ./ length(iters)
     η .+= file["η"] ./ length(iters)
end

jldsave("mean_temp.jld2", Tm = tm, ηm = η)