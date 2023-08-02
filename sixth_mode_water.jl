using DataDeps

include("utils.jl")
include("oceanscalingtest_faces.jl")

const λW, λE = (-81.5, -40)
const φS, φN = ( 26.0, 43.0)
const year   = 1995

arch = CPU()

@info "building grid"
z_faces = exponential_z_faces(100)
sixth_grid = LatitudeLongitudeGrid(arch; size = (2160, 900, 100),
                                 longitude = (-180, 180),
                                 latitude = (-75, 75),
                                 z = z_faces)

volumes = arch_array(arch, interior(VolumeField(sixth_grid)))

monthly_days() = [1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]

@info "what about files?"
const regex = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$";
files = readdir("/orcd/nese/raffaele/001/ssilvest/ocean_sixth/")
files = filter(x -> x[1:10] == "compressed", files)
iters = []
for file in files
    file   = file[1:end-5]
    string = ""
    i = length(file)
    while occursin(regex, "$(file[i])")
        string = file[i] * string
        i -= 1
    end
    push!(iters, string)
end
iters = sort(iters)

@info "my iterations" files iters

mode_water = []
days       = [] 

for (idx, iter) in enumerate(iters)
    day  = mod((idx - 1) * 10, 365)
    file = "compressed_iteration$(iter).jld2"
    @info "doing file $file and day $day"
    T = load_and_setup_data(file, sixth_grid)
    V = calc_mode_water(T, sixth_grid, volumes)
    push!(mode_water, V)
    push!(days, day)
end

jldsave("sixth_mode_water.jl", days = days, mode_water = mode_water)