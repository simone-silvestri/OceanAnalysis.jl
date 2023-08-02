using DataDeps

include("ecco_z_faces.jl")

const λW, λE = (-81.5, -40)
const φS, φN = ( 26.0, 43.0)
const year   = 1995

arch = GPU()

z_faces = ECCO_z_faces()
ecco_grid = LatitudeLongitudeGrid(arch; size = (1440, 720, 50),
                                 longitude = (-180, 180),
                                 latitude = (-90, 90),
                                 z = z_faces)

volumes = arch_array(arch, interior(VolumeField(ecco_grid)))

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
    T = load_and_setup_data(file, ecco_grid)
    V = calc_mode_water(T, ecco_grid, volumes)
    push!(mode_water, V)
    push!(days, day)
end

jldsave("ecco_mode_water.jl", days = days, mode_water = mode_water)