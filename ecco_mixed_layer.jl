using DataDeps
using SeawaterPolynomials
using Oceananigans.BuoyancyModels: ρ′, buoyancy_perturbationᶜᶜᶜ

include("utils.jl")
include("ecco_z_faces.jl")

arch = CPU()

@info "building grid"
z_faces = ECCO_z_faces()
ecco_grid = LatitudeLongitudeGrid(arch; size = (1440, 720, 50),
                                 longitude = (-180, 180),
                                 latitude = (-90, 90),
                                 z = z_faces)

T = CenterField(ecco_grid)
S = CenterField(ecco_grid)

eos = SeawaterPolynomials.TEOS10EquationOfState()

ρ = Field(KernelFunctionOperation{Center, Center, Center}(ρ′, ecco_grid, eos, T, S) + eos.reference_density)

monthly_days() = [1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30, 1:31, 1:30, 1:31]

all_iters  = []
all_Tfiles = []
all_Sfiles = []

@info "what about files?"
for month in 1:12
    for day in monthly_days()[month]
        push!(all_iters, string(1995) * string(month, pad=2) * string(day, pad=2))
        push!(all_Tfiles, "eccodata/THETA.1440x720x50." * all_iters[end] * ".nc")
        push!(all_Sfiles, "eccodata/SALT.1440x720x50." * all_iters[end] * ".nc")
    end
end

iters  = [all_iters[i]  for i in 1:3:length(all_iters)]
Tfiles = [all_Tfiles[i] for i in 1:3:length(all_Tfiles)]
Sfiles = [all_Sfiles[i] for i in 1:3:length(all_Sfiles)]

@info "my iterations" Tfiles Sfiles iters

mixed_layer = []
days        = [] 

for (Tfile, Sfile, day) in zip(Tfiles, Sfiles, iters)
    @info "doing files $Tfile $Sfile"
    Td = load_and_setup_ecco_data(Tfile, ecco_grid, "THETA")
    Sd = load_and_setup_ecco_data(Sfile, ecco_grid, "SALT")
    
    set!(T, Td)
    set!(S, Sd)

    compute!(ρ)
    @info "extrema of b" extrema(ρ)

    h = compute_buoyancy_mixed_layer(ρ, ecco_grid)
    @info "extrema of h" extrema(h)
    push!(mixed_layer, h)
    push!(days, day)
end

jldsave("ecco_mixed_layer.jld2", days = days, mixed_layer = mixed_layer)
