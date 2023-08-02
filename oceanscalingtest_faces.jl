
@inline exponential_profile(z; Lz, h) = (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) 

function exponential_z_faces(Nz, Depth=5244.5; h = Nz / 4.5)

    z_faces = exponential_profile.((1:Nz+1); Lz = Nz, h)

    # Normalize
    z_faces .-= z_faces[1]
    z_faces .*= - Depth / z_faces[end]
    
    z_faces[1] = 0.0

    return reverse(z_faces)
end