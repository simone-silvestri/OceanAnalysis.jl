using DataDeps


function ECCO_z_faces()

    path = "https://github.com/CliMA/OceananigansArtifacts.jl/raw/ss/new_hydrostatic_data_after_cleared_bugs/quarter_degree_near_global_input_data/"

    dh = DataDep("quarter_degree_near_global_lat_lon",
      "Forcing data for global latitude longitude simulation",
       path * "z_faces-50-levels.jld2"
    )

    DataDeps.register(dh)

    datadep"quarter_degree_near_global_lat_lon"
    
    datadep_path = @datadep_str "quarter_degree_near_global_lat_lon/z_faces-50-levels.jld2"
    file_z_faces = jldopen(datadep_path)
    
    # Stretched faces taken from ECCO Version 4 (50 levels in the vertical)
    return file_z_faces["z_faces"];
end