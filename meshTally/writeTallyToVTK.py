import openmc
statepoint = openmc.StatePoint("./statepoint.100.h5")

mesh = statepoint.meshes[1]

print(statepoint.meshes)

flux_tally = statepoint.get_tally(id=30)

# Shape will depend on filters present, order of filters, etc.
# In this case, the mesh filter was the only one applied, and the mesh
# dimensions make up the first three axes of the tally (r, z, phi I believe)
flux_mean = flux_tally.get_reshaped_data("mean", expand_dims=True)[
    :, :, :, 0, 0
]
flux_std_dev = flux_tally.get_reshaped_data("std_dev", expand_dims=True)[
    :, :, :, 0, 0
]

fission_mean = flux_tally.get_reshaped_data("mean", expand_dims=True)[
    :, :, :, 0, 2
]

mesh_data = {
    "Flux (Mean) [cm]": flux_mean,
    "Flux (Std. Dev.) [cm]": flux_std_dev,
    "Fission (Mean) [cm]": fission_mean,
}

flux_mean = flux_tally.get_reshaped_data("mean", expand_dims=True)

print(flux_mean.shape)

mesh.write_data_to_vtk(
    filename=f"flux.vtk",
    datasets=mesh_data,
    volume_normalization=False,
)