# Matthew Nyberg
# Subcritical Systems Inc. Model

import openmc
import sys
import numpy as np

### Materials ### 

U_mat = openmc.Material(20, 'Uranium')
U_mat.add_nuclide('U235', 0.1, 'wo')
U_mat.add_nuclide('U238', 0.9, 'wo')

O_mat = openmc.Material(21, 'Oxygen')
O_mat.add_element('O', 1, 'ao')

# UO2_mat = openmc.Material(1, "UO2 10w% U-235", temperature=900)
UO2_mat = openmc.Material.mix_materials([U_mat, O_mat], [1/3,2/3], 'ao')
UO2_mat.set_density('g/cm3', 10.0)

Pb_mat = openmc.Material(2, "Lead-208", temperature=900)
Pb_mat.add_nuclide('Pb208', 1)
Pb_mat.set_density('g/cm3', 10.0)

Fe_mat = openmc.Material(3, 'Fe-56', temperature=900)
Fe_mat.add_nuclide('Fe56', 1)
Fe_mat.set_density('g/cm3', 8.0)

center_hole = openmc.Material(4, 'Central Void', temperature=900)
center_hole.add_element('O', 1)
center_hole.set_density('g/cm3', 0.1)


# Do math to determine weight ratios and mix materials
volume_fractions = {'UO2':0.4, 'Fe-56':0.1, 'Pb-208':0.5}
weights = [volume_fractions['UO2']*UO2_mat.density, volume_fractions['Pb-208']*Pb_mat.density, volume_fractions['Fe-56']*Fe_mat.density]
weight_fractions = np.divide(weights,sum(weights))

# Mix materials into homogenized material
homogenized_mat = openmc.Material.mix_materials([UO2_mat, Pb_mat, Fe_mat], weight_fractions, 'wo')
# TODO check density math
hom_density = sum(weights)
homogenized_mat.set_density('g/cm3', hom_density)
# TODO check temperature assumption
homogenized_mat.temperature = 900

# Add materials
all_mats = openmc.Materials()
all_mats += [UO2_mat, Pb_mat, Fe_mat, homogenized_mat, center_hole]
all_mats.export_to_xml()

###### GEOMETRIES ##########

# Given properties/dimensions
central_hole_diam = 30

### Outputs from search for keff 0.95 (0.94982 +/- 0.00027)
# Annular Height = 241.0842316614915
# Outer Core Radius = 180.8131737461186
ann_height = float(sys.argv[1])
fuel_outer_diam = ann_height*1.5

# Create planes for inner central hole
hole_outer = openmc.ZCylinder(r=central_hole_diam/2)
top = openmc.ZPlane(z0=ann_height/2, boundary_type='vacuum')
bottom = openmc.ZPlane(z0=-ann_height/2, boundary_type='vacuum')

hole_cell = openmc.Cell(name="Central Hole (currently Void)")
hole_region = -hole_outer & -top & +bottom
hole_cell.region = hole_region
# hole_cell.fill = center_hole

# Create homogenous fuel cylinder
fuel_outer = openmc.ZCylinder(r=fuel_outer_diam/2, boundary_type='vacuum')
bulkFuel = openmc.Cell(name='bulkFuel')
bulkFuel_region = +hole_outer & -fuel_outer & -top & +bottom
bulkFuel.fill = homogenized_mat
bulkFuel.region = bulkFuel_region

# Create void cube
# box = openmc.model.RectangularPrism(width=1000, height=1000, boundary_type='reflective')
box = openmc.model.RectangularPrism(width=2000, height=2000, boundary_type='vacuum')
void_region = -box & ~bulkFuel_region & ~hole_region
void = openmc.Cell(name='void')
void.region = void_region

# Create a universe using the fuel cyl+void cube
root_universe = openmc.Universe(cells=(hole_cell, bulkFuel, void))
geometry = openmc.Geometry()
geometry.root_universe = root_universe
geometry.export_to_xml()

###### TALLIES ######

### Mesh Used in Tallies ### 
cyl_mesh = openmc.CylindricalMesh(mesh_id=1, name="Full fuel mesh",z_grid=[], r_grid=[])
n_r_bins = 100
n_z_bins = 100
height = ann_height

cyl_mesh.r_grid = np.zeros(n_r_bins+1)
cyl_mesh.z_grid = np.zeros(n_z_bins+1)
for i in range(0,n_r_bins+1):
    cyl_mesh.r_grid[i] = 0+i*(fuel_outer_diam/2)/n_r_bins
    # cyl_mesh.r_grid[i] = radii[1]+i*width[2]/n_r_bins
for i in range(0,n_z_bins+1):
    cyl_mesh.z_grid[i] = -0.5*height+i*height/n_z_bins

meshVols = cyl_mesh.volumes
radMeshVols = np.zeros((100,1))
for i in range(0,n_r_bins):
    radMeshVols[i,0] = sum(sum(meshVols[i]))

cell_filter = openmc.CellFilter(bulkFuel)
tally1 = openmc.Tally(1)
tally1.filters = [cell_filter]
tally1.scores = ['nu-fission', 'scatter', 'total', 'fission', 'absorption']

# surface_filter = openmc.SurfaceFilter(mult_outer, fuel_outer, reflect_outer, top, bottom)
surface_filter = openmc.SurfaceFilter(hole_outer)
tally_current1 = openmc.Tally(3)
tally_current1.filters = [surface_filter]
tally_current1.scores = ['current']

# surface_filter = openmc.SurfaceFilter(mult_outer, fuel_outer, reflect_outer, top, bottom)
surface_filter = openmc.SurfaceFilter(bottom)
tally_current2 = openmc.Tally(4)
tally_current2.filters = [surface_filter]
tally_current2.scores = ['current']

# surface_filter = openmc.SurfaceFilter(mult_outer, fuel_outer, reflect_outer, top, bottom)
surface_filter = openmc.SurfaceFilter(top)
tally_current3 = openmc.Tally(5)
tally_current3.filters = [surface_filter]
tally_current3.scores = ['current']

# surface_filter = openmc.SurfaceFilter(mult_outer, fuel_outer, reflect_outer, top, bottom)
surface_filter = openmc.SurfaceFilter(fuel_outer)
tally_current4 = openmc.Tally(6)
tally_current4.filters = [surface_filter]
tally_current4.scores = ['current']

tally10 = openmc.Tally(10)
tally10.scores = ['heating', 'fission', 'absorption', 'scatter', '(n,2n)','(n,3n)', 'total']

mesh_tally = openmc.Tally(30)
mesh_tally.filters = [openmc.MeshFilter(cyl_mesh)]
mesh_tally.scores = ['flux', 'heating', 'fission', '(n,2n)']

# Export tallies
tallies = openmc.Tallies([tally1, tally_current1, tally_current2, tally_current3, tally_current4, tally10, mesh_tally])
tallies.export_to_xml()

# Set the source/k-code settings
settings=openmc.Settings()
settings.run_mode = 'fixed source'
# settings.run_mode = 'eigenvalue'
settings.temperature = {'method' : 'interpolation'}
# Set the source type and location

# Fixed source settings
if settings.run_mode == 'fixed source':
    # Assuming point spatial distribution
    point = openmc.stats.Point((0,0,0))
    # TODO more complex spatial distribution (based on https://pdf.sciencedirectassets.com/277348/1-s2.0-S1875389217X00059/1-s2.0-S1875389217301943/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEE0aCXVzLWVhc3QtMSJHMEUCIGXl3d9d%2BxHgMceS%2BkRXBeU0HNXcCbUdKmYqBfcbbeYPAiEAqjJXpU1p%2BPZ5Ixo82itgh06BINpZcSPXJjMPFVKmRqoqsgUIFhAFGgwwNTkwMDM1NDY4NjUiDPqlUbjiGlJ8PrX0HSqPBe8Hr0o6Y8UuuLX4FPCwAgQw0j6faHLIXtuAgnlTDQESbsKaZ1jp6jmDr0qwyquhef3Tmk%2FzZ0zbQmPH2yMmFR01xJ2Bst0Hw3JD%2BtpNCOPXrIMq6HFkW%2BabDt0TXk0Eu1i2Kz5v31ko7JQEDnX6jmIVR1WDuKHpgOh%2Bnxh1oRmEh3lqAj%2FVGSctYJw4dTtCllLIZrk8f3KYHeX5PnwfPT3mlGihZEZEzaMrH0ayiratbLtnRhJx3UqihXyahhBIYxqHtKfb5EEkaKoLAZNKUsVs1ld%2FzWJWmXaEXtYaYmEPx0ei1x5H0cUHwb3LSFzLHwkH9Ql0naAkLBmJ0HvHvqSIKshG7hfK1QXBciUnuXu1YGal1OQEU3pye8JnozTRyKmdxg%2BuY0GheJIk4icpt%2Fuyhwki8CH2OUcyXKzjfknd0ABKNJphKQix7ld%2B2t%2BjaXut9BOQNpf%2FYzOQpadpqLQyLJGk%2Fa3s%2Fs4g89MNWyIUr4j8sB9Fueiv7%2Bst8JQq5SbFvzHEN%2FOz25KprrY74tO4wi2w0OjBTkEVraQUQPSOWwJqWRlepwfjUyo3NcKKpgbosXReATcmRc7%2BakxY%2FFue0o%2BkRMDHV6omqXYTM4iZhwmOY2YK6AG16i8q9aXCTFc0HavECYmu8nnojF5FeaSjc37JLNTUQJnENRagdH5XOSvEpEgWSzdG2f%2FdnOyNCyh1Q39qVcPeI5wvb56aJtFdkxoNMJuOfeMQxlA%2FGlZX54arrjJNEikDP5YraMfBaZEz594GFXN6DPHmCjG9CWEY%2BPiNSKcLBWabMGExkOxAK3rRkK0Xi63LKJGCVt6Nv1wKzQV8U7ZtEf%2F0aj6l%2Bv65zrucBGnUKnNCyVZJ0vswqueSyAY6sQGWbhUsNeazMdIiprICW07Y1Pze0decigvo2OcimO%2B3nygYpAcAdqsOQhBJv2%2B07wcNoqa1YfcypmcAxQtKwoFDsMMHxKQbEjcq4GifEckFAVuBtJB0VkW8r67iEjF8XD%2Bwz5i1GR%2FJHQMd6VnMTvMZk%2BZj0A7y69%2B7IiyjupGCyO0Y3DVw5I9C4TRTPWFCRDANqnTAHPm87rZt3qhi6UvSHAo5%2FIRCWNUdOVlOvjGAeP4%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20251031T133121Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYWS7BLQON%2F20251031%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=11c54795491bf83177fea44769ade2e624123bfb71d6fd43b47a13f74ea50b41&hash=a2cff02c1b171afb38040c49c4b0b040e12600416ebc0fbf5641bd7d47959da9&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S1875389217301943&tid=spdf-0bbd295d-d7e8-4cf7-b7c0-cd0734f64b2f&sid=3650e0e45ea7044a272908859dd4a1d37c8agxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f165e5f50505e56565550&rr=9973832469abf861&cc=us):
    r=openmc.stats.Uniform(-5,5)
    phi=openmc.stats.Uniform(0, 6.283185307179586)
    z=openmc.stats.Tabular(x = [5, 2.5, 0, -2.5, -5], p = [100, 10, 1, 0.1, 0.01]) # , interpolation='log-log')
    cylASource = openmc.stats.CylindricalIndependent(r, phi, z, origin=(0.0, 0.0, 0.0))
    z=openmc.stats.PowerLaw(a=1, b=11, n=-1)
    # NOTE: This is flipped (source coming in from negative y). Doesn't matter now since the geometry is symmetric about z=0
    cylASource_powerLaw = openmc.stats.CylindricalIndependent(r, phi, z, origin=(0.0, 0.0, -6))
    # TODO non-isotropic angular dist:
    # IAEA resource cites 3-4 MeV: https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://inis.iaea.org/collection/NCLCollectionStore/_Public/43/099/43099436.pdf&ved=2ahUKEwip_tW4xM6QAxVFOTQIHQFMCKQQFnoECBsQAQ&usg=AOvVaw1lq7wNWCDT6EYFJmX8UFqy
    energyDist_simple = openmc.stats.Normal(mean_value=3.5e6, std_dev=3e5)
    # Further spectra (energy and energy dependent on angle) here: https://pdf.sciencedirectassets.com/277348/1-s2.0-S1875389217X00059/1-s2.0-S1875389217301943/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEE0aCXVzLWVhc3QtMSJHMEUCIGXl3d9d%2BxHgMceS%2BkRXBeU0HNXcCbUdKmYqBfcbbeYPAiEAqjJXpU1p%2BPZ5Ixo82itgh06BINpZcSPXJjMPFVKmRqoqsgUIFhAFGgwwNTkwMDM1NDY4NjUiDPqlUbjiGlJ8PrX0HSqPBe8Hr0o6Y8UuuLX4FPCwAgQw0j6faHLIXtuAgnlTDQESbsKaZ1jp6jmDr0qwyquhef3Tmk%2FzZ0zbQmPH2yMmFR01xJ2Bst0Hw3JD%2BtpNCOPXrIMq6HFkW%2BabDt0TXk0Eu1i2Kz5v31ko7JQEDnX6jmIVR1WDuKHpgOh%2Bnxh1oRmEh3lqAj%2FVGSctYJw4dTtCllLIZrk8f3KYHeX5PnwfPT3mlGihZEZEzaMrH0ayiratbLtnRhJx3UqihXyahhBIYxqHtKfb5EEkaKoLAZNKUsVs1ld%2FzWJWmXaEXtYaYmEPx0ei1x5H0cUHwb3LSFzLHwkH9Ql0naAkLBmJ0HvHvqSIKshG7hfK1QXBciUnuXu1YGal1OQEU3pye8JnozTRyKmdxg%2BuY0GheJIk4icpt%2Fuyhwki8CH2OUcyXKzjfknd0ABKNJphKQix7ld%2B2t%2BjaXut9BOQNpf%2FYzOQpadpqLQyLJGk%2Fa3s%2Fs4g89MNWyIUr4j8sB9Fueiv7%2Bst8JQq5SbFvzHEN%2FOz25KprrY74tO4wi2w0OjBTkEVraQUQPSOWwJqWRlepwfjUyo3NcKKpgbosXReATcmRc7%2BakxY%2FFue0o%2BkRMDHV6omqXYTM4iZhwmOY2YK6AG16i8q9aXCTFc0HavECYmu8nnojF5FeaSjc37JLNTUQJnENRagdH5XOSvEpEgWSzdG2f%2FdnOyNCyh1Q39qVcPeI5wvb56aJtFdkxoNMJuOfeMQxlA%2FGlZX54arrjJNEikDP5YraMfBaZEz594GFXN6DPHmCjG9CWEY%2BPiNSKcLBWabMGExkOxAK3rRkK0Xi63LKJGCVt6Nv1wKzQV8U7ZtEf%2F0aj6l%2Bv65zrucBGnUKnNCyVZJ0vswqueSyAY6sQGWbhUsNeazMdIiprICW07Y1Pze0decigvo2OcimO%2B3nygYpAcAdqsOQhBJv2%2B07wcNoqa1YfcypmcAxQtKwoFDsMMHxKQbEjcq4GifEckFAVuBtJB0VkW8r67iEjF8XD%2Bwz5i1GR%2FJHQMd6VnMTvMZk%2BZj0A7y69%2B7IiyjupGCyO0Y3DVw5I9C4TRTPWFCRDANqnTAHPm87rZt3qhi6UvSHAo5%2FIRCWNUdOVlOvjGAeP4%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20251031T133121Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYWS7BLQON%2F20251031%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=11c54795491bf83177fea44769ade2e624123bfb71d6fd43b47a13f74ea50b41&hash=a2cff02c1b171afb38040c49c4b0b040e12600416ebc0fbf5641bd7d47959da9&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S1875389217301943&tid=spdf-0bbd295d-d7e8-4cf7-b7c0-cd0734f64b2f&sid=3650e0e45ea7044a272908859dd4a1d37c8agxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f165e5f50505e56565550&rr=9973832469abf861&cc=us 
    energyList = np.multiply(1e6,[1, 2, 3, 4, 5, 10, 20, 100, 500])
    energyDist_tabular = openmc.stats.Tabular(x = energyList, p = [2e5, 1.5e5, 1e5, 5e4, 5e4, 1e4, 0.9e4, 0.3e4, 2e2]) # , interpolation='log-log')
    source = openmc.Source(space=cylASource, energy=energyDist_tabular)
    settings.source = source
    settings.batches = 100
    settings.particles = 10000
    # settings.create_fission_neutrons = False
    settings.create_fission_neutrons = True

# k-code settings
if settings.run_mode == 'eigenvalue':
    r=openmc.stats.Uniform(central_hole_diam/2,fuel_outer_diam/2)
    phi=openmc.stats.Uniform(0, 6.283185307179586)
    z=openmc.stats.Uniform(-ann_height/4,ann_height/4)
    cylFirstSource = openmc.stats.CylindricalIndependent(r, phi, z, origin=(0.0, 0.0, 0.0))
    source = openmc.Source(space=cylFirstSource)
    settings.batches = 250
    settings.inactive = 100
    settings.particles = 10000

settings.export_to_xml()

#### General Model Holder #### 
myModel = openmc.Model()
myModel.materials=all_mats
myModel.geometry=geometry
myModel.tallies=tallies
myModel.settings=settings
myModel.export_to_model_xml()

openmc.run()

