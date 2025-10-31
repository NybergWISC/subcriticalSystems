# Matthew Nyberg
# Subcritical Systems Inc. Model

import openmc
import numpy as np
import math

def build_model(ann_height, **model_args):

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
    # ann_height = 300
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
    # root_universe = openmc.Universe(cells=(hole_cell, bulkFuel, void))
    root_universe = openmc.Universe(cells=(hole_cell, bulkFuel))
    geometry = openmc.Geometry()
    geometry.root_universe = root_universe
    geometry.export_to_xml()

    ###### TALLIES ######

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

    # Export tallies
    tallies = openmc.Tallies([tally1, tally_current1, tally_current2, tally_current3, tally_current4, tally10])
    tallies.export_to_xml()

    # Set the source/k-code settings
    settings=openmc.Settings()
    # settings.run_mode = 'fixed source'
    settings.run_mode = 'eigenvalue'
    settings.temperature = {'method' : 'interpolation'}
    # Set the source type and location

    # Fixed source settings
    if settings.run_mode == 'fixed source':
        point = openmc.stats.Point((0,0,0))
        # TODO check this mean value
        energyDist = openmc.stats.Normal(mean_value=1e6, std_dev=1e4)
        source = openmc.Source(space=point, energy=energyDist)
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
        settings.batches = 300
        settings.inactive = 100
        settings.particles = 25000

    settings.export_to_xml()

    #### General Model Holder #### 
    myModel = openmc.Model()
    myModel.materials=all_mats
    myModel.geometry=geometry
    myModel.settings=settings

    # Calculate total uranium mass
    totalUO2Mass = weights[0]*math.pi*ann_height*((fuel_outer_diam/2)**2-(central_hole_diam/2)**2)
    U238atPerc = (0.9/238)/((0.1/235)+(0.9/238))
    U235atPerc = (0.1/235)/((0.1/235)+(0.9/238))
    avgUmm = U235atPerc*0.1+U238atPerc*0.9
    totalUMass = totalUO2Mass*avgUmm/16.0
    print(totalUMass)

    return myModel


startingGuess = 241.62210539432047
min_height = startingGuess - 20
max_height = startingGuess + 20

try:
    ann_height, guesses, keffs = openmc.search_for_keff(build_model, target=0.95, bracket=[min_height, max_height], bracketed_method = "brentq", tol=1e-4, print_iterations=True, model_args={"test":"TESTING PRINT"})
    print(f'Annular Height = {ann_height}')
    outerCoreRadius = (ann_height*1.5)/2
    print(f'Outer Core Radius = {outerCoreRadius}')
except ValueError: 
    print("Some Error Occured, defaulting to 300cm height")
    ann_height=300

# openmc.run()
