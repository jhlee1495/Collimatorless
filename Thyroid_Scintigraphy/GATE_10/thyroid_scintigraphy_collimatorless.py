#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import opengate as gate
# import opengate.contrib.spect.genm670 as spect_ge_nm670
import pathlib
from pathlib import Path
from scipy.spatial.transform import Rotation
from opengate.geometry.utility import get_grid_repetition, get_circular_repetition
import opengate.geometry.materials as gmat
import uproot
import pandas as pd
from scipy import io
import os

# help(gate.geometry.utility.get_grid_repetition)
# help(gate.geometry.utility.get_circular_repetition)
# help(gate.image.get_translation_between_images_center)
current_path = pathlib.Path(__file__).parent.resolve()
output_path = current_path / "output"
output_file = output_path / "output.mhd"

# colors
invisible = [0, 0, 0, 0]
red = [1, 0, 0, 1]
blue = [0, 0, 1, 1]
green = [0, 1, 0, 1]
yellow = [0.9, 0.9, 0.3, 1]
gray = [0.5, 0.5, 0.5, 1]
white = [1, 1, 1, 0.8]
purple = [0.502, 0, 0.502, 1]

file_number = 3

if __name__ == "__main__":
    # create the simulation
    sim = gate.Simulation()

    # main options
    sim.g4_verbose = False
    sim.visu = True
    # sim.visu_type = "gdml"
    sim.visu_type = "vrml"
    # sim.visu_type = "qt"
    sim.number_of_threads = 96
    sim.random_seed = "auto"

    # units
    m = gate.g4_units.m
    sec = gate.g4_units.second
    hour = gate.g4_units.hour
    days = 3600 * 24 * sec
    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    nm = gate.g4_units.nm
    MeV = gate.g4_units.MeV
    keV = gate.g4_units.keV
    Bq = gate.g4_units.Bq
    kBq = 1000 * Bq
    MBq = 1000 * kBq
    deg = gate.g4_units.deg
    g_cm3 = gate.g4_units.g_cm3

    # materials
    sim.volume_manager.add_material_database(
        current_path / "data/GateMaterials_XCAT.db"
    )

    # set the world size
    world = sim.world
    world.size = [3 * m, 3 * m, 3 * m]
    world.material = "Air"

        # CT image
    XCAT_phantom = sim.add_volume("Image", "XCAT_phantom")
    if sim.visu == True:
        # XCAT_phantom.image = current_path / "data/XCAT_phantom_visu.mhd"
        XCAT_phantom.image = current_path / "data/XCAT_phantom_visu.mhd"
    else:
        XCAT_phantom.image = current_path / "data/XCAT_phantom.mhd"
    XCAT_phantom.material = "Air"  # material used by default
    tol = 0.1 * g_cm3
    mat_table = current_path / "data/Schneider2000MaterialsTable.txt"
    density_table = current_path / "data/Schneider2000DensitiesTable.txt"
    (
        XCAT_phantom.voxel_materials,
        materials,
    ) = gmat.HounsfieldUnit_to_material(sim, tol, mat_table, density_table)
        # XCAT_phantom.dump_label_image = current_path / "output/labels.mhd"
    XCAT_phantom.translation = [0 * mm, 0 * mm, 0 * mm]

    # XCAT_phantom.dump_label_image = (
    #     current_path / "output/XCAT_phantom_HU_label.mhd"
    # )

    # add IEC phantom
    # iec_phantom = gate_iec.add_iec_phantom(sim, "iec")
    # iec_phantom.translation = [0 * cm, 0 * cm, 0 * cm]

    # Detector
    SPECT_outer = sim.add_volume("Box", "SPECT_outer")
    SPECT_outer.mother = world.name
    SPECT_outer.size = [30 * cm, 5 * cm, 21 * cm]
    SPECT_outer.rotation = None
    SPECT_outer.material = "Lead"
    SPECT_outer.color = green
    SPECT_outer.translation = [0, -44.925 * cm, -5 * cm]

    SPECT_inner = sim.add_volume("Box", "SPECT_inner")
    SPECT_inner.mother = SPECT_outer.name
    SPECT_inner.size = [29.7 * cm, 4.7 * cm, 20.7 * cm]
    SPECT_inner.rotation = None
    SPECT_inner.material = "Air"
    SPECT_inner.color = blue
    SPECT_inner.translation = [0, 0.2 * cm, 0 * cm]

    Detector = sim.add_volume("Box", "Detector")
    Detector.mother = SPECT_inner.name
    Detector.size = [28 * cm, 1 * cm, 19 * cm]
    Detector.rotation = None
    Detector.material = "NaITl"
    Detector.color = red
    Detector.translation = [0, 1.8 * cm, 0 * cm]

    # pixel = sim.add_volume("Box", "pixel")
    # pixel.mother = SPECT_inner.name
    # pixel.size = [4.75 * mm, 1 * cm, 4.667 * mm]
    # translations_grid = get_grid_repetition(size=[60, 1, 40], spacing=[4.75 * mm,  0, 4.667 * mm])
    # pixel.translation = translations_grid

    # Back-Compartment
    compartment = sim.add_volume("Box", "compartment")
    compartment.mother = SPECT_inner.name
    compartment.size = [28 * cm, 2.5 * cm, 19 * cm]
    compartment.rotation = None
    compartment.material = "Glass"
    compartment.color = gray
    compartment.translation = [0, 0.05 * cm, 0 * cm]

    Collimator_outer = sim.add_volume("Cons", "Collimator_outer")
    Collimator_outer.mother = world.name
    translations_Collimator_outer, rotations_Collimator_outer = (
        get_circular_repetition(
            number_of_repetitions=1,
            start_angle_deg=90,
            first_translation=[0 * cm, -5 * cm, -28.425 * cm],
            axis=[1, 0, 0],
        )
    )
    Collimator_outer.rmin1 = 0 * mm
    Collimator_outer.rmax1 = 3.45 * cm
    Collimator_outer.rmin2 = 0 * mm
    Collimator_outer.rmax2 = 20.3 * cm
    Collimator_outer.dz = 14 * cm
    Collimator_outer.sphi = 0 * deg
    Collimator_outer.dphi = 360 * deg
    Collimator_outer.translation = [0 * cm, -28.425 * cm, -5 * cm]
    Collimator_outer.rotation = Rotation.from_euler(
        "x", 90, degrees=True
    ).as_matrix()
    # Collimator_outer.rotation =  rotations_Collimator_outer
    Collimator_outer.material = "Lead"
    Collimator_outer.color = yellow

    Collimator_inner = sim.add_volume("Cons", "Collimator_inner")
    Collimator_inner.mother = Collimator_outer.name
    translations_Collimator_inner, rotations_Collimator_inner = (
        get_circular_repetition(
            number_of_repetitions=1,
            start_angle_deg=90,
            first_translation=[0 * cm, 0 * cm, 0 * cm],
            axis=[0, 0, 0],
        )
    )
    Collimator_inner.rmin1 = 0 * mm
    Collimator_inner.rmax1 = 3.35 * cm
    Collimator_inner.rmin2 = 0 * mm
    Collimator_inner.rmax2 = 20 * cm
    Collimator_inner.dz = 14 * cm
    Collimator_inner.sphi = 0 * deg
    Collimator_inner.dphi = 360 * deg
    Collimator_inner.translation = [0 * cm, 0 * cm, 0 * cm]
    Collimator_inner.rotation = None
    Collimator_inner.material = "Air"
    Collimator_inner.color = purple
    # spect, crystal = spect_ge_nm670.add_ge_nm67_spect_head(
    #     sim, "spect", collimator_type="megp", debug=sim.visu
    # )
    # spect.translation = [0, 0, -50 * cm]

    # spect digitizer channels
    channels = [
        {"name": f"spectrum", "min": 120 * keV, "max": 160 * keV},
    ]

    # spect digitizer : Hits + Adder + EneWin + Projection
    # Hits
    hc = sim.add_actor("DigitizerHitsCollectionActor", f"Hits_{Detector.name}")
    hc.mother = Detector.name
    hc.output = current_path / Path("output") / "scintigraphy_1.root"
    hc.attributes = [
        "PostPosition",
        "TotalEnergyDeposit",
        "PreStepUniqueVolumeID",
        "GlobalTime",
        "LocalTime",
        "StepLength",
        "TrackLength",
    ]
    # list of attributes :https://opengate-python.readthedocs.io/en/latest/user_guide.html#actors-and-filters

    # Singles
    sc = sim.add_actor("DigitizerAdderActor", f"Singles_{Detector.name}")
    sc.mother = hc.mother
    sc.input_digi_collection = hc.name
    sc.policy = "EnergyWinnerPosition"
    sc.output = hc.output

    # energy windows
    cc = sim.add_actor(
        "DigitizerEnergyWindowsActor", f"EnergyWindows_{Detector.name}"
    )
    cc.mother = sc.mother
    cc.input_digi_collection = sc.name
    cc.channels = channels
    cc.output = hc.output

    # projection image
    proj = sim.add_actor("DigitizerProjectionActor", f"Projection_{Detector.name}")
    proj.mother = cc.mother
    proj.input_digi_collections = [x["name"] for x in cc.channels]
    proj.spacing = [0.8906 * mm, 0.6914 * mm]
    proj.size = [314, 274]
    proj.output = current_path / Path("output") / "projection_1.mhd"
    proj.detector_orientation_matrix = Rotation.from_euler("x", 90).as_matrix()

    # Activity source from an image
    source = sim.add_source("VoxelsSource", "vox")
    source.mother = XCAT_phantom.name
    source.particle = "gamma"
    if sim.visu == True:
        source.activity = 10 * Bq / sim.number_of_threads
    else:
        source.activity = 5.18e6 * Bq / sim.number_of_threads
    # source.activity = 4013 * Bq / sim.number_of_threads
    source.image = current_path / Path("data") / "XCAT_phantom_act_1.mhd"
    source.direction.type = "iso"
    source.energy.mono = 140 * keV
    source.direction.theta = [60 * deg, 120 * deg]
    source.direction.phi = [60 * deg, 120 * deg]

    source.translation = [-19.2 * cm, -13.425 * cm, -22.5 * cm]
    source.half_life = 6 * hour
    # compute the translation to align the source with CT
    # (considering they are in the same physical space)
    # source.position.translation = gate.image.get_translation_between_images_center(
    #     XCAT_phantom.image, source.image
    # )

    print(f"Reading source image {source.image}")

    # # add stat actor
    s = sim.add_actor("SimulationStatisticsActor", "stats")
    s.track_types_flag = True
    s.output = current_path / Path("output") / "stats1.txt"

    # # phys
    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"
    sim.physics_manager.set_production_cut("world", "all", 1 * m)
    sim.physics_manager.set_production_cut("XCAT_phantom", "all", 0.1 * mm)

    # ---------------------------------------------------------------------
    # start simulation
    # sim.running_verbose_level = gate.EVENT
    sim.run_timing_intervals = [[0, 20 * sec]]
    sim.run()

    # print stats
    stats = sim.output.get_actor("stats")
    print(stats)

    # file_number = 1
    # number = 11
    # # print(os.getcwd())

    # SPECT_df = uproot.open(current_path / Path("output") / f"scintigraphy_{i}.root")

    # tree_list = SPECT_df.keys()
    # tree_name = tree_list[0]

    # # 트리에서 데이터 읽기
    # tree = SPECT_df[tree_name]
    # a = tree[0]
    # data = tree.arrays(library="pd")
    # #
    # # # 데이터프레임으로 변환
    # df = pd.DataFrame(data)

    # # # SPECT_df_total = df[['energy', 'globalPosX', 'globalPosY', 'globalPosZ']]
    # # SPECT_df_total = df[['edep', 'posX', 'posY', 'posZ']]
    # SPECT_df_total = df[
    #     ["TotalEnergyDeposit", "PostPosition_X", "PostPosition_Y", "PostPosition_Z"]
    # ]

    # SPECT_df_total = SPECT_df_total.to_numpy()

    # io.savemat(
    #     current_path / Path("output") / f"scintigraphy_{i}.mat",
    #     {"thyroid": SPECT_df_total},
    # )
