#!/usr/bin/env python3                                                                                        
# -*- coding: utf-8 -*-                                                                                       
import opengate as gate                                                                                       
import pathlib                                                                                                
from pathlib import Path                                                                                      
from scipy.spatial.transform import Rotation                                                                  
from opengate.geometry.utility import get_grid_repetition, get_circular_repetition                            
import opengate.geometry.materials as gmat                                                                    
import numpy as np                                                                                            
import uproot
import pandas as pd
from scipy import io
                                                                                                              
current_path = pathlib.Path(__file__).parent.resolve()                                                        
output_path = current_path / "output"                                                                         
                                                                                                              
# colors                                                                                                      
invisible = [0, 0, 0, 0]                                                                                      
red = [1, 0, 0, 1]                                                                                            
blue = [0, 0, 1, 1]                                                                                           
green = [0, 1, 0, 1]                                                                                          
                                                                                                              
yellow = [0.9, 0.9, 0.3, 1]                                                                                   
gray = [0.5, 0.5, 0.5, 1]                                                                                     
white = [1, 1, 1, 0.8]                                                                                        
purple = [0.502, 0, 0.502, 1]                                                                                 
                                                                                                              
                                                                                                              
if __name__ == "__main__":                                                                                    
    # create the simulation                                                                                   
    sim = gate.Simulation()                                                                                   
                                                                                                              
    # main options                                                                                            
    sim.g4_verbose = False                                                                                    
    sim.visu = False                                                                                          
    # sim.visu_type = "gdml"                                                                                  
    sim.visu_type = "vrml"                                                                                    
    # sim.visu_type = "qt"                                                                                    
    sim.number_of_threads = 25                                                                                
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
                                                                                                              
    # Detector                                                                                                
    SPECT_outer = sim.add_volume("Box", "SPECT_outer")                                                        
    SPECT_outer.mother = world.name                                                                           
    SPECT_outer.size = [15 * cm, 4 * cm, 4 * cm]                                                              
    SPECT_outer.rotation = None                                                                               
    SPECT_outer.material = "Lead"                                                                             
    SPECT_outer.color = green                                                                                 
    SPECT_outer.translation = [-5.45 * cm, 5 * cm, 0 * cm]                                                   
                                                                                                              
    SPECT_inner = sim.add_volume("Box", "SPECT_inner")                                                        
    SPECT_inner.mother = SPECT_outer.name                                                                     
    SPECT_inner.size = [14.8 * cm, 3.8 * cm, 3.8 * cm]                                                        
    SPECT_inner.rotation = None                                                                               
    SPECT_inner.material = "Air"                                                                              
    SPECT_inner.color = yellow                                                                                
    SPECT_inner.translation = [0 * cm, 0 * cm, 0 * cm]                                                        
                                                                                                              
    Detector_head = sim.add_volume("Box", "Detector_head")                                                    
    Detector_head.mother = SPECT_outer.name                                                                   
    Detector_head.size = [2.2 * cm, 0.2 * cm, 2.2 * cm]                                                       
    Detector_head.rotation = None                                                                             
    Detector_head.material = "Air"                                                                            
    Detector_head.color = white                                                                               
    Detector_head.translation = [5.45 * cm, -1.9 * cm, 0 * cm]                                                
                                                                                                              
    Detector = sim.add_volume("Box", "Detector")                                                              
    Detector.mother = Detector_head.name                                                                      
    Detector.size = [2 * cm, 0.1 * cm, 2 * cm]                                                                
    Detector.rotation = None                                                                                  
    Detector.material = "CdTe"                                                                                
    Detector.color = red                                                                                      
    Detector.translation = [0 * cm, 0.05 * cm, 0 * cm]                                                        
                                                                                                              
    # pixel = sim.add_volume("Box", "pixel")                                                                  
    # pixel.mother = Detector.name                                                                            
    # pixel.size = [2 * mm, 0.1 * cm, 2 * mm]                                                                 
    # translations_grid = get_grid_repetition(size=[10, 1, 10], spacing=pixel.size)                           
    # pixel.translation = translations_grid                                                                   
                                                                                                              
    # Back-Compartment                                                                                        
    compartment = sim.add_volume("Box", "compartment")                                                        
    compartment.mother = SPECT_outer.name                                                                     
    compartment.size = [2 * cm, 0.25 * cm, 2 * cm]                                                            
    compartment.rotation = None                                                                               
    compartment.material = "Glass"                                                                            
    compartment.color = blue                                                                                  
    compartment.translation = [5.45 * cm, -1.675 * cm, 0 * cm]                                                
                                                                                                              
    # spect digitizer channels                                                                                
    channels = [                                                                                              
        {"name": f"spectrum", "min": 100 * keV, "max": 140 * keV},                                            
    ]                                                                                                         
                                                                                                              
    # spect digitizer : Hits + Adder + EneWin + Projection                                                    
    # Hits                                                                                                    
    hc = sim.add_actor("DigitizerHitsCollectionActor", f"Hits_{Detector.name}")                               
    hc.mother = Detector.name                                                                                 
    hc.output = current_path / "output" / "point_1.root"
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
    cc = sim.add_actor("DigitizerEnergyWindowsActor", f"EnergyWindows_{Detector.name}")                       
    cc.mother = sc.mother                                                                                     
    cc.input_digi_collection = sc.name                                                                        
    cc.channels = channels                                                                                    
    cc.output = hc.output                                                                                     
                                                                                                              
    # projection image                                                                                        
    # proj = sim.add_actor("DigitizerProjectionActor", f"Projection_{Detector.name}")                         
    # proj.mother = cc.mother                                                                                 
    # proj.input_digi_collections = [x["name"] for x in cc.channels]                                          
    # proj.spacing = [2 * mm, 2 * mm]                                                                         
    # proj.size = [10, 10]                                                                                    
    # proj.output = current_path / Path("output") / "projection1.mhd"                                         
    # proj.detector_orientation_matrix = Rotation.from_euler("x", 90).as_matrix()                             
                                                                                                              
    # Activity source from an image                                                                           
    source = sim.add_source("GenericSource", "Co57")                                                          
    source.particle = "gamma"                                                                                 
    if sim.visu == True:                                                                                      
        source.activity = 1 * Bq / sim.number_of_threads                                                      
    else:                                                                                                     
        source.activity = 7.5 * 3.7e4 * Bq / sim.number_of_threads                                            
                                                                                                              
                                                                                                              
    translations = [-1.0 * cm, 0 * cm, -1.0 * cm]
                                                                                                              
    source.position.translation = translations                                                                
                                                                                                              
    source.direction.type = "iso"                                                                             
    source.energy.type = "gauss"                                                                            
    source.energy.mono = 122 * keV                                                                            
    source.energy.sigma_gauss = 12.2 * keV                                                                  
                                                                                                              
    # # add stat actor                                                                                        
    s = sim.add_actor("SimulationStatisticsActor", "stats")                                                   
    s.track_types_flag = True                                                                                 
    s.output = current_path / Path("output") / "stats1.txt"                                                   
                                                                                                              
    # # phys                                                                                                  
    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option4"                                     
    sim.physics_manager.set_production_cut("world", "all", 1 * m)                                             
    sim.physics_manager.set_production_cut("Detector", "all", 0.1 * mm)                                       
                                                                                                              
    # ---------------------------------------------------------------------                                   
    # start simulation                                                                                        
    # sim.running_verbose_level = gate.EVENT                                                                  
    sim.run_timing_intervals = [[0, 300 * sec]]                                                               
    sim.run()                                                                                                 
                                                                                                              
    # print stats                                                                                             
    stats = sim.output.get_actor("stats")                                                                     
    print(stats)                                                                                              
current_path = pathlib.Path(__file__).parent.resolve()
SPECT_df = uproot.open(current_path / Path("output") / "point_1.root")
tree_list = SPECT_df.keys()
tree_name = tree_list[2]
tree = SPECT_df[tree_name]
data = tree.arrays(library="pd")
df = pd.DataFrame(data)
SPECT_df_total = df[
    ["TotalEnergyDeposit", "PostPosition_X", "PostPosition_Y", "PostPosition_Z"]
]
SPECT_df_total = SPECT_df_total.to_numpy()
io.savemat(current_path / Path("output") / "point_1.mat", {"point": SPECT_df_total})
