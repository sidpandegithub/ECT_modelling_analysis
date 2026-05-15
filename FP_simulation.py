import os
from simnibs import sim_struct, run_simnibs

# Path to the folder containing all participant subfolders
base_folder = r"C:\Users\z5171263\Downloads\Augusta_nifti\simulations"

participant_folders = [f for f in os.listdir(base_folder) if os.path.isdir(os.path.join(base_folder, f))]

for participant in participant_folders:
    print(f"Running simulation for {participant}")

    participant_path = os.path.join(base_folder, participant)

    clean_name = participant.replace("m2m_", "")
    output_folder = os.path.join(base_folder, f"{clean_name}_FrontoParietal_simulation_results")

    # create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    # Create session object
    S = sim_struct.SESSION()
    S.subpath = participant_path
    S.pathfem = output_folder
    S.map_to_surf = False
    S.open_in_gmsh = False
    S.fields = "eEjJ"
    # Set up tDCS electrodes
    tdcslist = S.add_tdcslist()
    tdcslist.currents = [-1e-3, 1e-3]

    # Cathode
    cathode = tdcslist.add_electrode()
    cathode.channelnr = 1
    cathode.dimensions = [50, 50]
    cathode.shape = 'ellipse'
    cathode.thickness = [1.75, 1]
    cathode.centre = 'C2'

    # Anode
    anode = tdcslist.add_electrode()
    anode.channelnr = 2
    anode.dimensions = [50, 50]
    anode.shape = 'ellipse'
    anode.thickness = [1.75, 1]
    anode.centre = 'Fp2'

    
    # Run the simulation
    run_simnibs(S)

    print(f"Simulation finished for {participant}\n")