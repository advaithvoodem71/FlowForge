import os
import subprocess
import pandas as pd
import openvsp as vsp

# PARAMETERS
winglet_heights = [1.0, 1.5, 2.0]  # meters
aoa_values = [0, 5, 10, 15]

results = []

for h in winglet_heights:
    # GEOMETRY
    vsp.ClearVSPModel()
    wing_id = vsp.AddGeom("WING")
    vsp.SetParmVal(wing_id, "Span", "XSec_1", 35.0)  # 777 half-span approx
    vsp.SetParmVal(wing_id, "Aspect", "XSec_1", 10.0)

    # Add winglet
    winglet_id = vsp.AddGeom("WING", wing_id)
    vsp.SetParmVal(winglet_id, "Span", "XSec_1", h)
    vsp.SetParmVal(winglet_id, "Sweep", "XSec_1", 15.0)
    vsp.SetParmVal(winglet_id, "Dihedral", "XSec_1", 60.0)

    # Export STL
    fname = f"winglet_h{h}.stl"
    vsp.ExportFile(fname, vsp.SET_ALL, vsp.EXPORT_STL)

    # MESHING
    mesh_out = f"winglet_h{h}.su2"
    gmsh_cmd = f"gmsh {fname} -3 -format su2 -o {mesh_out}"
    subprocess.run(gmsh_cmd, shell=True)

    # RUN SU2
    for aoa in aoa_values:
        cfg_file = f"case_h{h}_aoa{aoa}.cfg"
        with open(cfg_file, "w") as f:
            f.write(f"MESH_FILENAME= {mesh_out}\n")
            f.write("SOLVER= RANS\n")
            f.write("TURB_MODEL= SST\n")
            f.write("MACH_NUMBER= 0.85\n")
            f.write("REYNOLDS_NUMBER= 1e7\n")
            f.write(f"AOA= {aoa}\n")
            f.write("ITER= 2000\n")
            f.write("CONV_RESIDUAL_MINVAL= -8\n")
            f.write(f"SURFACE_OUTPUT= forces_h{h}_aoa{aoa}.dat\n")

        subprocess.run(f"SU2_CFD {cfg_file}", shell=True)

        # EXTRACT RESULTS
        with open(f"forces_h{h}_aoa{aoa}.dat") as ff:
            lines = ff.readlines()
            cl = float(lines[-1].split()[1])
            cd = float(lines[-1].split()[2])
            ld = cl/cd
            results.append({"h": h, "aoa": aoa, "CL": cl, "CD": cd, "L/D": ld})

# Save dataset
df = pd.DataFrame(results)
df.to_csv("winglet_dataset.csv", index=False)
print(df)