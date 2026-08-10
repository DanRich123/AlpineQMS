import numpy as np
from pyevtk.hl import imageToVTK
import os
from scipy.io import FortranFile

# move to the current directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Simulation parameters
Nx = 70
Ny = 70
Nz = 70
time_steps = 1
use_time = 1
geom_filename = "geom.bin"
Hx_field_name = "Hx_all.bin"
Hy_field_name = "Hy_all.bin"
Hz_field_name = "Hz_all.bin"
use_B_instead = False # assumes microT and assume mu of free space everywhere - really only useful for metals

# Output file info
output_folder = "./Output"
field_output = "{}/H_Fields".format(output_folder)
if use_B_instead==True:
    field_output = "{}/B_Fields".format(output_folder)
geometry_output = "{}/Geometry".format(output_folder)

# Output folder for VTK files
os.makedirs(output_folder, exist_ok=True)

# some helper definitions
def read_fortran_records(file_path, shape):
    frames = []
    with FortranFile(file_path, 'r') as f:
        while True:
            try:
                rec = f.read_reals(dtype=np.float64)
            except Exception:
                break
            frames.append(rec.reshape(shape, order='F'))
    return np.array(frames)

def load_field_data(field_name):
    file_path = '{}'.format(field_name)
    return read_fortran_records(file_path, (Nz,Ny,Nx))

# import hx fields and process
hx = load_field_data(Hx_field_name)
hy = load_field_data(Hy_field_name)
hz = load_field_data(Hz_field_name)

hx=hx[use_time-1,:,:,:]
hy=hy[use_time-1,:,:,:]
hz=hz[use_time-1,:,:,:]

mu_0=4*np.pi*1E-7

# make image for paraview of fields as 1 unit
if use_B_instead==False:
    imageToVTK(
        field_output, 
        spacing=(1.0, 1.0, 1.0), 
        origin=(0.0, 0.0, 0.0), 
        cellData={"H_Field": (hx, hy, hz)}
    )
if use_B_instead==True:
    imageToVTK(
        field_output, 
        spacing=(1.0, 1.0, 1.0), 
        origin=(0.0, 0.0, 0.0), 
        cellData={"B_Field": (hx*mu_0*1E6, hy*mu_0*1E6, hz*mu_0*1E6)}
    )
# Load geometry
geo = np.fromfile(geom_filename, dtype=np.float64)
geo = geo.reshape((Nz, Ny, Nx)).transpose(2, 1, 0)

# Export as a single static VTK file
imageToVTK(
    geometry_output,
    spacing=(1.0, 1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    cellData={"Material": geo},
)