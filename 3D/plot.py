import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
import time as tp
import os
import subprocess
from scipy.io import FortranFile

# LICENSE FILE is included in the upper folder
# this file processes the quasi-magnetostatic 3D results

###############################################################################
#### SETUP ####################################################################

inputs_filename='inputs.txt'            # add inputs file name for loading info below

Fields = ['Hx', 'Hy', 'Hz', 'Hmag']     # fields to plot in a heat map, user can remove if any are unwanted

Hx_monitor=[40,40]                      # location to monitor and produce an fft for Hx fields
Hy_monitor=[40,40]                      # location to monitor and produce an fft for Hy fields
Hz_monitor=[40,40]                      # location to monitor and produce an fft for Hz fields

tsteps_to_ignore_for_fft=0              # number of time steps to ignore before fft if desired

scale = 'auto'                          # 'auto' or 'manual'
max_vals = [1, 1, 1 ,1]                 # if manual, these are upper bounds
min_vals = [0, 0, 0, 0]                 # if manual, these are lower bounds

log_mag = False                         # use log scale for magnitude - can't use with turn off sci if True
turn_off_sci = True                     # Use scientific notation or not

Normalize_mu  = True                    # Produces B field instead in µT - assumes mu of free space everywhere though - not useful for permeability materials

color_mapping = 'nipy_spectral'         # color mapping for the fields

use_geometry = False                    # show geometry of cross section or not
opacity = 0.3                           # opacity for the combination: 0.0 shows no EM fields, 1.0 shows no geometry
geom_colors = "Pastel1"                 # color mapping for the geometry
colors = ["white", "white", "black"]    # optional method for geometry colors - first entry will be background.
geom_colors = LinearSegmentedColormap.from_list("custom_div", colors)

###############################################################################
###############################################################################

###############################################################################
# Determine output writer: prefer ffmpeg (MP4), fall back to pillow (GIF)
def get_writer(still):
    if still:
        return None, '.png'
    if animation.FFMpegWriter.isAvailable():
        return animation.FFMpegWriter(fps=30), '.mp4'
    return animation.PillowWriter(fps=30), '.gif'

#setup ffmpeg settings
def make_ffmpeg_proc(out_path, fig, fps=30):
    w = int(fig.get_figwidth()  * fig.dpi)
    h = int(fig.get_figheight() * fig.dpi)
    # round to even (ffmpeg requirement)
    w += w % 2
    h += h % 2
    cmd = [
    'ffmpeg', '-y',
    '-f', 'rawvideo',   # <-- Fixed the double hyphen here!
    '-vcodec', 'rawvideo',
    '-s', f'{w}x{h}',
    '-pix_fmt', 'rgba',
    '-r', str(fps),
    '-i', 'pipe:0',
    '-vcodec', 'libopenh264',
    '-pix_fmt', 'yuv420p',
    out_path
    ]
    return subprocess.Popen(cmd, stdin=subprocess.PIPE), w, h

###############################################################################
start_time = tp.time()

sizing=np.loadtxt(inputs_filename, max_rows=1, delimiter=",", skiprows=0)
del_t=np.loadtxt(inputs_filename, max_rows=1, skiprows=2)
div=np.loadtxt(inputs_filename, max_rows=1, skiprows=4)
slice_axis=[]
with open(inputs_filename, 'r') as f:
    lines = f.readlines()
    raw_line = lines[5].strip().split(',')
    slice_axis = raw_line[0].strip()       # 'z' (str)
    slice_location = int(raw_line[1].strip())   # 25  (int)

x_size=int(sizing[0])
y_size=int(sizing[1])
z_size=int(sizing[2])

if slice_axis=='x':
    slice_type=0
    video2=y_size
    video1=z_size
if slice_axis=='y':
    slice_type=1
    video2=x_size
    video1=z_size
if slice_axis=='z':
    slice_type=2
    video2=x_size
    video1=y_size


###############################################################################
# Axis labels and limits based on slice type
slice_meta = {
    0: ('Cubes in Y direction', 'Cubes in Z direction', 'X',
        (0.5, y_size + 0.5), (0.5, z_size + 0.5)),
    1: ('Cubes in X direction', 'Cubes in Z direction', 'Y',
        (0.5, x_size + 0.5), (0.5, z_size + 0.5)),
    2: ('Cubes in X direction', 'Cubes in Y direction', 'Z',
        (0.5, x_size + 0.5), (0.5, y_size + 0.5)),
}
xlabel, ylabel, zlabel, xlim, ylim = slice_meta[slice_type]

###############################################################################
# Load geometry once (shared across all field components)
geom_rgba = None
if use_geometry:
    geom = np.fromfile('{}'.format('geom.bin'), dtype=np.float64)
    geom = geom.reshape((x_size, y_size, z_size), order='F')
    sl = int(slice_location) - 1
    if slice_type == 0:
        mat    = geom[sl, :, :]
        dims   = (y_size, z_size)
        s1_dir, s2_dir = 'y', 'z'
    elif slice_type == 1:
        mat    = geom[:, sl, :]
        dims   = (x_size, z_size)
        s1_dir, s2_dir = 'x', 'z'
    else:
        mat    = geom[:, :, sl]
        dims   = (x_size, y_size)
        s1_dir, s2_dir = 'x', 'y'

###############################################################################
# Helper function to load raw field data
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
    if field_name == 'Hmag':
        hx_data = read_fortran_records('Hx.bin', (video2, video1))
        hy_data = read_fortran_records('Hy.bin', (video2, video1))
        hz_data = read_fortran_records('Hz.bin', (video2, video1))
        return np.sqrt(hx_data**2 + hy_data**2 + hz_data**2)
    else:
        file_path = '{}.bin'.format(field_name)
        return read_fortran_records(file_path, (video2, video1))

###############################################################################
mu=1
if Normalize_mu==True:
    mu=4*np.pi*1E-7*1E6 #1E6 so we set microT

# Auto-scale pass: load all data first, compute bounds, then animate
if scale == 'auto':
    min_vals = []
    max_vals = []
    for field in Fields:
        d = load_field_data(field)
        min_vals.append(mu*float(np.real(d).min()))
        max_vals.append(mu*float(np.real(d).max()))

###############################################################################
# Main loop: one figure per field component
writer, ext = get_writer(False)

# used for monitor location
Hx_history=[]
Hy_history=[]
Hz_history=[]

for i, field in enumerate(Fields):

    # --- Load data ---
    all_frames = np.real(load_field_data(field))*mu
    vmin, vmax = min_vals[i], max_vals[i]
    n_frames   = all_frames.shape[0]

    # --- location monitoring ---
    if field=='Hx':
        Hx_history.append(all_frames)
    if field=='Hy':
        Hy_history.append(all_frames)
    if field=='Hz':
        Hz_history.append(all_frames)

    # --- Build figure ---
    fig, ax = plt.subplots()

    if use_geometry:
        geom_im = ax.imshow(mat.T,
                            origin='lower',
                            extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
                            aspect='equal',
                            cmap=geom_colors,
                            interpolation='nearest',
                            zorder=0)

    if field == 'Hmag' and log_mag==True:
        log_vmin = vmin if vmin > 0 else 1e-5 
        norm_type = LogNorm(vmin=log_vmin, vmax=vmax)
    else:
        norm_type = Normalize(vmin=vmin, vmax=vmax)

    im = ax.imshow(all_frames[0].T,
                origin='lower',
                extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
                aspect='equal',
                cmap=color_mapping,
                norm=norm_type,
                interpolation='nearest',
                alpha=opacity if use_geometry else 1.0,
                zorder=1)

    cb_label = 'A/m'
    if Normalize_mu==True:
        cb_label = 'µT'

    cb = fig.colorbar(im, ax=ax, label=cb_label)
    if turn_off_sci==True:
        cb.formatter.set_useOffset(False)
        cb.formatter.set_scientific(False)
    cb.update_ticks()

    fig.set_layout_engine('none')

    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.grid(True, color='lightgray', which='both', linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)

    title = 'Amplitude {}'.format(field)
    if Normalize_mu==True:
        if field=='Hx':
            fieldB='Bx'
        if field=='Hy':
            fieldB='By'
        if field=='Hz':
            fieldB='Bz'
        if field=='Hmag':
            fieldB='Bmag'
        title = 'Amplitude {}'.format(fieldB)
    ax.set_title(title, x=0.5, y=1.05, ha='center')

    fig.canvas.draw()
    background = fig.canvas.copy_from_bbox(ax.bbox)

    if ext == '.png':
        fig.savefig('{}{}'.format(field, ext), dpi=150)
        plt.close(fig)
        continue

    proc, fw, fh = make_ffmpeg_proc('{}{}'.format(field, ext), fig, fps=30)

    fig.canvas.draw()
    background = fig.canvas.copy_from_bbox(ax.bbox)

    for frame_idx in range(n_frames):
        fig.canvas.restore_region(background)
        
        # MODIFIED HERE: Added .T so the animated frames update with the rotation
        im.set_data(all_frames[frame_idx].T)
        
        ax.draw_artist(im)
        fig.canvas.blit(ax.bbox)
        buf = fig.canvas.buffer_rgba()
        proc.stdin.write(buf.tobytes())

    proc.stdin.close()
    proc.wait()
    plt.close(fig)

    if i < len(Fields) - 1:
        print('{} done. Starting {}.'.format(field, Fields[i + 1]))

# --- now point monitoring processing

# get the rest of the data I need
t=np.linspace(0,div*del_t*(len(Hx_history[0])-1),len(Hx_history[0]))
t=t[tsteps_to_ignore_for_fft:]
freq=np.fft.fftfreq(len(t),del_t*div)
Hx_fft=np.fft.fft(Hx_history[0][tsteps_to_ignore_for_fft:,Hx_monitor[0],Hx_monitor[1]], norm='forward')
Hy_fft=np.fft.fft(Hy_history[0][tsteps_to_ignore_for_fft:,Hy_monitor[0],Hy_monitor[1]], norm='forward')
Hz_fft=np.fft.fft(Hz_history[0][tsteps_to_ignore_for_fft:,Hz_monitor[0],Hz_monitor[1]], norm='forward')

# prep the plotting
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharex='col')
hx_time_data = Hx_history[0][tsteps_to_ignore_for_fft:, Hx_monitor[0], Hx_monitor[1]]
hy_time_data = Hy_history[0][tsteps_to_ignore_for_fft:, Hy_monitor[0], Hy_monitor[1]]
hz_time_data = Hz_history[0][tsteps_to_ignore_for_fft:, Hz_monitor[0], Hz_monitor[1]]
half_freq_idx = len(freq) // 2
freq_half = freq[:half_freq_idx]
hx_fft_db = 10 * np.log10(np.abs(Hx_fft)[:half_freq_idx])
hy_fft_db = 10 * np.log10(np.abs(Hy_fft)[:half_freq_idx])
hz_fft_db = 10 * np.log10(np.abs(Hz_fft)[:half_freq_idx])

# Hx
axes[0, 0].plot(t, hx_time_data, 'x-', color='blue', markersize=4, linewidth=1.5, label='$H_x$')
axes[0, 0].set_ylabel('Amplitude (A/m)', fontsize=11)
axes[0, 0].set_title(f'$H_x$ History at ({Hx_monitor[0]}, {Hx_monitor[1]})', fontsize=12, fontweight='bold')
axes[0, 0].grid(True, linestyle='--', alpha=0.6)

# Hy
axes[1, 0].plot(t, hy_time_data, 'x-', color='blue', markersize=4, linewidth=1.5, label='$H_y$')
axes[1, 0].set_xlabel('Time (s)', fontsize=11)
axes[1, 0].set_ylabel('Amplitude (A/m)', fontsize=11)
axes[1, 0].set_title(f'$H_y$ History at ({Hy_monitor[0]}, {Hy_monitor[1]})', fontsize=12, fontweight='bold')
axes[1, 0].grid(True, linestyle='--', alpha=0.6)

# Hz
axes[2, 0].plot(t, hz_time_data, 'x-', color='blue', markersize=4, linewidth=1.5, label='$H_z$')
axes[2, 0].set_xlabel('Time (s)', fontsize=11)
axes[2, 0].set_ylabel('Amplitude (A/m)', fontsize=11)
axes[2, 0].set_title(f'$H_z$ History at ({Hz_monitor[0]}, {Hz_monitor[1]})', fontsize=12, fontweight='bold')
axes[2, 0].grid(True, linestyle='--', alpha=0.6)

# Hx fft
axes[0, 1].plot(freq_half, hx_fft_db, 'x-', color='green', markersize=4, linewidth=1.5, label='$H_x$ FFT')
axes[0, 1].set_ylabel('Magnitude (dB)', fontsize=11)
axes[0, 1].set_title(f'$H_x$ FFT at ({Hx_monitor[0]}, {Hx_monitor[1]})', fontsize=12, fontweight='bold')
axes[0, 1].grid(True, linestyle='--', alpha=0.6)

# Hy fft
axes[1, 1].plot(freq_half, hy_fft_db, 'x-', color='green', markersize=4, linewidth=1.5, label='$H_y$ FFT')
axes[1, 1].set_xlabel('Frequency (Hz)', fontsize=11)
axes[1, 1].set_ylabel('Magnitude (dB)', fontsize=11)
axes[1, 1].set_title(f'$H_y$ FFT at ({Hy_monitor[0]}, {Hy_monitor[1]})', fontsize=12, fontweight='bold')
axes[1, 1].grid(True, linestyle='--', alpha=0.6)

# Hz fft
axes[2, 1].plot(freq_half, hz_fft_db, 'x-', color='green', markersize=4, linewidth=1.5, label='$H_z$ FFT')
axes[2, 1].set_xlabel('Frequency (Hz)', fontsize=11)
axes[2, 1].set_ylabel('Magnitude (dB)', fontsize=11)
axes[2, 1].set_title(f'$H_z$ FFT at ({Hy_monitor[0]}, {Hy_monitor[1]})', fontsize=12, fontweight='bold')
axes[2, 1].grid(True, linestyle='--', alpha=0.6)

# finalize and save
plt.tight_layout()
plt.savefig("H_fields.png", dpi=300, bbox_inches='tight')

###############################################################################
end_time = tp.time()
print('Total time: {:.2f} minutes'.format((end_time - start_time) / 60))