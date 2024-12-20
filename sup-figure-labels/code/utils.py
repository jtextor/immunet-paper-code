import numpy as np
from scipy.ndimage import distance_transform_edt
from pathlib import Path

DATA_FOLDER = Path("data")
IMAGE_FOLDER = Path("images")
ANNOTATIONS_FILE = DATA_FOLDER / "annotations.json.gz"
COMPOSITE_FILE = IMAGE_FOLDER / "composite.png"

def rnd(i):
    return int(round(i))


def background(c):
    return c["t"] in ["Tumor cell", "Other cell", "No cell"]


def map_phenotype(c, channels_num):
    if "positivity" not in c:
        return [0, 0, 0, 0, 0]

    positivity = c["positivity"][1:6]
    if len(positivity) < channels_num:
        positivity += list(np.ones(channels_num - len(positivity)).astype(int))

    return (np.array(positivity) - 1) / 4.0


def extract_labels(components, annotations, out_markers_num=5, cell_radius=5):
    h = components.shape[0]
    w = components.shape[1]

    # The first output channel is a prediction of a distance map,
    # the rest are predictions of phenotype marker expression maps
    out_channels_num = 1 + out_markers_num
    out = np.zeros((h, w, out_channels_num), np.float16)
    known_status = np.zeros((h, w), np.uint8)

    # 0 - no cell, i - cell number (just an order in which we got annotations)
    cell_at = np.zeros((h, w), np.uint32)
    if out_markers_num > 0:
        cell_ph = np.zeros((len(annotations), out_markers_num), np.float16)

    cell_known = np.zeros(len(annotations))
    cell_fg = np.zeros(len(annotations))
    for i, a in enumerate(annotations):
        x, y = map(rnd, [a["x"], a["y"]])
        if x >= w:
            x = w - 1
        if y >= h:
            y = h - 1
        cell_at[y, x] = i + 1
        if out_markers_num > 0:
            cell_ph[i, :] = map_phenotype(a, out_markers_num)
        cell_fg[i] = not background(a)
        cell_known[i] = True

    # distance to the nearest cell (will be false in cell_at == 0 matrix) and indises
    # cell_at == 0 means 1 at background, 0 - foreground
    # od - euclidean distance to the cell center (foreground points)
    # oi - indices of the closest cell center
    od, oi = distance_transform_edt(cell_at == 0, return_indices=True)
    in_cell = od <= cell_radius

    # for each x, y get cell index in order of cell_at (how we got annotations)
    which_cell = cell_at[oi[0][in_cell], oi[1][in_cell]] - 1
    # Sets 255 for all circles that represent cells, 0 otherwise
    known_status[in_cell] = cell_known[which_cell] * 255

    # Distances channel
    # for places inside cells set cell_radius - distance to cell so that max value were at the center of the cell
    # and -2 for background
    out[in_cell, 0] = np.where(cell_fg[which_cell], cell_radius - od[in_cell], -2)
    # Phenotypes channel - same ph value for all points inside the cell
    if out_markers_num > 0:
        out[in_cell, 1:] = cell_ph[which_cell, :]

    # places without annotations have unknown status
    out[known_status == 0, 0] = -1
    return out
