import os
import warnings
warnings.filterwarnings("ignore")
import tensorflow as tf
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False
import numpy as np
import tifffile
from csbdeep.utils import normalize
from skimage.feature import blob_log
from skimage import draw
from utils import PredictionHandler


class ImmuNetPredictionHandler(PredictionHandler):
    def __init__(self, model_path, images_path):
        self.images_path = images_path
        self.model = tf.keras.models.load_model(str(model_path), custom_objects={"pool": tf.nn.pool, "pad": tf.pad}, compile=False)

    def get_tile_path(self, dataset, slide, tile):
        return self.images_path / dataset / slide / tile / "components.tiff"

    def tile_exists(self, dataset, slide, tile):
        return self.get_tile_path(dataset, slide, tile).is_file()

    def get_prediction(self, dataset, slide, tile):
        tile_path = self.get_tile_path(dataset, slide, tile)
        return find_cells(tile_path, self.model)

    def phenotype_output(self, phenotype):
        return "\t".join(["%.2f" % i for i in phenotype])

    def result_columns(self):
        return ["tissue", "id", "ann_type", "dist", "pCD3", "pFOXP3", "pCD20", "pCD45RO", "pCD8", "DAPI", "CD3", "FOXP3", "CD20", "CD45RO", "CD8"]



def get_input(tile_path, in_channels_num):
    tile = tifffile.imread(tile_path, key=range(0, in_channels_num))
    return np.moveaxis(normalize(tile), 0, -1)


def scale_image(image, scalar, low_clip=0, up_clip=255):
    image = scalar * image
    image[image < low_clip] = low_clip
    image[image > up_clip] = up_clip
    return image


def predict(input, model, dist_scalar=50):
    try:
        dist_map, pheno_map = model.predict(np.array([input]))
    except tf.errors.ResourceExhaustedError:
        raise ValueError("Insufficient GPU memory to run inference for a single tile. Change the code to perform stitching")

    dist_map = dist_map[0][:, :, 0]
    scale_image(dist_map, dist_scalar)

    pheno_map = pheno_map[0]
    return dist_map, pheno_map


def find_cells(tile_path,
               model,
               log_threshold=0.07,
               min_log_std=3,
               max_log_std=5,
               dist_scalar=50):
    # Determine number of input channels from the model
    in_channels_num = model.layers[0].input_shape[0][3]
    input = get_input(tile_path, in_channels_num)

    dist_map, pheno_map = predict(input, model, dist_scalar)

    cell_centers = [(int(x[0]), int(x[1]))
                    for x in blob_log(dist_map, min_sigma=min_log_std, max_sigma=max_log_std, threshold=log_threshold)]

    out_markers_num = pheno_map.shape[2]
    predicted_cells = []
    for coords in cell_centers:
        rr, cc = draw.circle(coords[0], coords[1], 2, dist_map.shape)
        phenotype = [float(np.mean(pheno_map[rr, cc, j])) for j in range(out_markers_num)]
        predicted_cells.append([coords, phenotype])

    return predicted_cells
