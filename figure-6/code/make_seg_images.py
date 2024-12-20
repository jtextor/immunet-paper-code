import numpy as np
import scipy.ndimage as nd
from skimage.filters import threshold_local, gaussian
from skimage.morphology import remove_small_objects
from PIL import Image
import pandas as pd
from config import data_folder, cache_folder, images_folder

if not data_folder.is_dir():
    data_folder.mkdir()

if not images_folder.is_dir():
    images_folder.mkdir()


tumor_color = [158 / 255., 202 / 255., 225 / 255.]
tumor_border_color = [66 / 255., 146 / 255., 198 / 255.]
stroma_color = [0.9, 0.9, 0.9]

#stroma, tumor thresholds
thresholds = {"2020-02-12-phenotyping-paper-lung-bcell": [0, 0],
    "2020-01-31-phenotyping-paper-prostate": [0.1, 0],
    "2020-01-31-phenotyping-paper-melanoma": [0.01, 0.1],
    "2020-01-31-phenotyping-paper-bladder": [0.0001, 0.0001]}

df = pd.read_csv(data_folder / "low_high.csv")

for index, row in df.iterrows():
    dataset = row['dataset']
    slide = row['slide']

    stroma_threshold, tumor_threshold = thresholds[dataset]

    if (dataset == "2020-01-31-phenotyping-paper-bladder") and (slide == "bladder12"):
        raise ValueError("Cannot handle this dataset and slide, change those used in low_high.csv")
    if slide == "melanoma06" or slide == "melanoma11":
        raise ValueError("Cannot handle this dataset and slide, change those used in low_high.csv")

    slide_cache_path = cache_folder / dataset / (slide + ".npz")
    if slide_cache_path.is_file():
        cachefile = np.load(slide_cache_path)
        stroma = cachefile["tissue"]
        tumor = cachefile["tumor"]

        sigma = 1
        t = threshold_local(tumor.astype(np.float32), block_size=9, offset=-1)
        mask = tumor.astype(np.float32) > t
        mask = gaussian(mask[:, :], sigma=sigma) > tumor_threshold
        mask = nd.binary_fill_holes(mask)
        tumor_mask = remove_small_objects(mask, min_size=150)

        t = threshold_local(stroma.astype(np.float32), block_size=9, offset=-1)
        mask = stroma.astype(np.float32) > t
        mask = gaussian(mask[:, :], sigma=sigma) > stroma_threshold

        stroma_mask = nd.binary_fill_holes(mask)

        print(np.sum(stroma_mask))
        stroma_mask[tumor_mask] = False
        stroma_mask = remove_small_objects(stroma_mask, min_size=400)
        print(np.sum(stroma_mask))

        tumor_dist_map = nd.distance_transform_edt(tumor_mask == False)
        tumor_dist_map[stroma_mask == False] = float("Inf")
        tumor_dist_map[tumor_mask] = -nd.distance_transform_edt(stroma_mask == False)[tumor_mask]

        # Make image of tumor landscape without borders
        a = np.ones(tumor_mask.shape + (3,))
        a[tumor_mask] = tumor_color
        a[stroma_mask] = stroma_color
        im = Image.fromarray((a * 255).astype(np.uint8))
        im.save(images_folder / f"{dataset}.{slide}.segmentation.png")

        # Make image of tumor border
        a = np.ones(tumor_mask.shape + (4,))
        a[:, :, :] = [1, 1, 1, 0]
        tumor_border = (tumor_dist_map < 1.1) & (tumor_dist_map > -1.1)
        a[tumor_border] = (66 / 255., 146 / 255., 198 / 255., 1)
        im = Image.fromarray((a * 255).astype(np.uint8))
        im.save(images_folder / f"{dataset}.{slide}.border.png")
