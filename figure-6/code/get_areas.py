import numpy as np
import scipy.ndimage as nd
import pandas as pd
from skimage.filters import threshold_local, gaussian
from skimage.morphology import remove_small_objects
from config import data_folder, cache_folder

if not data_folder.is_dir():
    data_folder.mkdir()

thresholds = {"2020-02-12-phenotyping-paper-lung-bcell": [0, 0],
    "2020-01-31-phenotyping-paper-prostate": [0.1, 0],
    "2020-01-31-phenotyping-paper-melanoma": [0.01, 0.1],
    "2020-01-31-phenotyping-paper-bladder": [0.0001,0.0001]}

data = []
datasets = [f.name for f in cache_folder.iterdir() if f.is_dir()]

for dataset in datasets:
    stroma_threshold, tumor_threshold = thresholds[dataset]
    dataset_folder = cache_folder / dataset
    slide_cache_files = [x for x in dataset_folder.glob("*.npz")]

    for slide_cache_file in slide_cache_files:
        slide = slide_cache_file.stem
        print(dataset, slide)

        if slide_cache_file.is_file():
            cachefile = np.load(slide_cache_file)
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

            stroma_mask[tumor_mask] = False
            stroma_mask = remove_small_objects(stroma_mask, min_size=400)

            dist_map = nd.distance_transform_edt(tumor_mask == False)
            dist_map[stroma_mask == False] = float("Inf")
            dist_map[tumor_mask] = -nd.distance_transform_edt(stroma_mask == False)[tumor_mask]

            tumor_region = np.sum(tumor_mask) * 4 ** 2
            stroma_region = np.sum(stroma_mask) * 4 ** 2
            data.append({"dataset": dataset, "slide": slide, "distance": float("Inf"), "tumor": tumor_region,
                         "stroma": stroma_region})

            for distance in [100, 500, 1000]:
                tumor_region_dist = np.sum(np.logical_and(dist_map < 0, dist_map > -distance / 4., tumor_mask)) * 4 ** 2
                stroma_region_dist = np.sum(
                    np.logical_and(dist_map >= 0, dist_map < distance / 4., stroma_mask)) * 4 ** 2
                print(distance, tumor_region_dist, stroma_region_dist)
                data.append({"dataset": dataset, "slide": slide, "distance": distance, "tumor": tumor_region_dist,
                             "stroma": stroma_region_dist})
        else:
            print("Cache is missing for dataset {}, slide {}".format(dataset, slide))


df = pd.DataFrame(data)
df.to_csv(data_folder / "regions_per_distance.csv", index=False)
