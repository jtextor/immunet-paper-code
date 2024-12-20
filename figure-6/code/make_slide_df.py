import numpy as np
import scipy.ndimage as nd
import pandas as pd
from skimage.filters import threshold_local, gaussian
from skimage.morphology import remove_small_objects
from config import data_folder, slides_data_folder, cache_folder, prediction_folder

if not data_folder.is_dir():
    data_folder.mkdir()

if not slides_data_folder.is_dir():
    slides_data_folder.mkdir()

thresholds = {"2020-02-12-phenotyping-paper-lung-bcell": [0, 0],
    "2020-01-31-phenotyping-paper-prostate": [0.1, 0],
    "2020-01-31-phenotyping-paper-melanoma": [0.01, 0.1],
    "2020-01-31-phenotyping-paper-bladder": [0.0001,0.0001]}

datasets = [f.name for f in cache_folder.iterdir() if f.is_dir()]

# Set seed for reproducibility
np.random.seed(99)

for dataset in datasets:
    stroma_threshold, tumor_threshold = thresholds[dataset]

    dataset_folder = cache_folder / dataset
    slide_cache_files = [x for x in dataset_folder.glob("*.npz")]

    for slide_cache_file in slide_cache_files:
        slide = slide_cache_file.stem
        print(dataset, slide)
        slide_prediction_path = prediction_folder / dataset / (slide + ".csv")

        if slide_cache_file.is_file() and slide_prediction_path.is_file():
            # Extract segmentation
            slide_data_path = slides_data_folder / dataset
            if not slide_data_path.is_dir():
                slide_data_path.mkdir()

            cachefile = np.load(slide_cache_file)
            stroma = cachefile["tissue"]
            tumor = cachefile["tumor"]

            sigma = 1
            t = threshold_local(tumor.astype(np.float32), block_size=9, offset=-1)
            mask = tumor.astype(np.float32)>t
            mask = gaussian(mask[:,:], sigma=sigma)>tumor_threshold
            mask = nd.binary_fill_holes(mask)
            tumor_mask = remove_small_objects(mask, min_size=100)

            t = threshold_local(stroma.astype(np.float32), block_size=9, offset=-1)
            mask = stroma.astype(np.float32)>t
            mask = gaussian(mask[:,:], sigma=sigma)>stroma_threshold

            stroma_mask = nd.binary_fill_holes(mask)
            stroma_mask = remove_small_objects(stroma_mask, min_size=100)

            stroma_mask[tumor_mask] = False

            dist_map = nd.distance_transform_edt(tumor_mask==False)
            dist_map[tumor_mask] = -nd.distance_transform_edt(tumor_mask)[tumor_mask]

            tumor_region = np.sum(tumor_mask)
            stroma_region = np.sum(stroma_mask)

            # Extract prediction
            prediction_df = pd.read_csv(slide_prediction_path)
            prediction_df = prediction_df.reset_index()

            slide_prediction = []
            for index, row in prediction_df.iterrows():
                cell_data = row.to_dict()
                del cell_data["index"]

                tumor_expression = t[int(cell_data["y"] / 8), int(cell_data["x"] / 8)]
                cell_data["tumor"] = tumor_expression

                distance = dist_map[int(cell_data["y"] / 8), int(cell_data["x"] / 8)]

                random_distance = float("Inf")
                inside_stroma = False
                inside_tumor = False
                while not (inside_tumor or inside_stroma):
                    a = np.random.randint(dist_map.shape[0])
                    b = np.random.randint(dist_map.shape[1])
                    inside_tumor = tumor_mask[a, b]
                    inside_stroma = stroma_mask[a, b]
                    random_distance = dist_map[a, b]

                cell_data["distance"] = distance * 4
                cell_data["random_distance"] = random_distance * 4
                cell_data["tumor_area"] = tumor_region
                cell_data["stroma_area"] = stroma_region
                slide_prediction.append(cell_data)

            slide_prediction_path = slide_data_path / (slide + ".csv")
            df = pd.DataFrame(slide_prediction)
            df.to_csv(slide_prediction_path, index=False)
        else:
            print("Cache or prediction is missing for dataset {}, slide {}".format(dataset, slide))
