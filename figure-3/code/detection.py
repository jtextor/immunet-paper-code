import numpy as np
from tqdm import tqdm
from scipy.spatial import cKDTree
import tarfile
from utils import load_annotations, PredictionHandler, ImmuNetPredictionLoader, InFormPredictionLoader, data_path, \
    immunet_pred_path, inform_pred_path, images_path, model_path, pred_archive_path, pred_path
# ImmuNet: uncomment to use ImmuNet to find cells on tiles
#from inference import ImmuNetPredictionHandler

DS_TS = {"2020-01-31-phenotyping-paper-bladder": "bladder",
         "2020-01-27-phenotyping-paper-cytoagars": "cytoagars",
         "2020-02-12-phenotyping-paper-lung-bcell": "lung",
         "2020-01-31-phenotyping-paper-melanoma": "melanoma",
         "2020-01-31-phenotyping-paper-prostate": "prostate",
         "2020-01-27-phenotyping-paper-tonsils": "tonsils"}


# Tiles in training set lost in new inForm analysis
tile_ids_to_exclude = ["2020-01-31-phenotyping-paper-prostate/prostate06/55154,3953",
                       "2020-01-31-phenotyping-paper-prostate/prostate06/58499,4953",
                       "2020-01-31-phenotyping-paper-prostate/prostate03/61592,8426"]


def match_cells(annotations_path, prediction_loader: PredictionHandler, fout, prediction_handler=None):
    # Matches annotations with predictions and saves the results to a file

    # extract dictionary {tile_id: annotations}
    annotations_dict = load_annotations(annotations_path)

    with open(fout, "w+") as f:
        f.write("\t".join(prediction_loader.result_columns()))
        f.write("\n")
        for (tile_id, annotations) in tqdm(annotations_dict.items()):
            if tile_id in tile_ids_to_exclude:
                continue

            dataset, slide, tile = tile_id.split("/")
            tissue = DS_TS[dataset]

            if prediction_handler is not None and prediction_handler.tile_exists(dataset, slide, tile):
                prediction = prediction_handler.get_prediction(dataset, slide, tile)
            else:
                prediction = prediction_loader.get_prediction(dataset, slide, tile)

            if prediction is None:
                continue

            # If no cells is predicted
            if len(prediction) == 0:
                prediction = [[[-100000, -100000], [0, 0, 0, 0, 0]]]

            tr = cKDTree(np.array([cell[0] for cell in prediction]))
            for ann in annotations:
                f.write("\t".join((tissue, ann["id"], ann["type"])))
                nn = tr.query([ann["y"], ann["x"]])
                f.write("\t%.1f\t" % (nn[0]))
                f.write(prediction_loader.phenotype_output(prediction[nn[1]][1]))
                f.write("\t")
                positivity = (
                    ann["positivity"] if "positivity" in ann else [0, 0, 0, 0, 0, 0]
                )
                f.write("\t".join(["%g" % i for i in positivity]))
                f.write("\n")
                f.flush()


def match_cells_immunet(annotations_path, prediction_path, suffix=None):
    if suffix is None:
        fout = data_path / "prediction-immunet.tsv"
    else:
        fout = data_path / f"prediction-immunet-{suffix}.tsv"

    prediction_loader = ImmuNetPredictionLoader(prediction_path)
    match_cells(annotations_path, prediction_loader, fout)

    # ImmuNet: replace two lines above with these lines to use ImmuNet to find cells on tiles
    #prediction_loader = ImmuNetPredictionLoader(prediction_path)
    #prediction_handler = ImmuNetPredictionHandler(model_path, images_path)
    #match_cells(annotations_path, prediction_loader, fout, prediction_handler)


def match_cells_inform(annotations_path, prediction_path, suffix=None):
    if suffix is None:
        fout = data_path / "prediction-inform.tsv"
    else:
        fout = data_path / f"prediction-inform-{suffix}.tsv"

    prediction_loader = InFormPredictionLoader(prediction_path)
    match_cells(annotations_path, prediction_loader, fout)


if __name__ == "__main__":
    if not data_path.is_dir():
        data_path.mkdir()

    train_ann_path = data_path / "annotations_train.json.gz"
    val_ann_path = data_path / "annotations_val.json.gz"

    if not pred_path.is_dir():
        file = tarfile.open(pred_archive_path)
        # extracting file
        file.extractall(data_path)
        file.close()

    # Obtain prediction and match with annotations for training data
    match_cells_immunet(train_ann_path, immunet_pred_path, "train")
    # Obtain prediction and match with annotations for validation data
    match_cells_immunet(val_ann_path, immunet_pred_path, "val")

    # inForm
    # Obtain prediction and match with annotations for training data
    match_cells_inform(train_ann_path, inform_pred_path, "train")
    # Obtain prediction and match with annotations for validation data
    match_cells_inform(val_ann_path, inform_pred_path, "val")
