import pandas as pd
import tarfile
# ImmuNet: uncomment to use ImmuNet to find cells on tiles
#from inference import ImmuNetPredictionHandler
from rois import load_val_rois
from utils import DS_TS, ImmuNetPredictionLoader, InFormPredictionLoader, pred_path, pred_archive_path, immunet_pred_path, inform_pred_path, \
    rois_json_path, data_path, model_path, images_path


def analyse():
    results = []
    rois = load_val_rois(rois_json_path)

    immunet_prediction_loader = ImmuNetPredictionLoader(immunet_pred_path)
    inform_prediction_loader = InFormPredictionLoader(inform_pred_path)

    # ImmuNet: uncomment to use ImmuNet to find cells on tiles
    #prediction_handler = ImmuNetPredictionHandler(model_path, images_path)

    rois = [roi for roi in rois if roi.dataset != "2020-01-27-phenotyping-paper-cytoagars"]

    for roi in rois:
        print("Dataset {}, slide {}, tile{}". format(roi.dataset, roi.slide, roi.tile))

        if roi.slide == "melanoma20" and roi.tile == "58334,7066":
            continue

        if roi.slide == "prostate01" and roi.tile == "59982,17214":
            continue

        # Load immunet predictions
        roi.load_immunet_cells(immunet_prediction_loader)
        # ImmuNet: replace the line above with the following line to use ImmuNet to find cells on tiles
        #roi.load_immunet_cells(immunet_prediction_loader, prediction_handler)

        # Load inForm predictions
        roi.load_inform_cells(inform_prediction_loader)

        immunet_row = { "roi_id": roi.id,
                        "annotations": len(roi.ann_cells), "detections": len(roi.pred_cells_immunet),
                        "dataset": DS_TS[roi.dataset], "area": roi.area, "algorithm": "NN"}
        results.append(immunet_row)

        inform_row = {"roi_id": roi.id,
                       "annotations": len(roi.ann_cells), "detections": len(roi.pred_cells_inform),
                       "dataset": DS_TS[roi.dataset], "area": roi.area, "algorithm": "inform"}
        results.append(inform_row)

    return pd.DataFrame(results)


if not pred_path.is_dir():
    file = tarfile.open(pred_archive_path)
    # extracting file
    file.extractall(data_path)
    file.close()

df = analyse()
df.to_csv(data_path / "counts.csv")
