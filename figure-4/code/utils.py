from pathlib import Path
import pandas as pd
from abc import ABCMeta, abstractmethod

data_path = Path("data")
pred_path = data_path / "prediction"
immunet_pred_path = pred_path / "immunet"
inform_pred_path = pred_path / "inform"
rois_json_path = data_path / "val_rois.json.gz"
panels_file_path = data_path / "panels.json"

# Paths to use when running inference with a model
model_path = data_path/ "immunet.h5"
images_path = data_path / "tilecache"


DS_TS = {"2020-01-31-phenotyping-paper-bladder": "bladder",
         "2020-01-27-phenotyping-paper-cytoagars": "cytoagars",
         "2020-02-12-phenotyping-paper-lung-bcell": "lung",
         "2020-01-31-phenotyping-paper-melanoma": "melanoma",
         "2020-01-31-phenotyping-paper-prostate": "prostate",
         "2020-01-27-phenotyping-paper-tonsils": "tonsils"}

class PredictionHandler(metaclass=ABCMeta):
    @abstractmethod
    def get_prediction(self, dataset, slide, tile):
        raise NotImplementedError


class PredictionLoader(PredictionHandler):
    def __init__(self, prediction_path):
        self.prediction_path = prediction_path


class ImmuNetPredictionLoader(PredictionLoader):
    def get_prediction(self, dataset, slide, tile):
        tile_prediction_path = self.prediction_path / dataset / slide / tile / "prediction.csv"

        if not tile_prediction_path.is_file():
            return None

        try:
            df = pd.read_csv(tile_prediction_path)
            df = df.reset_index()
        except:
            print("ImmuNet: problem with a tile {} {} {}".format(dataset, slide, tile))
            return None

        prediction = []
        for index, row in df.iterrows():
            pred_pheno = row["phenotype"]
            pred_pheno = pred_pheno.strip("'[]").split(", ")
            pred_pheno = [float(x) for x in pred_pheno]
            prediction.append([[row['y'], row['x']], pred_pheno])

        return prediction


class InFormPredictionLoader(PredictionLoader):
    def get_prediction(self, dataset, slide, tile):
        tile_prediction_path = self.prediction_path / dataset / slide / tile / "prediction.csv"

        if not tile_prediction_path.is_file():
            return None

        try:
            df = pd.read_csv(tile_prediction_path)
            df = df.reset_index()
        except:
            print("Inform: problem with a tile {} {} {}".format(dataset, slide, tile))
            return None

        prediction = [[[row['y'], row['x']], row["phenotype"]] for _, row in df.iterrows()]
        return prediction
