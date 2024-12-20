from abc import ABCMeta, abstractmethod
import json
import gzip
from pathlib import Path
import pandas as pd

data_path = Path("data")
pred_archive_path = data_path / "prediction.tar.gz"
pred_path = data_path / "prediction"
immunet_pred_path = pred_path / "immunet"
inform_pred_path = pred_path / "inform"
panels_file_path = data_path / "panels.json"

# Paths to use when running inference with a model
model_path = data_path/ "immunet.h5"
images_path = data_path / "tilecache"

def load_annotations(annotations_path):
    annotations_dict = {}
    with gzip.open(annotations_path) as f:
        annotations = json.loads(f.read())

        for dataset in annotations:
            for slide in dataset["slides"]:
                for tile in slide["tiles"]:
                    ttile = f"{dataset['ds']}/{slide['slide']}/{tile['tile']}"
                    annotations_dict[ttile] = tile["annotations"]

    return annotations_dict


class PredictionHandler(metaclass=ABCMeta):
    @abstractmethod
    def get_prediction(self, dataset, slide, tile):
        raise NotImplementedError

    @abstractmethod
    def phenotype_output(self, phenotype):
        raise NotImplementedError

    @abstractmethod
    def result_columns(self):
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

    def phenotype_output(self, phenotype):
        return "\t".join(["%.2f" % i for i in phenotype])

    def result_columns(self):
        return ["tissue", "id", "ann_type", "dist", "pCD3", "pFOXP3", "pCD20", "pCD45RO", "pCD8", "CD3", "FOXP3", "CD20", "CD45RO", "CD8"]


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

    def phenotype_output(self, phenotype):
        return "{}".format(phenotype)

    def result_columns(self):
        return ["tissue", "id", "ann_type", "dist", "pred_type", "CD3", "FOXP3", "CD20", "CD45RO", "CD8"]
