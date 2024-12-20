import json
from shapely.geometry import Polygon, Point
from panels import load_panels
import gzip

panel = load_panels()["lymphocyte"]
pix_per_mm = 2

class Lymphocyte:
    def __init__(self, x, y, type, cd45ro):
        self.x = x
        self.y = y
        self.type = type
        self.cd45ro = cd45ro

    @property
    def types(self):
        return [self.type.main, self.type.detailed]


class ROI:
    def __init__(self, iden, dataset, slide, tile, polygon=None, ann_cells=None):
        self.id = iden
        self.dataset = dataset
        self.slide = slide
        self.tile = tile
        self.polygon = polygon
        self.ann_cells = ann_cells
        self.pred_cells_immunet = []
        self.pred_cells_inform = []

    @classmethod
    def from_dict(cls, dictionary):
        ann_cells = []
        for ann in dictionary["annotations"]:
            # Exclude DAPI
            ann_pheno = ann["positivity"][1:]
            ann_type = ann["type"]
            ann_cell_type = panel.annotation_phenotype(ann_type, ann_pheno)
            if ann_cell_type.main in panel.foreground_types:
                cd45ro_ann = (ann_pheno[3] - 1) / 4
                lymphocyte = Lymphocyte(ann["x"], ann["y"], ann_cell_type, cd45ro_ann)
                ann_cells.append(lymphocyte)

        roi = cls(
            dictionary["id"],
            dictionary["dataset"],
            dictionary["slide"],
            dictionary["tile"],
            Polygon(dictionary["vertices"]),
            ann_cells
        )
        return roi

    @property
    def area(self):
        """
        Area in millimeters^2
        """
        return self.polygon.area * 10 ** (-6) / pix_per_mm**2

    @property
    def density(self):
        return len(self.ann_cells) / self.area

    def load_immunet_cells(self, prediction_loader, prediction_handler=None):

        if prediction_handler is not None and prediction_handler.tile_exists(self.dataset, self.slide, self.tile):
            prediction = prediction_handler.get_prediction(self.dataset, self.slide, self.tile)
        else:
            prediction = prediction_loader.get_prediction(self.dataset, self.slide, self.tile)

        for cell in prediction:
            pred_pos = Point(cell[0][1], cell[0][0])
            if self.polygon.contains(pred_pos):
                pred_pheno = cell[1]
                pred_cell_type = panel.prediction_phenotype(0, pred_pheno)
                if pred_cell_type.main in panel.foreground_types:
                    cd45ro_pred = pred_pheno[3]
                    lymphocyte = Lymphocyte(pred_pos.x, pred_pos.y, pred_cell_type, cd45ro_pred)
                    self.pred_cells_immunet.append(lymphocyte)


    def load_inform_cells(self, prediction_loader):
        prediction = prediction_loader.get_prediction(self.dataset, self.slide, self.tile)

        for cell in prediction:
            pred_pos = Point(cell[0][1], cell[0][0])
            if self.polygon.contains(pred_pos):
                pred_pheno = cell[1]
                pred_cell_type = panel.inform_phenotype(0, pred_pheno)
                if pred_cell_type.main in panel.foreground_types:
                    cd45ro_pred = int("cd45ro+" in pred_pheno.lower())
                    lymphocyte = Lymphocyte(pred_pos.x, pred_pos.y, pred_cell_type, cd45ro_pred)
                    self.pred_cells_inform.append(lymphocyte)


def load_val_rois(rois_file_path):
    with gzip.open(rois_file_path) as f:
        rois_json = json.loads(f.read())

    rois = [ROI.from_dict(roi_dict) for roi_dict in rois_json]

    return rois
