import networkx as nx
import math
import numpy as np
import pandas as pd
from rois import load_val_rois, panel
from utils import DS_TS, ImmuNetPredictionLoader, InFormPredictionLoader, immunet_pred_path, inform_pred_path, \
    rois_json_path, data_path, model_path, images_path
# ImmuNet: uncomment to use ImmuNet to find cells on tiles
# from inference import ImmuNetPredictionHandler

distance_cutoff_in_pixels = 7.0
positivity_cutoff = 0.4

def match_dist(a, b):
    def dist(v, w):
        return math.sqrt((v[0]-w[0])**2 + (v[1]-w[1])**2)
    a = [[i[0],i[1]] for i in a]

    if len(a) == 0 or len(b) == 0:
        return None

    G = nx.Graph()
    leftside = ["c"+str(i) for i in range(len(a))]
    rightside  = ["s"+str(i) for i in range(len(b))]
    G.add_nodes_from(leftside, bipartite=0)
    G.add_nodes_from(rightside, bipartite=1)
    for (i,c) in enumerate(a):
        for (j,s) in enumerate(b):
            if dist(s,c) <= distance_cutoff_in_pixels:
                G.add_edge( "c"+str(i), "s"+str(j))
    m = nx.bipartite.maximum_matching(G, leftside)
    return m


def analyse_box(annotated_cells, detected_cells):
    detection_results = []
    cd45ro_results = []
    cell_types = set(panel.foreground_types).union(panel.foreground_subtypes)

    annotations_coords = np.array([[cell.x, cell.y] for cell in annotated_cells])
    detections_coords = np.array([[cell.x, cell.y] for cell in detected_cells])

    matching = match_dist(annotations_coords, detections_coords)
    for cell_type in cell_types:
        # Get indices of annotations and detections with a specific cell type
        annotations_inds = [i for i, cell in enumerate(annotated_cells) if cell_type in cell.types]
        detections_inds = [i for i, cell in enumerate(detected_cells) if cell_type in cell.types]

        # Get the corresponding name of vertices in graph
        leftside = ["c" + str(i) for i in annotations_inds]
        rightside = ["s" + str(i) for i in detections_inds]

        if matching is not None and len(leftside) > 0 and len(rightside) > 0:
            # All non-matched annotations are false negatives
            false_negatives = sum([int(l not in matching) for l in leftside])
            # All non-matched annotations are false positives
            false_positives = sum([int(l not in matching) for l in rightside])

            # For all matched annotations checked if they are linked to
            # a prediction with the same cell type
            annot_sub_graph = [(key, value) for (key, value) in matching.items() if
                          key[0] == "c" and int(key[1:]) in annotations_inds]
            for key, value in annot_sub_graph:
                # Count those which have different types as false negatives
                if value not in rightside:
                    false_negatives += 1

            # For all matched detections checked if they are linked to
            # an annotation with the same cell type
            detections_sub_graph = [(key, value) for (key, value) in matching.items() if
                               key[0] == "s" and int(key[1:]) in detections_inds]
            for key, value in detections_sub_graph:
                # Count those which have different types as false positives
                if value not in leftside:
                    false_positives += 1
        else:
            if len(leftside) == 0:
                false_negatives, false_positives = 0, len(rightside)

            if len(rightside) == 0:
                false_negatives, false_positives = len(leftside), 0

        true_positives = len(annotations_inds) - false_negatives

        detection_results.append((true_positives, false_negatives, false_positives, cell_type))

        if matching is not None and cell_type in panel.decorable_types:
            matched_ct = [(key, value) for (key, value) in matching.items() if key[0] == "c" and int(key[1:]) in annotations_inds]
            annot = [int(key[1:]) for key, _ in matched_ct]
            detect = [int(key[1:]) for _, key in matched_ct if key[0] == 's']

            annot_cd45ro = [(annotated_cells[i]).cd45ro for i in annot]
            detect_cd45ro = [(detected_cells[i]).cd45ro for i in detect]

            for ann_cd45ro, det_cd45ro in zip(annot_cd45ro, detect_cd45ro):
                cd45ro_results.append([ann_cd45ro, det_cd45ro])

    return detection_results, cd45ro_results


def analyse():
    f_score_results = []
    cd45ro_corr_results = []
    roi_data = []

    immunet_prediction_loader = ImmuNetPredictionLoader(immunet_pred_path)
    inform_prediction_loader = InFormPredictionLoader(inform_pred_path)

    # ImmuNet: uncomment to use ImmuNet to find cells on tiles
    #prediction_handler = ImmuNetPredictionHandler(model_path, images_path)

    rois = load_val_rois(rois_json_path)
    rois = [roi for roi in rois if roi.dataset != "2020-01-27-phenotyping-paper-cytoagars"]

    for roi in rois:
        print("Dataset {}, slide {}, tile{}".format(roi.dataset, roi.slide, roi.tile))

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

        annotated_cells = roi.ann_cells
        roi_info = {"roi_id": roi.id, "annotations": len(annotated_cells), 'dataset': DS_TS[roi.dataset], 'area': roi.area}
        roi_data.append(roi_info)

        for algo in ["NN", "inform"]:

            detected_cells = roi.pred_cells_immunet if algo == "NN" else roi.pred_cells_inform
            detection_results, cd45ro_results = analyse_box(annotated_cells, detected_cells)

            for result in detection_results:
                # ROIs csv separately
                tp, fn, fp, cell_type = result
                new_row = {"roi_id": roi.id, "tp": tp, "fn": fn, "fp": fp, "celltype": cell_type, "algorithm": algo}
                f_score_results.append(new_row)

            for result in cd45ro_results:
                # ROIs csv separately
                ann_cd45ro, det_cd45ro = result
                new_row = {"roi_id": roi.id, "ann_cd45ro": ann_cd45ro, "det_cd45ro": det_cd45ro, "algorithm": algo}
                cd45ro_corr_results.append(new_row)

    pd.DataFrame(roi_data).to_csv(data_path / "rois_info.csv", index=False)
    pd.DataFrame(f_score_results).to_csv(data_path / "detection_stat.csv", index=False)
    pd.DataFrame(cd45ro_corr_results).to_csv(data_path / "cd45ro_stat.csv", index=False)

analyse()
