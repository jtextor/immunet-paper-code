import pandas as pd
from panels import load_panels
from utils import data_path

panel = load_panels()["lymphocyte"]

DETECTED_COLUMN = "detected"
ERROR_RATE_COLUMN = "err_rate"
CASES_COLUMN = "cases"
ERRORS_COLUMN = "errors"


def evaluate_immunet(
    prediction_file_path,
    data,
    output_folder,
    subtypes=True
):
    df = pd.read_csv(prediction_file_path, sep="\t", index_col=False)
    df = df.assign(foreground=lambda x: x.ann_type.isin(["B cell", "T cell"]))

    if subtypes:
        df["ann_subtype"] = df.apply(lambda x: panel.annotation_phenotype(x["ann_type"],
                                                                          [x["CD3"], x["FOXP3"], x["CD20"], x["CD45RO"],
                                                                           x["CD8"]]).detailed, axis=1)
        df = df[~df.ann_subtype.isin(["Invalid", "Memory T cell"])]

    def detected(row):
        distance = row["dist"]

        if not row["foreground"]:
            return int(distance > 7)

        annotated_type = row["ann_type"]

        predicted_pos = [row["pCD3"], row["pFOXP3"], row["pCD20"], row["pCD45RO"], row["pCD8"]]
        predicted_pheno = panel.prediction_phenotype(distance, predicted_pos)
        if subtypes:
            annotated_pos = [row["CD3"], row["FOXP3"], row["CD20"], row["CD45RO"], row["CD8"]]
            annotated_pheno = panel.annotation_phenotype(annotated_type, annotated_pos)

            return int((distance <= 7) & (predicted_pheno.detailed == annotated_pheno.detailed))
        else:
            return int((distance <= 7) & (predicted_pheno.main == annotated_type))


    df[DETECTED_COLUMN] = df.apply(lambda x: detected(x), axis=1)

    if subtypes:
        df["ann_subtype"] = df.apply(lambda x: x["ann_subtype"] if x["foreground"] else "Other cell", axis=1)
    else:
        df["ann_type"] = df.apply(lambda x: x["ann_type"] if x["foreground"] else "Other cell", axis=1)

    groupby_type_col = "ann_subtype" if subtypes else "ann_type"
    group = df.groupby(["tissue", groupby_type_col])

    err_rates = 1 - group[DETECTED_COLUMN].mean()
    nums = group[DETECTED_COLUMN].count()
    err_cases = nums - group[DETECTED_COLUMN].sum()

    errors_df = pd.concat([err_rates.to_frame(), nums.to_frame(), err_cases.to_frame()], axis=1)
    errors_df.columns = [ERROR_RATE_COLUMN, CASES_COLUMN, ERRORS_COLUMN]
    errors_df = errors_df.reset_index()

    suff = "subtypes" if subtypes else "TB"
    errors_df.to_csv(output_folder / "immunet_perf_{}_{}.csv".format(data, suff))


def evaluate_inform(
    prediction_file_path,
    desc,
    output_folder,
    subtypes=True
):
    df = pd.read_csv(prediction_file_path, sep="\t", index_col=False)

    df = df.assign(foreground=lambda x: x.ann_type.isin(["B cell", "T cell"]))
    df["pred_type"] = df.apply(lambda x: x["pred_type"] if isinstance(x["pred_type"], str) else "Other cell", axis=1)

    if subtypes:
        df["ann_subtype"] = df.apply(lambda x: panel.annotation_phenotype(x["ann_type"],
                                                                          [x["CD3"], x["FOXP3"], x["CD20"], x["CD45RO"],
                                                                           x["CD8"]]).detailed, axis=1)
        df = df[~df.ann_subtype.isin(["Invalid", "Memory T cell"])]

    df["ann_type"] = df.apply(lambda x: x["ann_type"] if x["foreground"] else "Other cell", axis=1)

    if subtypes:
        df["ann_subtype"] = df.apply(lambda x: x["ann_subtype"] if x["foreground"] else "Other cell", axis=1)

    def detected(row):
        distance = row["dist"]
        annotated_type = row["ann_type"]
        predicted_pheno = row["pred_type"]
        predicted_type = panel.inform_phenotype(distance, predicted_pheno)

        if not row["foreground"]:
            return int((distance > 7) | (predicted_type.main == annotated_type))

        if subtypes:
            annotated_pos = [row["CD3"], row["FOXP3"], row["CD20"], row["CD45RO"], row["CD8"]]
            annotated_pheno = panel.annotation_phenotype(annotated_type, annotated_pos)

            return int((predicted_type.detailed == annotated_pheno.detailed) & (row["dist"] <= 7))
        else:
            return int((predicted_type.main == annotated_type) & (row["dist"] <= 7))

    df[DETECTED_COLUMN] = df.apply(lambda x: detected(x), axis=1)

    groupby_type_col = "ann_subtype" if subtypes else "ann_type"
    group = df.groupby(["tissue", groupby_type_col])

    err_rates = 1 - group[DETECTED_COLUMN].mean()
    nums = group[DETECTED_COLUMN].count()
    err_cases = nums - group[DETECTED_COLUMN].sum()

    errors_df = pd.concat([err_rates.to_frame(), nums.to_frame(), err_cases.to_frame()], axis=1)
    errors_df.columns = [ERROR_RATE_COLUMN, CASES_COLUMN, ERRORS_COLUMN]
    errors_df = errors_df.reset_index()

    suff = "subtypes" if subtypes else "TB"
    errors_df.to_csv(output_folder / "inform_perf_{}_{}.csv".format(desc, suff))


if __name__ == '__main__':
    # Evaluate ImmuNet on validation data
    immunet_prediction_val_path = data_path / "prediction-immunet-val.tsv"
    evaluate_immunet(immunet_prediction_val_path, "val", data_path, True)

    # Evaluate ImmuNet on training data
    immunet_prediction_train_path = data_path / "prediction-immunet-train.tsv"
    evaluate_immunet(immunet_prediction_train_path, "train", data_path, True)

    # Evaluate new InForm analysis on validation data
    inform_new_prediction_val_path = data_path / "prediction-inform-val.tsv"
    evaluate_inform(inform_new_prediction_val_path, "val", data_path, True)

    # Evaluate new InForm analysis on training data
    inform_new_prediction_train_path = data_path / "prediction-inform-train.tsv"
    evaluate_inform(inform_new_prediction_train_path, "train", data_path, True)
