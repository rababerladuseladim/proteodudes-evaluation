from pathlib import Path

import numpy as np
import pandas as pd

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake import script


TAX_LEVELS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "subspecies"]
TAXON_ABBREVIATIONS = {
    "superkingdom": "SK",
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
}

def read_ground_truth_file(ground_truth_file):
    df_t_taxids = pd.read_csv(
        ground_truth_file,
        usecols=TAX_LEVELS,
        sep="\t",
        dtype={
            c: pd.Int64Dtype()  # support nan values
            for c in TAX_LEVELS
        },
    )
    return df_t_taxids


def get_value_overlap(hits: pd.Series, ground_truth: pd.Series):
    """
    build series of values from both input series and annotate them with TP, FP, FN

    :param hits:
    :param ground_truth:
    :return: series with index: union of values from hits and ground_truth, values: TP if index in both,
        FP if only in hits, FN if only in ground_truth
    """
    left = {t for t in hits if not (pd.isna(t) or (t == "-"))}
    right = {t for t in ground_truth if not (pd.isna(t) or (t == "-"))}
    taxa = left.union(right)
    dct = dict()
    for t in taxa:
        found = t in left
        expected = t in right
        if found and expected:
            dct[t] = "TP"
        elif found:
            dct[t] = "FP"
        elif expected:
            dct[t] = "FN"
        else:
            raise KeyError(f"Taxon '{t}' missing in both Series!")
    return pd.Series(dct)


def get_diamond_hit_counts(f: Path | str, ground_truth_df_tax_ids: pd.DataFrame):
    df_dmnd = read_ground_truth_file(f)
    df = pd.DataFrame(columns=TAX_LEVELS, index=["TP", "FP", "FN"],)
    for t in TAX_LEVELS:
        s = get_value_overlap(df_dmnd[t], ground_truth_df_tax_ids[t])
        df[t] = s.value_counts()
    df = df.fillna(0)
    df["eval"] = df.index
    return df


def get_unipept_hit_counts(f, ground_truth_df_tax_ids):
    columns = [f"{t}_id" for t in TAX_LEVELS]
    columns = ["domain_id" if item == "superkingdom_id" else item for item in columns]
    df_uni = pd.read_csv(
        f,
        usecols=columns,
        dtype=pd.Int64Dtype()
    )
    df_uni.rename(columns={"domain_id": "superkingdom_id"}, inplace=True)
    df_uni.rename(lambda x: x.strip("_id"), axis="columns", inplace=True)
    df = pd.DataFrame(columns=TAX_LEVELS, index=["TP", "FP", "FN"],)
    for t in TAX_LEVELS:
        s = get_value_overlap(df_uni[t], ground_truth_df_tax_ids[t])
        df[t] = s.value_counts()
    df = df.fillna(0)
    df["eval"] = df.index
    return df


def get_dudes_hit_counts(f, ground_truth_df_tax_ids):
    df_dudes = pd.read_csv(f, sep="\t", skiprows=5)
    df_dudes = df_dudes.rename(lambda c: c.lower().strip("@"), axis=1)
    df = pd.DataFrame(columns=TAX_LEVELS, index=["TP", "FP", "FN"])
    for t in TAX_LEVELS:
        s_dudes = (
            df_dudes.loc[df_dudes["rank"] == t, "taxpath"]
            .map(lambda x: x.split("|")[-1])
            .astype(int)
        )
        s = get_value_overlap(s_dudes, ground_truth_df_tax_ids[t])
        s_classification = s.value_counts()
        for k in ["TP", "FP", "FN"]:
            if k not in s_classification:
                s_classification[k] = 0
        df[t] = s_classification
    df = df.fillna(0)
    df["eval"] = df.index
    return df


def calculate_sensitivity(df: pd.DataFrame) -> float:
    true_positives = df[df["eval"] == "TP"]["value"].values[0]
    false_negatives = df[df["eval"] == "FN"]["value"].values[0]
    if true_positives + false_negatives == 0:
        return np.nan
    return true_positives * 100 / (true_positives + false_negatives)


def calculate_precision(df: pd.DataFrame) -> float:
    true_positives = df[df["eval"] == "TP"]["value"].values[0]
    false_positives = df[df["eval"] == "FP"]["value"].values[0]
    if true_positives + false_positives == 0:
        return np.nan
    return true_positives * 100 / (true_positives + false_positives)


def calculate_fdr(df: pd.DataFrame) -> float:
    true_positives = df[df["eval"] == "TP"]["value"].values[0]
    false_positives = df[df["eval"] == "FP"]["value"].values[0]
    if true_positives + false_positives == 0:
        return np.nan
    return false_positives * 100 / (true_positives + false_positives)


def calc_eval_metrics(df_hits: pd.DataFrame) -> pd.DataFrame:
    df_hits_long = df_hits.melt(id_vars=["method", "eval"], var_name="taxon_level")
    df_hits_long["value"] = df_hits_long["value"].apply(int)
    grouped = df_hits_long.groupby(["method", "taxon_level"], sort=False)

    df_sens_g = grouped.apply(calculate_sensitivity)
    df_precision_g = grouped.apply(calculate_precision)
    df_f_one_g = (2 * df_precision_g * df_sens_g) / (df_precision_g + df_sens_g)
    df_fdr_g = grouped.apply(calculate_fdr)
    df_eval = pd.DataFrame(
        {
            "Sensitivity": df_sens_g,
            "Precision": df_precision_g,
            "F1-score": df_f_one_g,
            "FDR": df_fdr_g,
        }
    )
    df_eval = df_eval.reset_index()
    return df_eval


def calculate_metrics(
    ground_truth_file: str,
    unipept_file: str,
    dudes_file_list: list[str],
    output: str,
    diamond_file: str | None = None
):
    hits: list[pd.DataFrame] = []
    # read ground truth
    df_gt_taxids = read_ground_truth_file(ground_truth_file)
    # get diamond hits
    if diamond_file:
        df_dmnd = get_diamond_hit_counts(diamond_file, df_gt_taxids)
        df_dmnd["method"] = "DIAMOND"
        hits.append(df_dmnd)
    # get unipept hits
    df_uni = get_unipept_hit_counts(unipept_file, df_gt_taxids)
    df_uni["method"] = "Unipept"
    hits.append(df_uni)
    # get proteodudes hits
    for dudes_file in dudes_file_list:
        method_name = "ProteoDUDes"
        if len(dudes_file_list) > 1:
            method_name += " on " + Path(dudes_file).parent.name
        df_dudes = get_dudes_hit_counts(
                dudes_file, df_gt_taxids
            )
        df_dudes["method"] = method_name
        hits.append(df_dudes)
    # build TP/FP/FN dataframe
    df_hits = pd.concat(hits, ignore_index=True)
    # build evaluation dataframe, containing sensitivity, precision, f1-score, fdr
    df_eval = calc_eval_metrics(df_hits)
    # drop subspecies, as no method has any true hits on this level
    df_plt = df_eval[df_eval["taxon_level"] != "subspecies"]
    df_plt["taxon_abbreviation"] = df_plt["taxon_level"].map(TAXON_ABBREVIATIONS)
    df_plt.to_csv(output, sep="\t")
    return df_plt


def test_calculate_metrics(tmpdir):
    import test
    test_data_directory = Path(test.__file__).parent / "scripts" / "plot_qc_taxonomic_profiling" / "test_data"
    ground_truth = test_data_directory / "ground_truth.csv"
    unipept_file = test_data_directory / "unipept.tsv"
    diamond_file = test_data_directory / "diamond.tsv"
    dudes_file = test_data_directory / "proteodudes.tsv"
    output_tsv = tmpdir / "qc_metrics.tsv"
    print(output_tsv)
    calculate_metrics(
        ground_truth_file=ground_truth.as_posix(),
        diamond_file=diamond_file.as_posix(),
        unipept_file=unipept_file.as_posix(),
        dudes_file_list=[dudes_file.as_posix()],
        output=output_tsv,
    )


snakemake: "script.Snakemake"
if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        calculate_metrics(
            ground_truth_file=snakemake.input.ground_truth,
            diamond_file=snakemake.input.get("diamond_result"),
            unipept_file=snakemake.input.unipept_result,
            dudes_file_list=snakemake.input.proteodudes_results,
            output=snakemake.output[0],
        )
