from pathlib import Path

import pandas as pd

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake import script


try:
    import seaborn as sns
    import matplotlib.pyplot as plt
except ImportError:
    pass


def plot_qc(df_plt):
    # context keywords: notebook, talk, paper, poster
    sns.set_theme(
        context="paper",
        style="whitegrid",
        rc={
            "axes.grid.axis": "y",
            "figure.figsize": [16 * 0.4, 9 * 0.4],
        },
    )
    metrics = ["Sensitivity", "Precision", "F1-score"]
    f, axes = plt.subplots(ncols=3, sharex="all", sharey="all")
    for metric, ax in zip(metrics, axes):
        sns.lineplot(
            data=df_plt,
            x="taxon_abbreviation",
            y=metric,
            hue="method",
            ax=ax,
            # legend only in first plot
            legend=bool(metric == metrics[0]),
            marker="s",
        )
        ax.set_title(metric)
        ax.set_ylim(-5, 105)
        ax.set_xlabel("Taxonomic rank")
        # x-axis label only in second plot
        ax.xaxis.label.set_visible(metric == metrics[1])
    axes[0].yaxis.label.set_visible(False)
    axes[-1].tick_params(labelright=True)
    sns.despine(left=True, bottom=True)
    f.tight_layout()
    return f


def qc_plots(
    qc_data_file: str,
    output: str,
):
    df_qc_data = pd.read_csv(
        qc_data_file,
        sep="\t",
    )
    f = plot_qc(df_qc_data)
    f.savefig(output)
    return f


def test_qc_plots(tmpdir):
    import test
    test_data_directory = Path(test.__file__).parent / "rules" / "qc_taxonomic_profiling" / "plot_qc_taxonomic_profiling" / "data"
    qc_data_file = test_data_directory / "results" / "qc" / "qc-foo.tsv"
    output_plot = tmpdir / "qc_plot.svg"
    print(output_plot)
    f = qc_plots(
        qc_data_file=qc_data_file.as_posix(),
        output=output_plot,
    )
    f.show()


snakemake: "script.Snakemake"
if snakemake := globals().get("snakemake"):
    with open(snakemake.log[0], "w") as log_handle:
        LOG_HANDLE = log_handle
        qc_plots(
            qc_data_file=snakemake.input.qc_data,
            output=snakemake.output[0],
        )
