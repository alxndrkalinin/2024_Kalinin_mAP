import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def get_umap_palette(n=3):
    tab10 = sns.color_palette("tab10", n_colors=2)
    pastel_grey = sns.color_palette("pastel")[-3]
    umap_palette = [pastel_grey] + tab10
    return umap_palette


def color_category(row):
    if row["Gene p < 0.05"] and row["MC p < 0.05"]:
        return "Both"
    elif row["Gene p < 0.05"]:
        return "Gene only"
    elif row["MC p < 0.05"]:
        return "MC only"
    else:
        return "None"


def plot_gene_vs_ms(
    aps,
    title=None,
    width=3,
    height=3,
    point_size=10,
    point_alpha=0.75,
    title_font_size=16,
    legend_markersize=5,
):
    aps["p < 0.05"] = aps.apply(color_category, axis=1)
    aps = aps.sort_values(
        by="p < 0.05",
        key=lambda x: x.map({"None": 4, "Gene only": 3, "MC only": 2, "Both": 1}),
    )

    palette = list(sns.color_palette())

    color_map = {
        "None": palette[0],  # Both False
        "Gene only": palette[4],  # Gene p < 0.05 is True, MC p < 0.05 is False
        "MC only": palette[5],  # MC p < 0.05 is True, Gene p < 0.05 is False
        "Both": palette[1],  # Both are True
    }

    g = sns.jointplot(
        data=aps,
        x="AP, Gene",
        y="AP, Morphological Class (MC)",
        hue="p < 0.05",
        palette=color_map,
        s=point_size,
        alpha=point_alpha,
        marginal_kws={"cut": 0, "fill": False},
    )

    g.ax_joint.plot([0, 1], [0, 1], color="grey", linestyle="--")

    handles, labels = g.ax_joint.get_legend_handles_labels()
    label_order = np.argsort(labels)
    labels = np.asarray(labels)[label_order].tolist()
    handles = np.asarray(handles)[label_order].tolist()

    for handle in handles:
        handle.set_markersize(legend_markersize)

    g.ax_joint.get_legend().remove()
    g.ax_joint.legend(
        handles=handles,
        labels=labels,
        title="p < 0.05",
        loc="center left",
        bbox_to_anchor=(1.2, 0.5),
        ncol=1,
        frameon=False,
    )

    g.ax_joint.text(
        0.05,
        1.0,
        f"mAP: {aps['AP, Morphological Class (MC)'].mean():.2f}\nRetrieved: {aps['MC p < 0.05'].mean():.0%}",
        transform=g.ax_joint.transAxes,
        ha="left",
        va="top",
    )
    g.ax_joint.text(
        0.98,
        0.005,
        f"mAP: {aps['AP, Gene'].mean():.2f}\nRetrieved: {aps['Gene p < 0.05'].mean():.0%}",
        transform=g.ax_joint.transAxes,
        ha="right",
        va="bottom",
    )

    if title is not None:
        plt.suptitle(title, fontsize=title_font_size)
    plt.show()

    fig = plt.gcf()
    fig.set_size_inches(width, height)
    return fig


def plot_mc_cp_vs_dp(
    aps,
    hue=None,
    legend=True,
    width=3,
    height=3,
    point_size=20,
    point_alpha=0.75,
    legend_markersize=5,
):
    mean_aps = aps[["AP, CellProfiler features", "AP, DeepProfiler features"]].mean()
    retrieved = aps[["CellProfiler retrieved", "DeepProfiler retrieved"]].mean()

    fig, ax = plt.subplots(figsize=(width, height))
    g = sns.scatterplot(
        aps,
        x="AP, CellProfiler features",
        y="AP, DeepProfiler features",
        hue=hue,
        s=point_size,
        alpha=point_alpha,
        ax=ax,
    )

    if legend:
        g.legend(
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            ncol=1,
            title="Morphological Class",
            frameon=False,
        )
    else:
        plt.legend([], [], frameon=False)

    ax.plot([0, 1], [0, 1], color="grey", linestyle="--")

    handles, labels = g.get_legend_handles_labels()
    label_order = np.argsort(labels)
    labels = np.asarray(labels)[label_order].tolist()
    handles = np.asarray(handles)[label_order].tolist()

    for handle in handles:
        handle.set_markersize(legend_markersize)

    ax.text(
        0.01,
        0.99,
        f"mAP: {mean_aps['AP, DeepProfiler features']:.2f}, Retrieved: {retrieved['DeepProfiler retrieved']:.0%}",
        transform=ax.transAxes,
        ha="left",
        va="top",
    )
    ax.text(
        0.99,
        0.00,
        f"mAP: {mean_aps['AP, CellProfiler features']:.2f}, Retrieved: {retrieved['CellProfiler retrieved']:.0%}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
    )
    return fig
