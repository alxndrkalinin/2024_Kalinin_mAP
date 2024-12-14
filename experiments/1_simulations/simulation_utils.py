import seaborn as sns
from matplotlib import pyplot as plt


def get_figure2_palette():
    tab10 = sns.color_palette("tab10")
    return [tab10[i] for i in [1, 0, 2, 4]]


def get_figureS5_palette():
    tab20c = sns.color_palette("tab20c")
    tab10 = sns.color_palette("tab10", n_colors=3)
    extra_oranges = [tab20c[i] for i in [4, 6]]
    return sns.color_palette([tab10[1]] + extra_oranges)


def plot_simulation_results(
    df,
    metrics,
    palette=None,
    save_path=None,
    y_label="Recall",
    title_x="Number of controls",
    title_y="Number of features",
    ylim=(0, 1),
):
    """
    Plotting function for simulation results.

    Parameters:
    - df: DataFrame containing the data.
    - metrics: List of metrics to plot.
    - palette: Color palette for the plot.
    - save_path: Path to save the plot (optional).
    - y_label: Label for the y-axis.
    - title_x: Title for the x-axis.
    - title_y: Title for the y-axis.
    - ylim: Tuple defining the y-axis limits.
    """
    plt.rcParams.update({"font.size": 14, "axes.labelsize": 14, "axes.titlesize": 14})
    palette = sns.color_palette("tab10") if palette is None else palette

    df["features_differ"] = df["features_differ"].astype(str)
    df_melted = df.melt(
        id_vars=["features_differ", "# replicates", "n_controls", "n_feats"],
        value_vars=metrics,
        var_name="metric",
        value_name="value",
    )

    g = sns.FacetGrid(
        df_melted, row="n_feats", col="n_controls", height=3, aspect=1.0, ylim=ylim
    )

    plt.subplots_adjust(bottom=0.33)
    g.fig.text(
        0.45, 0.01, "% features that differ from control", ha="center", fontsize="large"
    )

    g.map_dataframe(
        sns.lineplot,
        x="features_differ",
        y="value",
        hue="metric",
        style="# replicates",
        markers=True,
        dashes=True,
        markersize=10,
        palette=palette,
    )

    g.add_legend(bbox_to_anchor=(1.1, 0.5), loc="center right", fontsize="large")
    g.set_axis_labels("", y_label)
    g.set_titles("# controls={col_name}", fontsize=14)

    for ax in g.axes.flatten():
        ax.tick_params(labelsize=12, labelrotation=25)

    for ax, title in zip(g.axes[:, -1], g.row_names):
        ax.text(
            1.1,
            0.55,
            f"# features={title}",
            rotation=270,
            verticalalignment="center",
            horizontalalignment="right",
            transform=ax.transAxes,
            fontsize=14,
        )

    plt.text(0.36, 1.0, title_x, fontsize=16, transform=plt.gcf().transFigure)
    plt.text(
        0.83, 0.45, title_y, fontsize=16, rotation=-90, transform=plt.gcf().transFigure
    )

    if save_path is not None:
        plt.savefig(f"figures/{save_path}.png", bbox_inches="tight")
        plt.savefig(f"figures/{save_path}.svg")

    plt.show()
