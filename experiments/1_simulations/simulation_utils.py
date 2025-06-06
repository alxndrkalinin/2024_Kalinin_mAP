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
    y_label="Recall",
    title_x="Number of controls",
    title_y="Number of features",
    height=3,
    ylim=(0, 1),
    font_size=16,
    markersize=10,
    labelrotation=25,
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
    # plt.rcParams.update({"font.size": font_size, "axes.labelsize": font_size, "axes.titlesize": font_size})
    palette = sns.color_palette("tab10") if palette is None else palette

    df["features_differ"] = df["features_differ"].astype(str)
    df_melted = df.melt(
        id_vars=["features_differ", "# replicates", "n_controls", "n_feats"],
        value_vars=metrics,
        var_name="metric",
        value_name="value",
    )

    g = sns.FacetGrid(
        df_melted, row="n_feats", col="n_controls", height=height, aspect=1.0, ylim=ylim
    )

    plt.subplots_adjust(bottom=0.4)
    g.figure.text(
        0.45,
        0.00,
        "% features that differ from control",
        ha="center",
        fontsize=font_size,
    )

    g.map_dataframe(
        sns.lineplot,
        x="features_differ",
        y="value",
        hue="metric",
        style="# replicates",
        markers=True,
        dashes=True,
        markersize=markersize,
        palette=palette,
    )

    g.add_legend(bbox_to_anchor=(1.1, 0.5), loc="center right", fontsize=font_size)
    g.set_axis_labels("", y_label)
    g.set_titles("# controls={col_name}")

    for ax in g.axes.flatten():
        ax.title.set_fontsize(font_size)
        ax.set_ylabel(ax.get_ylabel(), fontsize=font_size)
        ax.tick_params(labelsize=font_size, labelrotation=labelrotation)

    for ax, title in zip(g.axes[:, -1], g.row_names):
        ax.text(
            1.1,
            0.55,
            f"# features={title}",
            rotation=270,
            verticalalignment="center",
            horizontalalignment="right",
            transform=ax.transAxes,
            fontsize=font_size,
        )

    fig = g.figure
    plt.text(0.36, 1.0, title_x, fontsize=font_size, transform=fig.transFigure)
    plt.text(
        0.83, 0.45, title_y, fontsize=font_size, rotation=-90, transform=fig.transFigure
    )

    fig.set_size_inches(len(g.col_names) * height * 1.2, len(g.row_names) * height)
    plt.show()

    return fig
