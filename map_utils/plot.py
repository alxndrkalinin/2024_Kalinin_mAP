from functools import partial
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import font_manager
from matplotlib.ticker import FuncFormatter, MultipleLocator
from pkg_resources import resource_filename


def x_axis_formatter(x, pos, n_labels=5):
    if pos == 1 or pos == n_labels - 1:
        return f"{int(x)}"
    else:
        return f"{x:.2f}"


def remove_inner_ticklabels(fig: plt.Figure):
    for ax in fig.axes:
        try:
            ax.label_outer()
        except AttributeError:
            pass


def add_corner_text_annotations(
    ax, data, top="exclusive", bottom="dropped", prefix="", h_offset=0.05, v_offset=0.05
):
    threshold = -np.log10(0.05)
    dropped_percent = (data[bottom] > threshold).mean()
    exclusive_percent = (data[top] > threshold).mean()
    ax.text(
        h_offset,
        1 - v_offset,
        f"{prefix}{exclusive_percent:.0%}",
        transform=ax.transAxes,
        ha="left",
        va="top",
    )
    ax.text(
        1 - 2 * h_offset,
        v_offset,
        f"{prefix}{dropped_percent:.0%}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
    )


def set_plotting_style():
    font_path = resource_filename(__name__, "fonts/OpenSans-Regular.ttf")
    font_manager.fontManager.addfont(font_path)
    _ = font_manager.FontProperties(fname=font_path)
    plt.rcParams["font.family"] = ["Open Sans"]
    plt.rcParams["font.size"] = 16


def format_fig3a(ax_scatter, pr_x, pr_y, frs, palette, fr_total):
    """Apply Fig3A specific formatting to the mAP scatter plot."""
    ax_scatter.text(pr_x - 0.33, pr_y, "Retrieved: ", transform=ax_scatter.transAxes)
    ax_scatter.text(
        0.82,
        0.12,
        "p=0.05",
        transform=ax_scatter.transAxes,
        color="grey",
        fontsize=12,
        fontstyle="italic",
    )
    for i, (hue_value, fr) in enumerate(frs.items()):
        ax_scatter.text(
            pr_x + i * 0.15,
            pr_y,
            f"{fr:.0%}",
            transform=ax_scatter.transAxes,
            color=palette[hue_value],
        )
    ax_scatter.text(
        pr_x + len(frs) * 0.15,
        pr_y,
        f"({fr_total:.0%})",
        transform=ax_scatter.transAxes,
    )


def format_fig3b(ax_scatter, pr_x, pr_y, frs, palette, fr_total):
    """Apply Fig3B specific formatting to the mAP scatter plot."""
    ax_scatter.text(pr_x, pr_y + 0.07, "Retrieved:", transform=ax_scatter.transAxes)
    ax_scatter.text(
        0.01,
        pr_y - 0.03,
        "p=0.05",
        transform=ax_scatter.transAxes,
        color="grey",
        fontsize=12,
        fontstyle="italic",
    )
    for i, (hue_value, fr) in enumerate(frs.items()):
        ax_scatter.text(
            (pr_x - 0.2) + i * 0.16,
            pr_y - 0.03,
            f"{fr:.0%}",
            transform=ax_scatter.transAxes,
            color=palette[hue_value],
        )
    ax_scatter.text(
        (pr_x - 0.2) + len(frs) * 0.16,
        pr_y - 0.03,
        f"({fr_total:.0%})",
        transform=ax_scatter.transAxes,
    )


def format_fig3e(ax_scatter, pr_x, pr_y, frs, palette, fr_total):
    """Apply Fig3E specific formatting to the mAP scatter plot."""
    ax_scatter.text(pr_x, pr_y + 0.1, "Retrieved:", transform=ax_scatter.transAxes)
    ax_scatter.text(
        0.82,
        pr_y - 0.45,
        "p=0.05",
        transform=ax_scatter.transAxes,
        color="grey",
        fontsize=12,
        fontstyle="italic",
    )
    for i, (hue_value, fr) in enumerate(frs.items()):
        ax_scatter.text(
            pr_x + i * 0.16,
            pr_y,
            f"{fr:.0%}",
            transform=ax_scatter.transAxes,
            color=palette[hue_value],
        )
    ax_scatter.text(
        pr_x + len(frs) * 0.16,
        pr_y,
        f"({fr_total:.0%})",
        transform=ax_scatter.transAxes,
    )


def apply_figure_formatting(
    figure, ax_scatter, pr_x, pr_y, frs, palette, fr_total, row_value=None, y_label=None
):
    """Apply formatting specific to different figure types."""
    if figure == "Fig3A":
        format_fig3a(ax_scatter, pr_x, pr_y, frs, palette, fr_total)
        if y_label is not None:
            ax_scatter.set_ylabel(f"Preprocessing: {row_value}\n-log10(mAP p-value)")
    elif figure == "Fig3B":
        format_fig3b(ax_scatter, pr_x, pr_y, frs, palette, fr_total)
    elif figure == "Fig3E":
        format_fig3e(ax_scatter, pr_x, pr_y, frs, palette, fr_total)


def plot_map_scatter(
    df,
    col,
    title,
    metric="mAP",
    row=None,
    hue_col=None,
    palette=None,
    move_legend="lower right",
    y_label=None,
    aspect=None,
    adjust=None,
    pr_x=0.61,
    pr_y=0.35,
    l_x=1.05,
    l_y=0.575,
    m_x=0.52,
    m_y=0.01,
    kde_y=0.4,
    figure=None,
    save_path=None,
):
    """
    Unified function for plotting mAP data with optional hue grouping.

    Parameters:
    -----------
    df : pandas.DataFrame
        The data to plot
    col : str
        Column to use for subplot columns
    title : str
        Plot title
    metric : str, default="mAP"
        Metric to plot on x-axis
    row : str, optional
        Column to use for subplot rows
    hue_col : str, optional
        Column to use for hue grouping. If None, "p < 0.05" is used as hue
    palette : str or dict, optional
        Color palette to use for hue groups
    move_legend : str, default="lower right"
        Position for the legend
    y_label : str, optional
        Custom y-axis label
    aspect : float, optional
        Aspect ratio of the subplots
    adjust : dict, optional
        Parameters for fig.subplots_adjust()
    pr_x, pr_y : float
        Position for retrieved percentage annotation
    l_x, l_y : float
        Position for legend
    m_x, m_y : float
        Position for metric label
    kde_y : float
        Vertical position for KDE annotations
    figure : str, optional
        Figure identifier for specific formatting (e.g., "Fig3A")
    save_path : str, optional
        Path to save the figure

    Returns:
    --------
    fig : matplotlib.figure.Figure
        The figure object
    """
    unique_col_values = df[col].unique()
    num_cols = len(unique_col_values)
    num_rows = 1
    unique_row_values = [None]

    # Determine if we're using hue mode
    use_hue_mode = hue_col is not None

    # Set default hue column if in hue mode
    if use_hue_mode:
        hue_col = hue_col if hue_col is not None else "p < 0.05"
        palette = palette if palette is not None else "Set1"

    # Process row information if provided
    if row is not None and row in df.columns:
        unique_row_values = sorted(df[row].unique())[::-1]
        num_rows = len(unique_row_values)

    # Create figure and axes
    fig, axes = plt.subplots(
        num_rows,
        num_cols,
        figsize=(num_cols * 4, num_rows * 4),
        layout="constrained",
        sharex=True,
        sharey=True,
    )

    for row_i, row_value in enumerate(unique_row_values):
        row_df = df[df[row] == row_value] if row is not None else df
        row_axes = axes[row_i] if num_rows > 1 else axes

        for col_i, col_value in enumerate(unique_col_values):
            sub_df = row_df[row_df[col] == col_value]
            mean_map = sub_df[metric].mean()
            fr_total = sub_df["p < 0.05"].mean()

            # Set up the scatter plot axis and the KDE inset axis
            ax_scatter = row_axes[col_i] if num_cols > 1 else axes
            if aspect is not None:
                ax_scatter.set(aspect=aspect)
            ax_kde = ax_scatter.inset_axes([0, 1, 1, 0.15], sharex=ax_scatter)

            # Create scatter plot with appropriate styling based on mode
            if use_hue_mode:
                frs = sub_df.groupby(hue_col)["p < 0.05"].mean()
                sns.scatterplot(
                    ax=ax_scatter,
                    data=sub_df,
                    x=metric,
                    y=f"-log10({metric} p-value)",
                    style="markers",
                    hue=hue_col,
                    palette=palette,
                    markers={"p<0.05": "o", "p>=0.05": "s"},
                    s=50,
                )
                # Add horizontal line for p=0.05
                ax_scatter.axhline(-np.log10(0.05), color="grey", linestyle="--")

                # Apply figure-specific formatting if applicable
                if figure in ["Fig3A", "Fig3B", "Fig3E"]:
                    apply_figure_formatting(
                        figure,
                        ax_scatter,
                        pr_x,
                        pr_y,
                        frs,
                        palette,
                        fr_total,
                        row_value,
                        y_label,
                    )
            else:
                sns.scatterplot(
                    ax=ax_scatter,
                    data=sub_df,
                    x=metric,
                    y=f"-log10({metric} p-value)",
                    hue="p < 0.05",
                    s=50,
                )
                ax_scatter.set_title(f"{col_value}", fontsize=16, pad=20)

                # Simple retrieved percentage for non-hue mode
                if col_i == 1:
                    ax_scatter.text(
                        pr_x,
                        0.02,
                        f"Retrieved: {fr_total:.0%}",
                        transform=ax_scatter.transAxes,
                    )
                else:
                    ax_scatter.text(
                        pr_x,
                        pr_y,
                        f"Retrieved: {fr_total:.0%}",
                        transform=ax_scatter.transAxes,
                    )

            # Common scatter plot formatting
            ax_scatter.set_xlabel("")
            ax_scatter.set_xlim(-0.05, 1.05)
            ax_scatter.set_ylim(-0.1, max(sub_df[f"-log10({metric} p-value)"]) + 0.25)
            ax_scatter.xaxis.set_major_locator(MultipleLocator(base=0.25))
            ax_scatter.get_xaxis().set_major_formatter(
                FuncFormatter(partial(x_axis_formatter, n_labels=6))
            )

            # Remove legend from subplot and store handles and labels for main legend
            handles, labels = ax_scatter.get_legend_handles_labels()
            if use_hue_mode:
                labels = [label.replace("markers", "") for label in labels]
            ax_scatter.get_legend().remove()

            # Create KDE plots at the top of each scatter plot
            mmaps = {}
            max_kde_y = 0

            if use_hue_mode:
                for hue_val in sorted(sub_df[hue_col].unique()):
                    sns.kdeplot(
                        ax=ax_kde,
                        data=sub_df[(sub_df[hue_col] == hue_val)],
                        x=metric,
                        label=str(hue_val),
                        cut=1,
                        color=palette[hue_val],
                    )
                    mmaps[hue_val] = sub_df[(sub_df[hue_col] == hue_val)][metric].mean()
                    max_kde_y = max(max_kde_y, max(ax_kde.lines[-1].get_ydata()))
            else:
                for p_value in sorted(sub_df["p < 0.05"].unique()):
                    sns.kdeplot(
                        ax=ax_kde,
                        data=sub_df[sub_df["p < 0.05"] == p_value],
                        x=metric,
                        label=str(p_value),
                    )
                    max_kde_y = max(max_kde_y, max(ax_kde.lines[-1].get_ydata()))

                # Add mean line and label for non-hue mode
                ax_kde.axvline(mean_map, color="grey", linestyle="--")
                if mean_map < 0.5:
                    ax_kde.text(
                        mean_map + 0.05,
                        kde_y,
                        f"Mean {metric}: {mean_map:.2f}",
                        transform=ax_kde.transAxes,
                    )
                else:
                    ax_kde.text(
                        mean_map - 0.2,
                        kde_y,
                        f"Mean {metric}: {mean_map:.2f}",
                        transform=ax_kde.transAxes,
                    )

            # Format the KDE subplot
            ax_kde.get_yaxis().set_visible(False)
            ax_kde.get_xaxis().set_visible(False)
            sns.despine(ax=ax_kde, left=True, bottom=True)

    # Add metric label to the figure
    fig.text(m_x, m_y, metric, ha="center", va="center")

    # Handle layout and legend
    plt.tight_layout()
    fig.legend(
        handles,
        labels,
        title="" if use_hue_mode else "p < 0.05",
        loc="upper center",
        bbox_to_anchor=(l_x, l_y),
        frameon=False,
    )

    if adjust is not None:
        fig.subplots_adjust(**adjust)

    # Save figure if path provided
    if save_path is not None and figure is not None:
        fig.savefig(Path(save_path) / f"{figure}.png", dpi=300)
        fig.savefig(Path(save_path) / f"{figure}.svg")

    plt.show()
    return fig
