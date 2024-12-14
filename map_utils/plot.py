from pathlib import Path
from functools import partial
from pkg_resources import resource_filename

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator, FuncFormatter


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


def plot_map_x3(
    df,
    col,
    title,
    metric="mAP",
    row=None,
    move_legend="lower right",
    aspect=None,
    adjust=None,
    pr_x=0.61,
    pr_y=0.35,
    l_x=1.05,
    l_y=0.575,
    m_x=0.52,
    m_y=0.01,
    kde_y=0.4,
):
    unique_col_values = df[col].unique()
    num_cols = len(unique_col_values)
    num_rows = 1
    unique_row_values = [None]

    if row is not None and row in df.columns:
        # unique_row_values = df[row].unique()
        unique_row_values = sorted(df[row].unique())[::-1]
        num_rows = len(unique_row_values)

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
            fr = sub_df["p < 0.05"].mean()

            ax_scatter = row_axes[col_i] if num_cols > 1 else axes
            if aspect is not None:
                ax_scatter.set(aspect=aspect)
            ax_kde = ax_scatter.inset_axes([0, 1, 1, 0.15], sharex=ax_scatter)

            sns.scatterplot(
                ax=ax_scatter,
                data=sub_df,
                x=metric,
                y=f"-log10({metric} p-value)",
                hue="p < 0.05",
                s=50,
            )
            ax_scatter.set_xlabel("")
            ax_scatter.set_xlim(-0.05, 1.05)
            ax_scatter.set_ylim(-0.1, max(sub_df[f"-log10({metric} p-value)"]) + 0.25)
            ax_scatter.xaxis.set_major_locator(MultipleLocator(base=0.25))
            ax_scatter.get_xaxis().set_major_formatter(
                FuncFormatter(partial(x_axis_formatter, n_labels=6))
            )

            if col_i == 1:
                ax_scatter.text(
                    pr_x, 0.02, f"Retrieved: {fr:.0%}", transform=ax_scatter.transAxes
                )
            else:
                ax_scatter.text(
                    pr_x, pr_y, f"Retrieved: {fr:.0%}", transform=ax_scatter.transAxes
                )

            # if col_i == num_cols - 1 and row_i == 0:
            #     sns.move_legend(ax_scatter, "upper left", bbox_to_anchor=(1., .5), frameon=False)
            # else:
            #     ax_scatter.get_legend().remove()
            handles, labels = ax_scatter.get_legend_handles_labels()
            ax_scatter.get_legend().remove()

            ax_scatter.set_title(f"{col_value}", fontsize=16, pad=20)

            max_kde_y = 0
            for p_value in sorted(sub_df["p < 0.05"].unique()):
                sns.kdeplot(
                    ax=ax_kde,
                    data=sub_df[sub_df["p < 0.05"] == p_value],
                    x=metric,
                    label=str(p_value),
                )
                max_kde_y = max(max_kde_y, max(ax_kde.lines[-1].get_ydata()))

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

            ax_kde.get_yaxis().set_visible(False)
            ax_kde.get_xaxis().set_visible(False)
            sns.despine(ax=ax_kde, left=True, bottom=True)

    fig.text(m_x, m_y, metric, ha="center", va="center")

    plt.tight_layout()
    fig.legend(
        handles,
        labels,
        title="p < 0.05",
        loc="upper center",
        bbox_to_anchor=(l_x, l_y),
        frameon=False,
    )
    if adjust is not None:
        fig.subplots_adjust(**adjust)
    plt.show()


def plot_map_x3_hue(
    df,
    col,
    title,
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
    figure=None,
    save_path=None,
):
    unique_col_values = df[col].unique()
    num_cols = len(unique_col_values)
    num_rows = 1
    unique_row_values = [None]

    hue_col = hue_col if hue_col is not None else "p < 0.05"
    palette = palette if palette is not None else "Set1"

    if row is not None and row in df.columns:
        # unique_row_values = df[row].unique()
        unique_row_values = sorted(df[row].unique())[::-1]
        num_rows = len(unique_row_values)

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
            mean_map = sub_df["mAP"].mean()  # noqa: F841
            fr_total = sub_df["p < 0.05"].mean()
            frs = sub_df.groupby(hue_col)["p < 0.05"].mean()

            ax_scatter = row_axes[col_i] if num_cols > 1 else axes
            if aspect is not None:
                ax_scatter.set(aspect=aspect)
            ax_kde = ax_scatter.inset_axes([0, 1, 1, 0.15], sharex=ax_scatter)

            sns.scatterplot(
                ax=ax_scatter,
                data=sub_df,
                x="mAP",
                y="-log10(mAP p-value)",
                style="markers",
                # hue='p < 0.05',
                hue=hue_col,
                palette=palette,
                markers={"p<0.05": "o", "p>=0.05": "s"},
                s=50,
            )
            ax_scatter.set_xlabel("")
            ax_scatter.set_xlim(-0.05, 1.05)
            ax_scatter.set_ylim(-0.1, max(sub_df["-log10(mAP p-value)"]) + 0.25)
            ax_scatter.xaxis.set_major_locator(MultipleLocator(base=0.25))
            ax_scatter.get_xaxis().set_major_formatter(
                FuncFormatter(partial(x_axis_formatter, n_labels=6))
            )

            ax_scatter.axhline(-np.log10(0.05), color="grey", linestyle="--")

            # Specific for Fig 3A
            if y_label is not None and figure == "Fig3A":
                ax_scatter.set_ylabel(
                    f"Preprocessing: {row_value}\n-log10(mAP p-value)"
                )

            # if col_i == 1:
            #     ax_scatter.text(pr_x-0.2, 0.02, f"Retrieved: ", transform=ax_scatter.transAxes)
            #     for i, (hue_value, fr) in enumerate(frs.items()):
            #         ax_scatter.text(pr_x+i*0.1, 0.02, f"{fr:.0%}", transform=ax_scatter.transAxes, color=palette[hue_value])
            # else:

            # Fig 3A
            if figure == "Fig3A":
                ax_scatter.text(
                    pr_x - 0.33, pr_y, "Retrieved: ", transform=ax_scatter.transAxes
                )
                ax_scatter.text(
                    0.82,
                    0.15,
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

            # Fig 3B
            if figure in ["Fig3B"]:
                ax_scatter.text(
                    pr_x, pr_y + 0.1, "Retrieved:", transform=ax_scatter.transAxes
                )
                ax_scatter.text(
                    0.8,
                    pr_y - 0.15,
                    "p=0.05",
                    transform=ax_scatter.transAxes,
                    color="grey",
                    fontsize=12,
                    fontstyle="italic",
                )
                for i, (hue_value, fr) in enumerate(frs.items()):
                    ax_scatter.text(
                        (pr_x - 0.2) + i * 0.16,
                        pr_y,
                        f"{fr:.0%}",
                        transform=ax_scatter.transAxes,
                        color=palette[hue_value],
                    )
                ax_scatter.text(
                    (pr_x - 0.2) + len(frs) * 0.16,
                    pr_y,
                    f"({fr_total:.0%})",
                    transform=ax_scatter.transAxes,
                )

            # Fig 3E
            if figure in ["Fig3E"]:
                ax_scatter.text(
                    pr_x, pr_y + 0.1, "Retrieved:", transform=ax_scatter.transAxes
                )
                ax_scatter.text(
                    0.82,
                    pr_y - 0.4,
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

            # if col_i == num_cols - 1 and row_i == 0:
            #     sns.move_legend(ax_scatter, "upper left", bbox_to_anchor=(1., .5), frameon=False)
            # else:
            #     ax_scatter.get_legend().remove()
            handles, labels = ax_scatter.get_legend_handles_labels()
            labels = [label.replace("markers", "") for label in labels]
            ax_scatter.get_legend().remove()

            mmaps = {}
            max_kde_y = 0
            for hue_val in sorted(sub_df[hue_col].unique()):
                sns.kdeplot(
                    ax=ax_kde,
                    data=sub_df[(sub_df[hue_col] == hue_val)],
                    x="mAP",
                    label=str(hue_val),
                    cut=1,
                    color=palette[hue_val],
                )
                mmaps[hue_val] = sub_df[(sub_df[hue_col] == hue_val)]["mAP"].mean()
                # sns.kdeplot(
                #     ax=ax_kde,
                #     data=sub_df[(sub_df[hue_col] == hue_val) & (sub_df['p < 0.05'] == False)],
                #     x='mAP',
                #     label=str(hue_val),
                #     linestyle="--"
                # )
                max_kde_y = max(max_kde_y, max(ax_kde.lines[-1].get_ydata()))

            # ax_kde.axvline(mean_map, color='grey')
            # for hue_value, mmap in mmaps.items():
            #     ax_kde.axvline(mmap, color=palette[hue_value])

            # if figure == "Fig3A":

            #     ax_scatter.text(pr_x-0.25, 0.02, f"mmAP: ", transform=ax_scatter.transAxes)
            #     for i, (hue_value, mmap) in enumerate(mmaps.items()):
            #         ax_scatter.text(pr_x+i*0.15, 0.02, f"{mmap:.2f}", transform=ax_scatter.transAxes, color=palette[hue_value],)
            #     ax_scatter.text(pr_x+len(mmaps)*0.15, 0.02, f"({mean_map:.2f})", transform=ax_scatter.transAxes)

            # else:
            #     if mean_map < 0.5:
            #         ax_kde.text(mean_map + 0.1, 0.4, f"Mean mAP: {mean_map:.2f}", transform=ax_kde.transAxes)
            #     else:
            #         ax_kde.text(mean_map -0.22, 0.4, f"Mean mAP: {mean_map:.2f}", transform=ax_kde.transAxes)

            ax_kde.get_yaxis().set_visible(False)
            ax_kde.get_xaxis().set_visible(False)
            sns.despine(ax=ax_kde, left=True, bottom=True)

    fig.text(m_x, m_y, "mAP", ha="center", va="center")

    plt.tight_layout()
    fig.legend(
        handles,
        labels,
        title="",
        loc="upper center",
        bbox_to_anchor=(l_x, l_y),
        frameon=False,
    )
    if adjust is not None:
        fig.subplots_adjust(**adjust)
    # fig.subplots_adjust(right=0.85)

    if save_path is not None:
        fig.savefig(Path(save_path) / f"{figure}.png", dpi=300)
        fig.savefig(Path(save_path) / f"{figure}.svg")

    plt.show()
