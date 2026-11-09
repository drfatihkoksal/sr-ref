"""Figures 1-4 from results/paper_numbers.json (aggregates only), as PDF and PNG.

Figure 1: records to pairs by source. Figure 2: discordance of the V2TR and TZ index calls (Wilson 95% CI).
Figure 3: V2TR discordance by time between the 2 ECGs (exploratory). Figure 4: paired differences between the
arms of the call comparison. At least 7 inches wide, 10-12 pt lettering, dot plots with 95% CIs.

    python -m src.figures   -> results/figures/Figure1-4.pdf, .png
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figures"
N = json.loads((ROOT / "results" / "paper_numbers.json").read_text())
S = ("mgh", "euh", "mimic")
NAME = {"mgh": "MGH", "euh": "Emory", "mimic": "MIMIC-IV-ECG"}
LONG = {"mgh": "MGH", "euh": "Emory", "mimic": "MIMIC-IV-ECG\n(preregistered)"}
COLOR = {"mgh": "#2a78d6", "euh": "#eb6834", "mimic": "#1baf7a"}
MARK = {"mgh": "o", "euh": "s", "mimic": "^"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
W = 7.5

plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 10.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.8, "pdf.fonttype": 42,
                     "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.facecolor": "white"})


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)


def style(ax, grid_axis="x"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis=grid_axis, color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)


def legend(ax, **kw):
    handles = [plt.Line2D([], [], color=COLOR[s], marker=MARK[s], linestyle="", markersize=8, label=NAME[s]) for s in S]
    ax.legend(handles=handles, frameon=False, **kw)


def figure1():
    f = N["flow"]
    rows = [("Records processed", lambda s: f"{f[s]['records']:,} records\n{f[s]['patients']:,} patients"),
            ("Records excluded", None),
            ("Eligible records with an\noutflow tract-like PVC", lambda s: f"{f[s]['ot_like']:,}"),
            ("Patients with ≥2 such records", lambda s: f"{f[s]['patients_2ot']:,}"),
            ("Primary pairs: matched PVC,\n≤365 days, 1 per patient", lambda s: f"{N['pairs'][s]['primary']:,}"),
            ("Both ECGs in the\npublished scope", lambda s: f"V2TR {N['discordance'][s]['v2tr']['both']['n']:,}\n"
                                                           f"TZ index {N['discordance'][s]['tzi']['both']['n']:,}")]
    heights = [0.62, 1.22, 0.45, 0.45, 0.45, 0.62]
    gap = 0.32
    total = sum(heights) + gap * (len(rows) - 1) + 0.55
    fig = plt.figure(figsize=(W, total))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, total)
    ax.axis("off")
    x0, cw = 1.95, (W - 2.0) / 3
    excl = (("Lead quality", "signal"), ("No SR template", "no_sr"), ("No usable PVC", "no_pvc"), ("Multiform PVCs", "multiform"),
            ("PVC not OT-like", "not_ot"))
    top = total - 0.1
    for k, s in enumerate(S):
        cx = x0 + cw * (k + 0.5)
        ax.text(cx, top - 0.12, LONG[s].replace("\n", " "), ha="center", va="center", fontsize=11, fontweight="bold", color=INK)
    y = top - 0.45
    centers = []
    for (label, txt), h in zip(rows, heights):
        excluded = label.startswith("Records excluded")
        ax.text(0.05, y - h / 2, label, ha="left", va="center", fontsize=10, color=INK2 if excluded else INK)
        for k, s in enumerate(S):
            cx = x0 + cw * (k + 0.5)
            bw = cw - 0.16
            box = FancyBboxPatch((cx - bw / 2, y - h), bw, h, boxstyle="round,pad=0,rounding_size=0.06",
                                 facecolor="#f4f4f1" if excluded else "white", edgecolor=MUTED if excluded else INK2, linewidth=0.9)
            ax.add_patch(box)
            if excluded:
                for i, (lab, key) in enumerate(excl):
                    yy = y - 0.16 - i * 0.22
                    ax.text(cx - bw / 2 + 0.07, yy, lab, ha="left", va="center", fontsize=10, color=INK2)
                    ax.text(cx + bw / 2 - 0.07, yy, f"{f[s][key]:,}", ha="right", va="center", fontsize=10, color=INK2)
            else:
                ax.text(cx, y - h / 2, txt(s), ha="center", va="center", fontsize=10.5, color=INK, linespacing=1.35)
        centers.append((y, h))
        y -= h + gap
    for (ya, ha), (yb, hb) in zip(centers[:-1], centers[1:]):
        for k in range(3):
            cx = x0 + cw * (k + 0.5)
            ax.add_patch(FancyArrowPatch((cx, ya - ha), (cx, yb), arrowstyle="-|>", mutation_scale=10, color=INK2, linewidth=0.9))
    save(fig, "Figure1")


def dots(ax, rows, key_f, xmax):
    for i, (s, cell) in enumerate(rows):
        f, lo, hi = 100 * cell["f"], 100 * cell["lo"], 100 * cell["hi"]
        ax.plot([lo, hi], [i, i], color=COLOR[s], linewidth=2, solid_capstyle="round")
        ax.plot(f, i, marker=MARK[s], color=COLOR[s], markersize=9, markeredgecolor="white", markeredgewidth=1.2, linestyle="")
        ax.text(hi + 1.2, i, f"{f:.1f}%", va="center", ha="left", fontsize=10, color=INK,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.0})
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([LONG[s] for s, _ in rows])
    ax.set_ylim(len(rows) - 0.5, -1.0)
    ax.set_xlim(0, xmax)


def figure2():
    fig, axes = plt.subplots(1, 2, figsize=(W, 3.5), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, (c, title) in zip(axes, (("v2tr", "V2 transition ratio"), ("tzi", "Transition zone index"))):
        style(ax)
        for ref, lab, side in ((10, "95%", "right"), (18, "91%", "left")):
            ax.axvline(ref, color=INK2, linestyle=(0, (4, 3)), linewidth=0.9)
            ax.text(ref + (0.8 if side == "left" else -0.8), -0.75, lab, fontsize=10, color=INK2, ha=side, va="center", linespacing=1.1)
        dots(ax, [(s, N["discordance"][s][c]["both"]) for s in S], None, 52)
        ax.set_title(title, fontsize=11, color=INK, loc="left", pad=6)
        ax.set_xlabel("Pairs whose call differed, % (95% CI)")
    save(fig, "Figure2")


def figure3():
    bins = (("lt1h", "<1 hour"), ("h1to24", "1-24 hours"), ("d1to30", "1-30 days"), ("d31to365", "31-365 days"))
    fig, ax = plt.subplots(figsize=(8.4, 3.7))
    style(ax, "y")
    for k, s in enumerate(S):
        off = (k - 1) * 0.2
        for j, (b, _) in enumerate(bins):
            c = N["interval"][s][b]
            ax.plot([j + off] * 2, [100 * c["lo"], 100 * c["hi"]], color=COLOR[s], linewidth=2, solid_capstyle="round")
            ax.plot(j + off, 100 * c["f"], marker=MARK[s], color=COLOR[s], markersize=9, markeredgecolor="white", markeredgewidth=1.2)
    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels([f"{lab}\n(n = " + ", ".join(str(N["interval"][s][b]["n"]) for s in S) + ")" for b, lab in bins])
    ax.set_ylim(0, 70)
    ax.set_ylabel("V2TR call differed, % (95% CI)")
    ax.set_xlabel("Time between the 2 ECGs (pairs: MGH, Emory, MIMIC-IV-ECG)")
    legend(ax, loc="upper left", ncols=3, handletextpad=0.3, columnspacing=1.2)
    save(fig, "Figure3")


def figure4():
    comps = (("raw_minus_shared", "Own SR for each ECG\nminus shared SR"), ("own_minus_shared", "Own SR-derived transform\nminus shared SR"),
             ("own_minus_donor", "Own transform\nminus donor transform"))
    fig, axes = plt.subplots(1, 2, figsize=(W, 3.9), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, (c, title) in zip(axes, (("v2tr", "V2 transition ratio"), ("tzi", "Transition zone index"))):
        style(ax)
        ax.axvline(0, color=INK2, linewidth=0.9)
        for j, (key, _) in enumerate(comps):
            for k, s in enumerate(S):
                y = j + (k - 1) * 0.22
                d = N["arms"][s][c][key]
                ax.plot([100 * d["lo"], 100 * d["hi"]], [y, y], color=COLOR[s], linewidth=2, solid_capstyle="round")
                ax.plot(100 * d["m"], y, marker=MARK[s], color=COLOR[s], markersize=8, markeredgecolor="white", markeredgewidth=1.2)
        ax.set_yticks(range(len(comps)))
        ax.set_yticklabels([lab for _, lab in comps])
        ax.set_ylim(len(comps) - 0.5, -0.5)
        ax.set_xlim(-15, 15)
        ax.set_xticks([-10, -5, 0, 5, 10])
        ax.set_title(title, fontsize=11, color=INK, loc="left", pad=6)
        ax.set_xlabel("Difference in discordance,\npercentage points (95% CI)")
    legend(axes[0], loc="lower left", bbox_to_anchor=(0, -0.36), ncols=3, handletextpad=0.3, columnspacing=1.2)
    save(fig, "Figure4")


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    figure4()
    print("figures:", ", ".join(sorted(p.name for p in OUT.glob("*.png"))))
