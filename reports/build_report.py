#!/usr/bin/env python3
"""Builds the AIIPS progress report (figures with matplotlib, PDF with reportlab).

Usage: build_report.py REPO OUT_PDF [BUILD_DIR]

Every table and figure is computed from the result files in REPO/benchmarks.
Figures are written to BUILD_DIR (default: a temporary directory). Needs
matplotlib, numpy, Pillow and reportlab (not dependencies of aiips itself),
and the DejaVu fonts in /usr/share/fonts/truetype/dejavu.
"""
import csv
import io
import os
import re
import statistics as st
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (CondPageBreak, Image, KeepTogether, PageBreak, Paragraph,
                                SimpleDocTemplate, Spacer, Table, TableStyle)
from reportlab.platypus.tableofcontents import TableOfContents
from reportlab.lib.fonts import addMapping

REPO, OUT = sys.argv[1:3]
FIG = sys.argv[3] if len(sys.argv) > 3 else __import__("tempfile").mkdtemp(prefix="aiips-report-")
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------- fonts
FD = "/usr/share/fonts/truetype/dejavu/"
for name, file in [("Serif", "DejaVuSerif.ttf"), ("Serif-Bold", "DejaVuSerif-Bold.ttf"),
                   ("Serif-Italic", "DejaVuSerif-Italic.ttf"), ("Serif-BoldItalic", "DejaVuSerif-BoldItalic.ttf"),
                   ("Sans", "DejaVuSans.ttf"), ("Sans-Bold", "DejaVuSans-Bold.ttf"),
                   ("Sans-Oblique", "DejaVuSans-Oblique.ttf"), ("Sans-BoldOblique", "DejaVuSans-BoldOblique.ttf"),
                   ("Mono", "DejaVuSansMono.ttf"), ("Mono-Bold", "DejaVuSansMono-Bold.ttf")]:
    pdfmetrics.registerFont(TTFont(name, FD + file))
for fam, r, b, i, bi in [("Serif", "Serif", "Serif-Bold", "Serif-Italic", "Serif-BoldItalic"),
                         ("Sans", "Sans", "Sans-Bold", "Sans-Oblique", "Sans-BoldOblique")]:
    addMapping(fam, 0, 0, r); addMapping(fam, 1, 0, b); addMapping(fam, 0, 1, i); addMapping(fam, 1, 1, bi)
addMapping("Mono", 0, 0, "Mono"); addMapping("Mono", 1, 0, "Mono-Bold")
addMapping("Mono", 0, 1, "Mono"); addMapping("Mono", 1, 1, "Mono-Bold")

ACCENT = colors.HexColor("#1f4e79")
styles = {
    "title": ParagraphStyle("title", fontName="Sans-Bold", fontSize=20, leading=25, alignment=TA_CENTER, textColor=ACCENT, spaceAfter=6),
    "subtitle": ParagraphStyle("subtitle", fontName="Sans", fontSize=11, leading=15, alignment=TA_CENTER, textColor=colors.HexColor("#444444")),
    "h1": ParagraphStyle("h1", fontName="Sans-Bold", fontSize=14.5, leading=19, spaceBefore=14, spaceAfter=6, textColor=ACCENT, keepWithNext=1),
    "h2": ParagraphStyle("h2", fontName="Sans-Bold", fontSize=11.5, leading=15, spaceBefore=10, spaceAfter=4, textColor=ACCENT, keepWithNext=1),
    "h3": ParagraphStyle("h3", fontName="Sans-Bold", fontSize=10, leading=13, spaceBefore=7, spaceAfter=2, keepWithNext=1),
    "body": ParagraphStyle("body", fontName="Serif", fontSize=9.6, leading=13.4, alignment=TA_JUSTIFY, spaceAfter=5),
    "bullet": ParagraphStyle("bullet", fontName="Serif", fontSize=9.6, leading=13.2, alignment=TA_LEFT, leftIndent=14, bulletIndent=4, spaceAfter=2.5),
    "caption": ParagraphStyle("caption", fontName="Serif-Italic", fontSize=8.4, leading=11, alignment=TA_JUSTIFY, spaceBefore=3, spaceAfter=10, textColor=colors.HexColor("#333333")),
    "cell": ParagraphStyle("cell", fontName="Serif", fontSize=7.6, leading=9.4, alignment=TA_LEFT),
    "cellc": ParagraphStyle("cellc", fontName="Serif", fontSize=7.6, leading=9.4, alignment=TA_CENTER),
    "cellh": ParagraphStyle("cellh", fontName="Sans-Bold", fontSize=7.4, leading=9.2, alignment=TA_CENTER, textColor=colors.white),
    "code": ParagraphStyle("code", fontName="Mono", fontSize=7.8, leading=10, leftIndent=10, backColor=colors.HexColor("#f3f3f3"), borderPadding=4, spaceBefore=3, spaceAfter=8),
    "ref": ParagraphStyle("ref", fontName="Serif", fontSize=8.6, leading=11, leftIndent=18, firstLineIndent=-18, spaceAfter=2),
    "toc1": ParagraphStyle("toc1", fontName="Sans-Bold", fontSize=9.5, leading=12.5, leftIndent=10, spaceBefore=2),
    "toc2": ParagraphStyle("toc2", fontName="Sans", fontSize=8.6, leading=10.8, leftIndent=26, textColor=colors.HexColor("#333333")),
}

story = []
counters = {"fig": 0, "tab": 0}


class Heading(Paragraph):
    pass


def h1(text):
    story.append(CondPageBreak(2.0 * inch))
    p = Heading(text, styles["h1"]); p._toc = (0, text); story.append(p)


def h2(text):
    story.append(CondPageBreak(1.5 * inch))
    story.append(Paragraph(text, styles["h2"]))


def h3(text):
    story.append(Paragraph(text, styles["h3"]))


def para(text):
    story.append(Paragraph(text, styles["body"]))


def bullets(items):
    for item in items:
        story.append(Paragraph(item, styles["bullet"], bulletText="•"))
    story.append(Spacer(1, 3))


def code(text):
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\n", "<br/>").replace("  ", "&nbsp;&nbsp;")
    story.append(Paragraph(text, styles["code"]))


def equation(tex, size=11.5, number=None):
    """Renders a mathtext equation to a PNG and adds it centered."""
    fig = plt.figure(figsize=(0.01, 0.01))
    text = fig.text(0, 0, f"${tex}$", fontsize=size)
    buf = io.BytesIO()
    fig.savefig(buf, dpi=300, bbox_inches="tight", pad_inches=0.04, transparent=True)
    plt.close(fig)
    buf.seek(0)
    from PIL import Image as PILImage
    w, h = PILImage.open(io.BytesIO(buf.getvalue())).size
    buf.seek(0)
    img = Image(buf, width=w / 300 * 72, height=h / 300 * 72)
    if number is not None:
        row = [[img, Paragraph(f"({number})", styles["cellc"])]]
        t = Table(row, colWidths=[6.0 * inch, 0.5 * inch])
        t.setStyle(TableStyle([("ALIGN", (0, 0), (0, 0), "CENTER"), ("VALIGN", (0, 0), (-1, -1), "MIDDLE")]))
        story.append(t)
    else:
        img.hAlign = "CENTER"
        story.append(img)
    story.append(Spacer(1, 4))


eq_no = [0]


def eqn(tex, size=11.5):
    eq_no[0] += 1
    equation(tex, size, eq_no[0])
    return eq_no[0]


def figure(path, width_in, caption):
    counters["fig"] += 1
    from PIL import Image as PILImage
    w, h = PILImage.open(path).size
    img = Image(path, width=width_in * inch, height=width_in * inch * h / w)
    story.append(KeepTogether([img, Paragraph(f"<b>Figure {counters['fig']}.</b> {caption}", styles["caption"])]))
    return counters["fig"]


def table(rows, col_widths, caption, header_rows=1, zebra=True, align_center_from=1, spans=()):
    counters["tab"] += 1
    data = []
    for r_i, row in enumerate(rows):
        out = []
        for c_i, cell in enumerate(row):
            if r_i < header_rows:
                out.append(Paragraph(str(cell), styles["cellh"]))
            else:
                out.append(Paragraph(str(cell), styles["cellc"] if c_i >= align_center_from else styles["cell"]))
        data.append(out)
    t = Table(data, colWidths=[w * inch for w in col_widths], repeatRows=header_rows)
    style = [("BACKGROUND", (0, 0), (-1, header_rows - 1), ACCENT),
             ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
             ("LINEBELOW", (0, header_rows - 1), (-1, header_rows - 1), 0.6, ACCENT),
             ("LINEBELOW", (0, -1), (-1, -1), 0.6, ACCENT),
             ("TOPPADDING", (0, 0), (-1, -1), 2.2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2),
             ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3)]
    for (c0, r0, c1, r1) in spans:
        style.append(("SPAN", (c0, r0), (c1, r1)))
    if zebra:
        for i in range(header_rows, len(rows)):
            if (i - header_rows) % 2 == 1:
                style.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#eef3f8")))
    t.setStyle(TableStyle(style))
    cap = Paragraph(f"<b>Table {counters['tab']}.</b> {caption}", styles["caption"])
    story.append(KeepTogether([cap, t, Spacer(1, 10)]))
    return counters["tab"]


# ---------------------------------------------------------------- data
def load(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


R10 = load(os.path.join(REPO, "benchmarks/simulated-10mb.tsv"))
R5 = load(os.path.join(REPO, "benchmarks/simulated-5mb.tsv"))
GAP = load(os.path.join(REPO, "benchmarks/simulated-5mb-gap-sweep.tsv"))
KSWEEP = load(os.path.join(REPO, "benchmarks/simulated-10mb-ksweep.tsv"))
TRUTH10 = load(os.path.join(REPO, "benchmarks/simulated-10mb-truth.tsv"))

TOOLS = [("skani", "skani 0.4.0"), ("committed", "aiips, previous commit"), ("aiips", "aiips, this commit"),
         ("aiips:both", "aiips, both flags")]
COLORS = {"skani": "#7f7f7f", "committed": "#9ecae1", "aiips": "#2171b5", "aiips:both": "#e6550d"}


def index(rows):
    by = {}
    for r in rows:
        by.setdefault(r["scenario"], {})[r["tool"]] = r
    return by


BY10 = index(R10)


def summary(rows):
    have = [r for r in rows if r["ani"] != "NA"]
    ae = [100 * (float(r["ani"]) - float(r["true_ani"])) for r in have]
    af = [100 * abs(float(r[k]) - float(r["true_" + k])) for r in have for k in ("af_a", "af_b")]
    out = dict(n=len(rows), missing=len(rows) - len(have),
               ani=st.mean(map(abs, ae)) if ae else float("nan"),
               bias=st.mean(ae) if ae else float("nan"),
               ani_max=max(map(abs, ae)) if ae else float("nan"),
               af=st.mean(af) if af else float("nan"))
    for k in ("wall_s", "user_s", "sys_s"):
        if k in rows[0]:
            out[k] = st.mean(float(r[k]) for r in rows) * 1000
    out["rss"] = st.mean(int(r["max_rss_kb"]) / 1024 for r in rows)
    return out


def group_of(s):
    a = float(BY10[s]["aiips"]["true_ani"])
    if s.startswith("sub_"):
        return ("Substitutions", "ANI ≥ 99%" if a >= 0.989 else ("ANI 90–97%" if a >= 0.895 else "ANI 80–85%"))
    if s.startswith("partial_"):
        return ("Partial overlap", f"{float(re.search(r'f([0-9.]+)', s).group(1)) * 100:.0f}% of A shared")
    if s.startswith("draft_"):
        n = int(re.search(r"_n(\d+)", s).group(1))
        return ("Indels + draft", {0: "one record", 50000: "contigs ≈ 50 kb", 10000: "contigs ≈ 10 kb", 3000: "contigs ≈ 3 kb"}[n])
    return ("Repeats + rearr.", "all three")


GROUP_ORDER = [("Substitutions", "ANI ≥ 99%"), ("Substitutions", "ANI 90–97%"), ("Substitutions", "ANI 80–85%"),
               ("Partial overlap", "20% of A shared"), ("Partial overlap", "50% of A shared"), ("Partial overlap", "90% of A shared"),
               ("Indels + draft", "one record"), ("Indels + draft", "contigs ≈ 50 kb"), ("Indels + draft", "contigs ≈ 10 kb"),
               ("Indels + draft", "contigs ≈ 3 kb"), ("Repeats + rearr.", "all three")]
GROUPS = {}
for s in BY10:
    GROUPS.setdefault(group_of(s), []).append(s)


def gsum(tool, g):
    return summary([BY10[s][tool] for s in (GROUPS[g] if g != "ALL" else list(BY10))])


ALL = {t: gsum(t, "ALL") for t, _ in TOOLS}
# Underestimation range on indel scenarios at 90-95% ANI, per tool (points).
IND = {}
for t, _ in TOOLS:
    errs = [100 * (float(float(BY10[s][t]["true_ani"]) - float(BY10[s][t]["ani"]))) for s in BY10
            if s.startswith("draft_") and float(BY10[s][t]["true_ani"]) < 0.97 and BY10[s][t]["ani"] != "NA"]
    IND[t] = (min(errs), max(errs))

# ---------------------------------------------------------------- figures
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def fig_schematic():
    fig, ax = plt.subplots(figsize=(7.2, 2.9))
    ax.set_xlim(-1, 52); ax.set_ylim(-0.6, 3.4); ax.axis("off")
    ya, yb = 2.4, 0.6
    for y, name in [(ya, "genome A"), (yb, "genome B")]:
        ax.plot([0, 50], [y, y], color="#555555", lw=1.2, zorder=1)
        ax.text(-0.6, y, name, ha="right", va="center", fontsize=8.5)
    a_anchors = [3, 13, 23, 33, 43]
    b_anchors = [3, 13, 33, 43, 28]  # 23 mutated away in B, 28 created by a mutation (B only)
    rng = np.random.default_rng(3)
    probes_a = sorted(set(np.round(rng.uniform(0.3, 49.5, 22), 1)) - set())
    probes_a = [p for p in probes_a if all(abs(p - a - 0.5) > 1.1 for a in a_anchors)]
    shared = [p for i, p in enumerate(probes_a) if i % 2 == 0 and not (22 < p < 25)]
    probes_b = shared + [p + 0.4 for i, p in enumerate(probes_a) if i % 4 == 1]
    probes_b = [p for p in probes_b if all(abs(p - a - 0.5) > 1.1 for a in b_anchors)]
    probes_b += [23.6]  # B keeps a probe where A had an anchor
    for x in a_anchors:
        ax.add_patch(plt.Rectangle((x, ya - 0.17), 1.0, 0.34, color="#2171b5", zorder=3))
    for x in b_anchors:
        col = "#2171b5" if x != 28 else "#9ecae1"
        ax.add_patch(plt.Rectangle((x, yb - 0.17), 1.0, 0.34, color=col, zorder=3))
    for p in probes_a:
        ax.plot([p, p], [ya - 0.12, ya + 0.12], color="#e6550d", lw=1.6, zorder=2)
    for p in probes_b:
        ax.plot([p, p], [yb - 0.12, yb + 0.12], color="#e6550d", lw=1.6, zorder=2)
    for p in shared:
        if p in probes_b:
            ax.plot([p, p], [ya - 0.14, yb + 0.14], color="#e6550d", lw=0.5, alpha=0.45, zorder=1)
    for x in [3, 13, 33, 43]:
        ax.annotate("", xy=(x + 0.5, yb + 0.2), xytext=(x + 0.5, ya - 0.2),
                    arrowprops=dict(arrowstyle="-", color="#2171b5", lw=1.2, ls="--"))
    ax.text(23.5, yb - 0.45, "anchor mutated away: the islands of A and B\nno longer coincide (one island pairs with several)", ha="center", va="top", fontsize=7)
    ax.plot(23.5, yb, marker="x", color="black", ms=7, zorder=4)
    ax.text(28.5, yb + 0.45, "B-only anchor", ha="center", fontsize=7, color="#3182bd")
    for i, (l, r) in enumerate(zip([4] + [14, 24, 34, 44], [13, 23, 33, 43, 50])):
        ax.text((l + r) / 2, ya + 0.42, f"island {i + 1}", ha="center", fontsize=7, color="#444444")
    ax.annotate("", xy=(49, 3.15), xytext=(4, 3.15), arrowprops=dict(arrowstyle="->", color="#31a354", lw=1.6))
    ax.text(26, 3.25, "walk from the seed anchor pair at the left", ha="center", fontsize=7.5, color="#31a354")
    ax.add_patch(plt.Rectangle((0.6, yb - 0.17), 0.7, 0.34, color="#31a354", zorder=3))
    ax.text(0.9, yb - 0.6, "promoted\nanchor", ha="center", va="top", fontsize=6.5, color="#31a354")
    handles = [plt.Rectangle((0, 0), 1, 1, color="#2171b5"), plt.Line2D([0], [0], color="#e6550d", lw=2),
               plt.Line2D([0], [0], color="#2171b5", ls="--"), plt.Line2D([0], [0], color="#e6550d", lw=0.6)]
    ax.legend(handles, ["anchor (k1-mer, sparse)", "probe (k2-mer, dense)", "matched anchors", "shared probe"],
              loc="lower right", fontsize=6.8, frameon=False, ncol=2, bbox_to_anchor=(1.0, -0.2))
    path = os.path.join(FIG, "schematic.png"); fig.savefig(path, dpi=220, bbox_inches="tight"); plt.close(fig)
    return path


def fig_likelihood():
    K1, N1, n1 = 21, 39640, 13410
    K2, N2, n2 = 15, 186093, 83194
    def ll(q, k, n, N):
        return n * k * np.log(q) + (N - n) * np.log1p(-q ** k)
    q = np.linspace(0.93, 0.965, 3000)
    la, lp = ll(q, K1, n1, N1), ll(q, K2, n2, N2)
    lj = la + lp
    qf = np.linspace(0.93, 0.965, 700001)
    mle = qf[np.argmax(ll(qf, K1, n1, N1) + ll(qf, K2, n2, N2))]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7))
    ax = axes[0]
    ax.plot(q, la - la.max(), label=f"anchors (k1={K1})", color="#2171b5")
    ax.plot(q, lp - lp.max(), label=f"probes (k2={K2})", color="#e6550d")
    ax.plot(q, lj - lj.max(), label="joint", color="#31a354", lw=2.2)
    ax.axvline(0.949949, color="grey", ls=":", lw=1, label="truth")
    ax.axhline(-1.92, color="red", lw=0.7)
    ax.set_xlim(0.94, 0.955); ax.set_ylim(-40, 2)
    ax.set_xticks([0.940, 0.945, 0.950, 0.955])
    ax.set_xlabel("per-base identity q"); ax.set_ylabel("log L(q) − max")
    ax.legend(fontsize=6.2, frameon=False, loc="lower left")
    ax = axes[1]
    w = (q > 0.944) & (q < 0.952)
    ax.plot(q[w], lj[w] - lj.max(), color="#31a354", lw=2.2)
    ci = q[lj >= lj.max() - 1.92]
    ax.axvspan(ci.min(), ci.max(), color="#31a354", alpha=0.15, label="Δlog L ≥ −1.92")
    ax.axvline(mle, color="black", ls="--", lw=1, label="MLE")
    ax.axvline(0.949949, color="grey", ls=":", lw=1, label="truth")
    ax.set_ylim(-60, 3); ax.set_xlabel("per-base identity q")
    ax.legend(fontsize=6.2, frameon=False, loc="upper left")
    fig.tight_layout()
    path = os.path.join(FIG, "likelihood.png"); fig.savefig(path, dpi=220, bbox_inches="tight"); plt.close(fig)
    return path, mle, (ci.min(), ci.max())


def fig_sweep():
    subs = sorted([s for s in BY10 if s.startswith("sub_")], key=lambda s: -float(BY10[s]["aiips"]["true_ani"]))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    for t, label in TOOLS:
        xs, ys, afs = [], [], []
        for s in subs:
            r = BY10[s][t]
            if r["ani"] == "NA":
                continue
            xs.append(100 * float(r["true_ani"])); ys.append(100 * (float(r["ani"]) - float(r["true_ani"])))
            afs.append(float(r["af_a"]))
        style = dict(marker="o", ms=3.5, lw=1.3, color=COLORS[t], label=label)
        if t == "committed":
            style.update(ls="--", marker="s", ms=3)
        axes[0].plot(xs, ys, **style); axes[1].plot(xs, afs, **style)
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set_xlabel("true ANI (%)"); axes[0].set_ylabel("estimated − true ANI (points)")
    axes[1].set_xlabel("true ANI (%)"); axes[1].set_ylabel("estimated AF of A (true 1.0)")
    for ax in axes:
        ax.invert_xaxis()
    axes[1].legend(fontsize=6.6, frameon=False, loc="lower left")
    fig.tight_layout()
    path = os.path.join(FIG, "sweep.png"); fig.savefig(path, dpi=220, bbox_inches="tight"); plt.close(fig)
    return path


def fig_groups():
    labels = [f"{k}\n{c}" for k, c in GROUP_ORDER]
    short = ["Sub ≥99", "Sub 90–97", "Sub 80–85", "Part 20%", "Part 50%", "Part 90%",
             "Draft 1 rec", "Draft 50 kb", "Draft 10 kb", "Draft 3 kb", "Rep+rearr"]
    tools = [("skani", "skani"), ("aiips", "aiips"), ("aiips:both", "aiips, both flags")]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 3.5), sharex=True)
    x = np.arange(len(GROUP_ORDER)); w = 0.27
    for i, (t, lab) in enumerate(tools):
        ani = [gsum(t, g)["ani"] for g in GROUP_ORDER]
        af = [gsum(t, g)["af"] for g in GROUP_ORDER]
        axes[0].bar(x + (i - 1) * w, ani, w, color=COLORS[t], label=lab)
        axes[1].bar(x + (i - 1) * w, af, w, color=COLORS[t], label=lab)
    axes[0].set_ylabel("ANI mean abs. error\n(points)"); axes[1].set_ylabel("AF mean abs. error\n(points)")
    axes[1].set_xticks(x); axes[1].set_xticklabels(short, rotation=35, ha="right", fontsize=7.5)
    axes[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    path = os.path.join(FIG, "groups.png"); fig.savefig(path, dpi=220, bbox_inches="tight"); plt.close(fig)
    return path


def fig_cost():
    cats = [("sub", "Substitutions"), ("partial", "Partial overlap"), ("draft", "Indels + drafts"), ("rearr", "Repeats + rearr."), ("ALL", "All 32")]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7))
    x = np.arange(len(cats)); w = 0.2
    for i, (t, lab) in enumerate(TOOLS):
        walls, rss = [], []
        for c, _ in cats:
            rows = [r for r in R10 if r["tool"] == t and (c == "ALL" or r["category"] == c)]
            walls.append(st.mean(float(r["wall_s"]) for r in rows) * 1000)
            rss.append(st.mean(int(r["max_rss_kb"]) / 1024 for r in rows))
        axes[0].bar(x + (i - 1.5) * w, walls, w, color=COLORS[t], label=lab)
        axes[1].bar(x + (i - 1.5) * w, rss, w, color=COLORS[t], label=lab)
    for ax, yl in zip(axes, ["wall time per pair (ms)", "peak memory (MiB)"]):
        ax.set_xticks(x); ax.set_xticklabels([c[1] for c in cats], rotation=25, ha="right", fontsize=7.5); ax.set_ylabel(yl)
    fig.legend(*axes[0].get_legend_handles_labels(), fontsize=7, frameon=False, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    path = os.path.join(FIG, "cost.png"); fig.savefig(path, dpi=220, bbox_inches="tight"); plt.close(fig)
    return path


f_schematic = fig_schematic()
f_like, mle, ci = fig_likelihood()
f_sweep = fig_sweep()
f_groups = fig_groups()
f_cost = fig_cost()


# ---------------------------------------------------------------- document
def fmt(x, d=2):
    return "—" if x != x else f"{x:.{d}f}"


def sgn(x, d=2):
    if x != x:
        return "—"
    return f"{0.0 if abs(x) < 0.5 * 10 ** -d else x:+.{d}f}"


story.append(Paragraph("AIIPS: Anchor-Indexed Island-Probe Sketching", styles["title"]))
story.append(Paragraph("Progress report — positional two-level FracMinHash sketches for ANI and alignment fraction", styles["subtitle"]))
story.append(Spacer(1, 6))
story.append(Paragraph("17 November 2026 · repository <font face='Mono'>aiips</font>, branch <font face='Mono'>aiips-v2</font>, commit <font face='Mono'>1ee7ba0</font>", styles["subtitle"]))
story.append(Spacer(1, 8))

toc = TableOfContents()
toc.dotsMinLevel = 1
toc.levelStyles = [styles["toc1"], styles["toc2"]]
story.append(Paragraph("Contents", ParagraphStyle("contents", parent=styles["h2"], keepWithNext=0)))
story.append(toc)
story.append(PageBreak())

# ================================================================ 1 Summary
h1("1. Summary")
para("This report describes AIIPS, a sketching method for estimating the <b>average nucleotide identity (ANI)</b> "
     "between two genomes and the <b>alignment fraction (AF)</b> of each, i.e. the fraction of its bases that have a "
     "homologous counterpart in the other. Classic MinHash and FracMinHash sketches estimate ANI from k-mer set overlap "
     "but discard positions, so they cannot tell a short, highly similar shared region from a long, diverged one and "
     "cannot estimate AF. Alignment-based tools (FastANI, skani) recover positions at a higher cost. AIIPS keeps positions "
     "inside a FracMinHash-style sketch with two levels of sampling:")
bullets([
    "<b>Anchors</b>: long k-mers (k<sub>1</sub> = 21) sampled sparsely (1 in 500) with one hash function. They are rare, "
    "specific, and serve as seeds to locate homologous regions.",
    "<b>Probes</b>: shorter k-mers (k<sub>2</sub> = 15) sampled densely (1 in 100) with an independent hash function, "
    "excluding any k-mer that overlaps an anchor. They lie in the <b>islands</b> between consecutive anchors.",
    "<b>Comparison</b>: starting from each matched anchor, a walk moves island by island, compares the probes of "
    "corresponding islands, and keeps extending while islands share probes. The aligned islands give AF directly; "
    "the matched and unmatched anchors and probes inside them give ANI through a joint maximum-likelihood estimate.",
])
para("The method is implemented in Rust as <font face='Mono'>aiips sketch|compare -a aiips</font> (one external "
     "dependency, needletail, for FASTA/FASTQ parsing) and was benchmarked against skani 0.4.0 on 32 simulated pairs of "
     "10 Mbp genomes with known ground truth, covering substitutions (ANI 80–99.9%), partial overlap, small indels with "
     "draft assemblies of varying fragmentation, and repeats with genome rearrangements. Main results (Section 6):")
s_, a_, b_ = ALL["skani"], ALL["aiips"], ALL["aiips:both"]
bullets([
    f"<b>Accuracy.</b> Mean absolute ANI error over the 32 pairs is {a_['ani']:.2f} percentage points for AIIPS with "
    f"default settings and {b_['ani']:.2f} with the two optional corrections (<font face='Mono'>--symmetric-probes "
    f"--bidirectional</font>), against {s_['ani']:.2f} for skani (which also reported no value for the 80% pair). Mean "
    f"absolute AF error is {a_['af']:.1f} / {b_['af']:.1f} points against {s_['af']:.1f} for skani; the difference is "
    "largest below 90% ANI, where skani's AF falls to about 0.4 for genomes that are entirely homologous.",
    f"<b>Cost.</b> Per genome pair (single thread, sketching included) AIIPS takes {a_['wall_s']:.0f} ms and "
    f"{a_['rss']:.1f} MiB of peak memory, against {s_['wall_s']:.0f} ms and {s_['rss']:.1f} MiB for skani; the previous "
    f"commit took {ALL['committed']['wall_s']:.0f} ms and {ALL['committed']['rss']:.1f} MiB.",
    f"<b>Weak spots.</b> Both tools underestimate ANI when small indels are frequent (at 90–95% ANI by "
    f"{IND['aiips:both'][0]:.1f}–{IND['aiips'][1]:.1f} points for AIIPS, {IND['skani'][0]:.1f}–{IND['skani'][1]:.1f} for "
    "skani), and both lose AF on very fragmented assemblies (3 kbp contigs) at 90% ANI.",
])
para("Section 2 gives the background and theory, Section 3 the method, Section 4 the implementation and the "
     "performance work, Sections 5–6 the experiments and results, and Section 7 the limitations and next steps.")

# ================================================================ 2 Theory
h1("2. Background and theory")
h2("2.1 k-mers and canonical k-mers")
para("A <i>k-mer</i> is a substring of length k of a DNA sequence over {A, C, G, T}. Because a genome may be sequenced "
     "or assembled from either strand, a k-mer and its reverse complement (reversed, with A↔T and C↔G) represent the "
     "same locus. With the 2-bit encoding A = 0, C = 1, G = 2, T = 3, a k-mer with k ≤ 32 fits in a 64-bit integer; its "
     "<i>canonical</i> form is the smaller of the codes of the k-mer and of its reverse complement, so a canonical k-mer "
     "is the same on both strands. Both codes are maintained in O(1) per base while scanning a sequence. Characters other "
     "than A/C/G/T (e.g. N) break the run of valid k-mers.")
h2("2.2 FracMinHash, containment and ANI")
para("<b>FracMinHash</b> [3, 4] keeps every distinct canonical k-mer whose hash falls below a fixed threshold. With a "
     "hash h that maps k-mers to (approximately uniform) 64-bit values and a scale factor S,")
eqn(r"\mathrm{FRAC}_S(A)=\{\,x\in A:\ h(x)\leq \lfloor (2^{64}-1)/S\rfloor \,\}")
para("so each k-mer is kept with probability about 1/S, independently of the others, and the sketch size is "
     "proportional to the genome size (unlike MinHash [1, 2], whose size is fixed). Because the threshold is the same for "
     "every genome, the sketch of the intersection is the intersection of the sketches, and the <b>containment</b> of A "
     "in B, C(A, B) = |A ∩ B| / |A|, is estimated by |FRAC(A) ∩ FRAC(B)| / |FRAC(A)| [3, 8].")
para("To relate containment to identity, consider the <i>simple mutation model</i> [5, 4]: each base of A is substituted "
     "independently with probability p = 1 − q, where q is the per-base identity (ANI), and spurious matches and repeats "
     "are ignored. A k-mer of A survives in B only if none of its k bases mutated, so its probability of being shared is "
     "q<super>k</super>, and")
eqn(r"\mathbb{E}[C]\approx q^{k}\qquad\Longrightarrow\qquad \hat q = C^{1/k}")
para("Hera, Pierce-Ward and Koslicki [4] derive this point estimate and confidence intervals for FracMinHash. Mash "
     "[2] uses the related transform of the Jaccard index. Two limitations matter here: a global containment mixes "
     "identity and coverage (a containment of 0.5 can mean half the genome is shared at 100% identity, or all of it at "
     "a lower identity), so it cannot give AF; and k-mers spanning an indel are lost just like k-mers spanning a "
     "substitution, so indels lower the estimate.")
h2("2.3 Alignment fraction and positional methods")
para("The <b>alignment fraction</b> AF of A with respect to B is the fraction of A's bases that align to (have a homologous "
     "counterpart in) B. ANI tools report ANI over the aligned part only, together with AF in both directions. FastANI "
     "[7] maps fragments of the query with a minimizer-based mapper; skani [6] chains sparse k-mer seeds (FracMinHash-like "
     "sampling, compression factor 125) and computes identity from chained seeds, with an optional learned regression "
     "correction. Both need positions. AIIPS keeps positions inside the sketch itself.")
h2("2.4 The AIIPS idea")
para("AIIPS samples two k-mer sizes with two independent hash functions, both by the FracMinHash rule (Figure 1). "
     "<b>Anchors</b> are canonical k<sub>1</sub>-mers kept at rate 1/S<sub>1</sub>; <b>probes</b> are canonical "
     "k<sub>2</sub>-mers (k<sub>2</sub> &lt; k<sub>1</sub>) kept at rate 1/S<sub>2</sub>, excluding every k<sub>2</sub>-mer "
     "that overlaps an anchor. The stretch between consecutive anchors is an <b>island</b>; each probe belongs to one "
     "island. Anchors are long and rare, so a matched anchor is strong evidence of homology and a natural seed; probes are "
     "short and dense, so they still match at moderate divergence and measure identity locally. The expected sketch "
     "sizes for a genome of length L are")
eqn(r"|\mathrm{anchors}|\approx \frac{L}{S_1},\qquad |\mathrm{probes}|\approx \frac{L}{S_2}\left(1-\frac{k_1+k_2-1}{S_1}\right)")
para("since a probe start position is excluded when it lies within k<sub>2</sub> − 1 bases before or k<sub>1</sub> − 1 "
     "bases after an anchor start. With the defaults (k<sub>1</sub> = 21, S<sub>1</sub> = 500, k<sub>2</sub> = 15, "
     "S<sub>2</sub> = 100) a 10 Mbp genome has about 20,000 anchors and 93,000 probes (measured: 19,795 and 92,906), "
     "islands average 500 bp, and an island holds about")
eqn(r"m\approx \frac{S_1-(k_1+k_2-1)}{S_2}\approx 4.65\ \mathrm{probes}.")
para("Two consequences drive the design of the comparison:")
bullets([
    "<b>Anchors rarely match at moderate divergence.</b> An anchor is shared with probability q<super>21</super>: 0.81 "
    "at 99% ANI, 0.34 at 95%, 0.11 at 90%, 0.03 at 85%. At 90% ANI matched anchors are on average "
    "S<sub>1</sub>/q<super>21</super> ≈ 4.6 kbp apart, so a walk must cross many unmatched anchors.",
    "<b>Single islands are noisy.</b> Under the mutation model the probability that a homologous island shares no probe "
    "is about (1 − q<super>15</super>)<super>m</super>: 0.06 at 95% ANI, 0.35 at 90%, 0.65 at 85%. A walk that stopped at "
    "the first island without a shared probe would end early and, by keeping only the lucky islands, overestimate ANI. "
    "AIIPS therefore tolerates runs of low islands (Section 3.2).",
])
fig_schem = figure(f_schematic, 6.6,
    "Schematic of an AIIPS comparison. Anchors (blue) split each genome into islands holding probes (orange ticks). The "
    "walk starts at a matched anchor pair (seed) and moves right island by island, each genome keeping its own islands; "
    "shared probes (thin lines) show that paired islands are homologous. Where an anchor of A was mutated away in B (×) "
    "and a mutation created an anchor present only in B (light blue), the islands of the two genomes no longer coincide "
    "and an island of one genome pairs with several of the other's. "
    "The first probe of a contig that precedes its first anchor is promoted to an anchor (green), so every island has an "
    "anchor before it.")

# ================================================================ 3 Methods
h1("3. Methods")
h2("3.1 Sketch")
para("For each record (contig) of a FASTA/FASTQ file, the sketch stores the record length and, in position order:")
bullets([
    "every <b>anchor</b>: a canonical k<sub>1</sub>-mer x with h<sub>1</sub>(x) = fmix64(x ⊕ s<sub>1</sub>) ≤ "
    "⌊(2<super>64</super> − 1)/S<sub>1</sub>⌋, with its start position, its strand (whether the canonical form is the "
    "reverse complement) and its hash;",
    "every <b>probe</b>: a canonical k<sub>2</sub>-mer with h<sub>2</sub>(x) = fmix64(x ⊕ s<sub>2</sub>) ≤ "
    "⌊(2<super>64</super> − 1)/S<sub>2</sub>⌋ that overlaps no anchor, with position, strand and hash, grouped by the "
    "island it lies in;",
    "the <b>promoted anchor</b>: the first probe before the record's first anchor (or the first probe of a record with no "
    "anchor) also becomes an anchor, so the first sampled k-mer of every record is an anchor. It remains a probe of its "
    "island, so containment counts it like the other genome's copy of the same k<sub>2</sub>-mer.",
])
para("fmix64 is the 64-bit finalizer of MurmurHash3 [10], a bijection on 64-bit integers; s<sub>1</sub> and "
     "s<sub>2</sub> are the first two outputs of the SplitMix64 generator [9] seeded with a user seed (default 42), so the "
     "two samplings use independent hash functions. Island lengths follow from anchor positions and need not be stored.")
h2("3.2 Comparison: seeds and walks")
h3("Seeds")
para("Seeds are anchor pairs with equal hash: k<sub>1</sub> anchors of A with k<sub>1</sub> anchors of B, and a promoted "
     "anchor of either genome with any sampled k<sub>2</sub>-mer of the other. Hashes occurring more than 32 times in either "
     "genome (highly repetitive k-mers) are skipped. The strands of the two copies give the relative orientation.")
h3("Walks")
para("From a seed, the walk compares the island after the anchor in A with the island after (same orientation) or before "
     "(opposite orientation) the matched k-mer in B, each genome keeping its own islands. At each step:")
bullets([
    "<b>Pairing.</b> The probes of the current island of A are compared with those of the current island of B; shared "
    "probes are counted.",
    "<b>Crossing.</b> If the next anchors of both genomes are the same k-mer, both are crossed and the walk re-anchors "
    "there. Otherwise the genome whose island ends first crosses its next anchor, so an island of one genome can pair "
    "with several islands of the other (Figure 1). Island ends are compared in A's coordinates after correcting for the "
    "offset between the genomes measured at the last shared probe, which indels change between re-anchorings.",
    "<b>Judging.</b> An island is judged when the walk leaves it, from all the probes it shared with the islands it was "
    "paired with: it is <i>high</i> when the shared fraction of its probes is at least X<super>k<sub>2</sub></super> "
    "(X = 0.8 by default, so in practice one shared probe suffices). When both genomes leave an island at the same step, "
    "either island being high suffices (\"containment in either genome\"); matched next anchors always count as high. "
    "Islands without probes are neutral.",
    "<b>Stopping.</b> The walk continues through low islands and stops after more than G consecutive low judgements "
    "(G = 50 by default) or at a record end, where the other genome's island is clipped. It is then cut back to its last "
    "high judgement and trimmed to its last shared k-mer, except that a walk that reached a record end with no low island "
    "since its last high one keeps its islands up to that end (the record end, not divergence, stopped it). Low islands "
    "inside a walk are kept: dropping them would select for conserved islands and bias ANI upward.",
])
para("A seed walks forward from its anchor in each genome whose anchor it is: one walk when the strands agree, two when "
     "they differ. A seed is skipped when an earlier walk with the same orientation crossed both of its anchors, since it "
     "lies on that walk. Overlapping walks (from repeats) are resolved best first, by score (shared probes plus matched "
     "anchors): an island's probe counts come from the first walk claiming it, and its aligned bases from the union of "
     "the spans claimed by all walks.")
h2("3.3 Estimators")
h3("ANI: joint maximum likelihood over anchors and probes")
para("Inside the claimed walks let N<sub>1</sub> be the number of k<sub>1</sub> anchors crossed and n<sub>1</sub> the "
     "number that matched (same k-mer at the corresponding boundary of the other genome), and N<sub>2</sub> the number of "
     "probes in aligned islands and n<sub>2</sub> the number shared; counts from both genomes are pooled. Seed anchors are "
     "excluded because they were selected for matching. Assuming, as in Section 2.2, that each k-mer of length k matches "
     "independently with probability q<super>k</super>, the log-likelihood of the per-base identity q is")
e_ll = eqn(r"\log L(q)=\sum_{t\in\{1,2\}}\left[\,n_t\,k_t\ln q+(N_t-n_t)\ln\left(1-q^{k_t}\right)\right]")
para("Writing u = ln q, each term n<sub>t</sub>k<sub>t</sub>u + (N<sub>t</sub> − n<sub>t</sub>) ln(1 − e<super>k<sub>t</sub>u</super>) "
     "is concave in u, so log L has a unique maximum, the root of")
eqn(r"\sum_{t}k_t\left[\,n_t-(N_t-n_t)\,\frac{q^{k_t}}{1-q^{k_t}}\right]=0")
para("whose left side decreases from Σ k<sub>t</sub>n<sub>t</sub> &gt; 0 at q → 0 to −∞ at q → 1; it is found by "
     "bisection. With a single k-mer size the root is the containment estimate q̂ = (n/N)<super>1/k</super> of Section 2.2, "
     "which the implementation is tested against. If every k-mer matched, q̂ = 1; if none did, q̂ = 0. Anchors and probes "
     "never overlap, so they share no bases; overlapping probes (rare at density 1/100) are treated as independent, which "
     "makes the likelihood's curvature, and any interval derived from it, overconfident. AIIPS reports the pooled estimate "
     "as <i>ani</i> and the estimates from each genome's counts alone as <i>ani_a</i> and <i>ani_b</i>.")
h3("Alignment fraction")
eqn(r"\mathrm{AF}_A=\frac{\mathrm{bases\ of\ A\ inside\ aligned\ islands}}{\mathrm{all\ bases\ of\ A}},\qquad \mathrm{AF}_B\ \mathrm{likewise}.")
h2("3.4 Optional corrections")
para("Two flags change the comparison, not the sketch; they are off by default so that the default follows the original "
     "specification, and Section 6 reports both settings.")
bullets([
    "<font face='Mono'>--symmetric-probes</font>. Probes are excluded where they overlap an anchor of their own genome. "
    "Where an anchor of B was destroyed by a mutation in A (or created by one), the two genomes apply the exclusion in "
    "different places: the genome lacking the anchor keeps probes the other lacks, and those probes can never be shared. "
    "Approximating, the fraction of one genome's probes affected is about (1 − q<super>k<sub>1</sub></super>)"
    "(k<sub>1</sub> + k<sub>2</sub> − 1)/S<sub>1</sub>, 6% at 90% ANI with the defaults (derived here, not taken from the "
    "literature), which lowers the probe containment by that factor and the ANI estimate by roughly 0.4 points before "
    "pooling with anchors. The flag counts a probe only if it also avoids the other genome's anchors at the same walk "
    "offset.",
    "<font face='Mono'>--bidirectional</font>. Seeds also walk backward, recovering the part of each homologous stretch "
    "before its first matched anchor (contig starts, breakpoints of rearrangements). Second-direction walks start only "
    "from seeds that started a walk.",
])
h2("3.5 Decisions beyond the original specification")
para("The algorithm was specified as: sample anchors and probes, walk from matched anchors through islands while the "
     "containment in either genome is high enough, store island lengths, and combine island results into ANI and AF. The "
     "following choices were made during implementation and are flagged as such: islands are judged when completed "
     "rather than per pair of islands (pairs overlapping by a few bases otherwise fail); the tolerance G for low islands; "
     "trimming to the last shared k-mer; the offset correction for indels; the seed exclusion from the likelihood; the "
     "repeat cap of 32; the promoted anchor also counting as a probe; and the two optional corrections.")
t_defaults = table([
    ["Parameter", "Option", "Default", "Role"],
    ["anchor length k<sub>1</sub>", "--anchor-k", "21", "seeds; specific k-mers"],
    ["anchor scale S<sub>1</sub>", "--anchor-scaled", "500", "1 in 500 k<sub>1</sub>-mers; islands ≈ 500 bp"],
    ["probe length k<sub>2</sub>", "--probe-k", "15", "local identity; must be below k<sub>1</sub>"],
    ["probe scale S<sub>2</sub>", "--probe-scaled", "100", "≈ 4.65 probes per island"],
    ["hash seed", "--seed", "42", "seeds both hash functions"],
    ["low-island tolerance G", "--max-gap", "50", "consecutive low judgements before a walk stops"],
    ["island threshold X", "--min-island-ani", "0.8", "high when shared fraction ≥ X<super>k<sub>2</sub></super>"],
    ["symmetric probes", "--symmetric-probes", "off", "Section 3.4"],
    ["bidirectional walks", "--bidirectional", "off", "Section 3.4"],
], [1.55, 1.35, 0.6, 2.9], "Parameters of AIIPS. The first five are recorded in the sketch; the others are options of "
   "<font face='Mono'>compare</font>.")

# ================================================================ 4 Implementation
h1("4. Implementation")
h2("4.1 Overview")
para("AIIPS is part of the <font face='Mono'>aiips</font> command-line tool (Rust 1.99, edition 2021), which also builds "
     "MinHash and FracMinHash sketches. Its only external dependency is needletail 0.7.4 for FASTA/FASTQ parsing, with "
     "gzip support. The AIIPS code is in <font face='Mono'>src/sketch/aiips.rs</font> (about 2,300 lines) with unit tests "
     "in <font face='Mono'>src/sketch/aiips/tests.rs</font> and command-line tests in "
     "<font face='Mono'>tests/aiips_cli.rs</font>. Usage:")
code("aiips sketch genome_a.fa -a aiips -o genome_a          # writes genome_a.aiips\n"
     "aiips compare genome_a.aiips genome_b.fa              # sketches B in memory\n"
     "aiips compare a.fa b.fa -a aiips --symmetric-probes --bidirectional")
para("The output of <font face='Mono'>compare</font> is a tab-separated table: input statistics, then "
     "<i>ani, af_a_to_b, af_b_to_a, ani_a, ani_b</i> and the counts n<sub>1</sub>, N<sub>1</sub>, n<sub>2</sub>, "
     "N<sub>2</sub> behind the estimate. A sketch file is text: a header line with every parameter and the number of "
     "records, anchors and probes; a summary line; then per record a line <font face='Mono'>c LENGTH</font> followed by "
     "its anchors (<font face='Mono'>a</font>, or <font face='Mono'>p</font> for a promoted anchor), each followed by the "
     "probes (<font face='Mono'>q</font>) of its island, every k-mer as position, strand and hash. Reading is strict "
     "(order, bounds, thresholds, overlap rules and counts are all checked).")
h2("4.2 Sketching engine")
para("Sketching dominates the run time: two 64-bit hashes per base of each genome. The engine is organised so that the "
     "inner loop does as little as possible:")
bullets([
    "<b>Streaming decode.</b> needletail delivers each record as raw bytes including line breaks. AIIPS decodes it in "
    "chunks of 65,536 bases, dropping line breaks on the fly: lines are cut where the first line's width predicts, the "
    "number of breaks consumed is checked against the record's count (computed in one vectorised pass), and the record "
    "is decoded again by search if they differ. Eight bases are converted to 2-bit codes per 64-bit word operation "
    "(validity by exact zero-byte tests; for A/C/G/T in either case, bits 1–2 of the ASCII code are 0, 1, 3, 2, mapped to "
    "0–3 by one Gray-code step). No copy of the record is made and no buffer grows with record length.",
    "<b>Eight lanes.</b> A chunk is cut into 8 segments scanned side by side, so the rolling forward and reverse codes, "
    "both canonical minima, both hashes and both threshold tests of 8 positions are one vector operation each. Each lane "
    "first warms up on the k<sub>1</sub> − 1 bases before its segment (carried over between chunks), so it samples exactly "
    "the k-mers a sequential scan would. Lanes load 8 bases at a time and an 8×8 byte transpose (word arithmetic) turns "
    "them into one code per lane per step.",
    "<b>Rare hits out of the loop.</b> About 9 steps in 10 sample nothing (rates 1/100 and 1/500 over 8 lanes) and cost "
    "only the threshold test. The others are written to a 256-entry buffer with their hashes and strands and emitted in "
    "bulk, so the loop's state stays in vector registers.",
    "<b>Runtime dispatch.</b> The same code is compiled for AVX-512, AVX2 and baseline x86-64 and chosen at run time; "
    "other architectures use the baseline. Per-record buffers are kept by the thread and reused by the next sketch.",
])
h2("4.3 Comparison engine")
para("Seeds are found by joining the sorted k<sub>1</sub> anchor hashes of both genomes; promoted anchors (one per contig) "
     "are matched against the other genome's probes by binary search rather than by sorting all probes. Seeds are walked "
     "upstream first, so later seeds are recognised as already covered. Each walk touches an island a bounded number of "
     "times, so the comparison is close to linear in the sketch sizes; it takes about 7–8 ms for two 10 Mbp genomes.")
h2("4.4 Verification")
para("All 204 tests of the repository pass (54 library unit tests, 21 unit tests of the binary, 129 integration tests), "
     "as do Clippy with warnings denied and a smoke script. The AIIPS-specific tests check:")
bullets([
    "<b>Sketch construction</b> against an independent brute-force builder (every window hashed from scratch, exclusion "
    "and promotion applied literally), on records with ambiguous bases, lane and chunk boundaries, short and empty "
    "records and line breaks (LF, CRLF, ragged and blank lines), for every instruction set the machine offers; the "
    "word encoder against the lookup table for every byte value; the 8×8 transpose against a naive transpose; decoding "
    "against plain line-break removal, including a deliberately mispredicted line width that must be caught.",
    "<b>The likelihood</b> against the closed form (n/N)<super>1/k</super> for one k-mer size and against a fine grid "
    "search for two.",
    "<b>Comparisons</b> on exact cases: identical genomes (ANI 1, AF ≥ 0.995), the reverse complement cut into contigs "
    "(ANI 1, AF ≥ 0.99), half of a genome shared with a breakpoint (AF 0.5 and 0.75 within 0.005), unrelated genomes "
    "(AF below 0.005), and 1 Mbp genomes substituted at 1% and 5% (ANI within 0.002 and 0.006).",
    "<b>Command line</b>: sketch format round trip and strict parsing, option validation, parameter-mismatch errors, and "
    "identical results from sketches and from sequence files.",
])
h2("4.5 Performance engineering")
para("The previous commit took 118 ms per 10 Mbp pair against 95 ms for skani. Measuring on this shared server required "
     "care: wall-clock timings inside a process varied up to threefold with the machine's load and could not be trusted, "
     "so per-phase costs were measured as CPU time (getrusage) over repeated runs, and tools were compared interleaved on "
     "one pinned, idle core. The changes, in order of effect:")
t_perf = table([
    ["Bottleneck found", "Change", "Measured effect"],
    ["Indels shifted the two genomes' walk coordinates, so walks missed the anchor pairs they should re-anchor on; on a "
     "10 Mbp pair with indels at 10% divergence about 47 walks each re-covered most of the genome (918,676 steps)",
     "offset correction from the last shared probe; skip seeds whose two anchors an earlier walk crossed",
     "1 walk; comparison 124 → 7 ms (389 → 11 ms with both flags)"],
    ["needletail's sequence accessor copies each record to strip line breaks", "decode raw bytes in chunks while dropping breaks",
     "parsing ≈ 3 ms per 10 Mbp genome; no per-record copy"],
    ["scalar transposition of lane inputs (≈ 40% of the scan); lane state kept in memory because the rare emit path "
     "indexed it; register spills on every sampled step", "word loads + 8×8 transpose; hit buffer emitted in bulk; "
     "counters moved out of the loop", "scan ≈ 21 ms of CPU per 10 Mbp genome"],
    ["seed generation sorted every probe of both genomes", "binary search for promoted anchors", "10 → 2.4 ms"],
    ["record-sized code buffer (10 MB) and page faults", "64 Ki-base chunks; per-thread buffer reuse",
     "peak memory 44.7 → 26.8 MB; page faults 15.1k → 10.7k (skani 17.2k)"],
], [2.55, 2.1, 1.75], "Performance work between the previous commit (2a05b3f) and this one (1ee7ba0), measured on "
   "10 Mbp genomes. Sketches are identical before and after; the comparison changes alter a few results slightly "
   "(Section 6.3).", align_center_from=3)
para("Changes that did not help were reverted: forcing the hash function to inline and compiling as one code-generation "
     "unit (both slowed the scan), and sizing needletail's buffer to the file size (more memory on draft assemblies). "
     "The remaining cost is dominated by the two fmix64 hashes per base, whose 64-bit vector multiplies compete for one "
     "execution port; skani computes one hash per base.")

# ================================================================ 5 Setup
h1("5. Experimental setup")
h2("5.1 Simulated genome pairs")
para("<font face='Mono'>scripts/simulate.py</font> (Python standard library, fixed seeds) generates a uniform random "
     "genome A of 10 Mbp and 32 derived genomes B in four scenario families (Table 3). Substitutions replace a base by "
     "one of the three others; indels occur at a tenth of the substitution rate, as insertions of random bases or "
     "deletions, with geometric lengths of mean 3 capped at 30. A second base genome carries 40 repeat families (2–5 "
     "copies of 1–7 kbp, each copy up to 1% diverged, half reverse-complemented) for the rearrangement scenarios.")
t_scen = table([
    ["Family", "Pairs", "Construction of B", "True ANI", "True AF"],
    ["Substitutions", "8", "A with uniform substitutions at rates 0.1–20%", "80.0–99.9%", "1"],
    ["Partial overlap", "9", "20, 50 or 90% of A in 20–100 kbp blocks (in order), substituted at 1, 5 or 10%, with foreign "
     "random sequence (half the unshared length) between blocks", "≈ 90, 95, 99%", "A: 0.2–0.9; B: 0.34–0.95"],
    ["Indels + draft contigs", "12", "substitutions at 1, 5 or 10% plus indels; B as one record or cut into contigs of mean "
     "50, 10 or 3 kbp, each in random orientation", "≈ 90, 95, 99%", "0.985–0.999"],
    ["Repeats + rearrangements", "3", "A with repeat families; B substituted at 1, 5 or 10%, cut into about 40 blocks, "
     "shuffled and partly reverse-complemented", "≈ 90, 95, 99%", "1"],
], [1.35, 0.45, 3.0, 0.85, 0.85], "Simulated scenarios (10 Mbp genomes).", align_center_from=1)
h2("5.2 Ground truth")
para("Truth is computed from the simulation itself. <b>ANI</b> is one minus the number of substitutions divided by the "
     "number of aligned columns, i.e. homologous base pairs, so indel gaps are excluded. <b>AF</b> of A is the number of "
     "A's bases with a homologous base in B divided by A's length, and likewise for B; deleted bases and foreign or "
     "inserted sequence have no homolog. Because truth excludes gaps while every k-mer method loses k-mers that span an "
     "indel, indel scenarios carry an expected downward bias for both tools.")
h2("5.3 Protocol")
para("<font face='Mono'>scripts/benchmark.py</font> runs each tool on the two FASTA files of every pair, so sketching is "
     "included: <font face='Mono'>aiips compare A B -a aiips</font> (with and without the two flags), the previous commit's "
     "binary with default options, and <font face='Mono'>skani dist -t 1 --min-af 0 -q B -r A</font> (skani 0.4.0 with "
     "default parameters, including its learned-ANI regression). All runs use one thread, pinned to one core of a shared "
     "two-socket Intel Xeon Gold 5416S server; each tool and pair is run three times, interleaved with the other tools, "
     "and the fastest run is kept. Wall time is measured around the process, user and system CPU time come from getrusage "
     "and peak resident memory from <font face='Mono'>/usr/bin/time</font>. Errors are reported in percentage points; AF "
     "errors average both directions.")

# ================================================================ 6 Results
h1("6. Results")
h2("6.1 Overall")
rows = [["Tool", "ANI MAE", "ANI bias", "max |ANI err|", "AF MAE", "Wall (ms)", "User (ms)", "System (ms)", "Memory (MiB)"]]
for t, lab in TOOLS:
    x = ALL[t]
    rows.append([lab + (" ¹" if x["missing"] else ""), fmt(x["ani"]), sgn(x['bias']), fmt(x["ani_max"]), fmt(x["af"], 1),
                 f"{x['wall_s']:.0f}", f"{x['user_s']:.0f}", f"{x['sys_s']:.0f}", f"{x['rss']:.1f}"])
t_overall = table(rows, [1.6, 0.6, 0.62, 0.75, 0.6, 0.62, 0.62, 0.7, 0.75],
    "All 32 pairs of 10 Mbp genomes: mean absolute error (MAE), signed bias and maximum absolute error of ANI, MAE of AF "
    "(percentage points), and mean cost per pair. ¹ skani reported no value for the 80% ANI pair (presumably its default "
    "screen at about 80% identity); its error statistics cover 31 pairs.")
para(f"With default settings AIIPS halves skani's ANI error ({a_['ani']:.2f} vs {s_['ani']:.2f} points) and reduces its AF "
     f"error from {s_['af']:.1f} to {a_['af']:.1f} points; the two corrections reduce them further to {b_['ani']:.2f} and "
     f"{b_['af']:.1f}. Both tools are biased downward by default, AIIPS less so; with the corrections the remaining bias "
     f"({b_['bias']:+.2f}) comes from the indel scenarios. AIIPS is faster than skani per pair "
     f"({a_['wall_s']:.0f} vs {s_['wall_s']:.0f} ms, equal user time, less system time) and uses about a third less memory. "
     "The previous commit gives the same accuracy as this one by default (identical sketches; Section 6.3) at higher cost.")
figure(f_cost, 6.6, "Wall time and peak memory per pair by scenario family (10 Mbp genomes, one thread).")
h2("6.2 Accuracy by scenario")
rows = [["Family", "Condition", "n", "skani", "aiips", "aiips, flags", "skani", "aiips", "aiips, flags"]]
for g in GROUP_ORDER:
    cells = [g[0], g[1], str(len(GROUPS[g]))]
    stats = {t: gsum(t, g) for t in ("skani", "aiips", "aiips:both")}
    for t in ("skani", "aiips", "aiips:both"):
        x = stats[t]
        cells.append(f"{x['ani']:.2f} ({sgn(x['bias'])})" + (" ¹" if x["missing"] else ""))
    for t in ("skani", "aiips", "aiips:both"):
        cells.append(fmt(stats[t]["af"], 1))
    rows.append(cells)
t_groups = table([["", "", "", "ANI: MAE (bias)", "", "", "AF: MAE", "", ""]] + rows,
    [1.05, 1.05, 0.25, 0.85, 0.8, 0.85, 0.5, 0.5, 0.65],
    "Accuracy by scenario and condition (percentage points; 10 Mbp genomes). Partial-overlap and draft conditions each "
    "pool the three divergences (≈ 99, 95, 90% ANI). ¹ one pair (80% ANI) not reported by skani.", header_rows=2,
    align_center_from=2, spans=[(3, 0, 5, 0), (6, 0, 8, 0), (0, 0, 2, 0)])
para("<b>Substitutions.</b> All tools are accurate above 99% ANI. Between 90 and 97% skani underestimates ANI by 0.36 "
     "points on average and loses 5 points of AF; AIIPS has no AF error and, with the symmetric-probe correction, an ANI "
     "error of 0.06. At 80–85% ANI skani's AF drops to about 0.4 (85%) or it reports nothing (80%), whereas AIIPS still "
     "aligns 90–100% of the genome and estimates ANI within 0.2 points (Figure 4). Figure 3 summarises Table 5.")
para("<b>Partial overlap.</b> skani's ANI error is largest here (0.6–1.3 points); AIIPS stays below 0.3 by default and "
     "below 0.15 with the corrections, and the bidirectional walks bring the AF error below 1 point.")
para("<b>Indels and draft assemblies.</b> Both tools are biased downward as expected (Section 5.2); AIIPS less so "
     "(condition means −0.37 to −0.68 vs −0.75 to −0.87 points). AF is accurate for contigs of 50 kbp or more, but both tools lose a large "
     "part of the genome on 3 kbp contigs at 90% ANI (AF 0.41–0.50 for AIIPS, 0.50 for skani): short contigs at low "
     "identity often contain no matched anchor at all. skani is slightly better there (AF error 22.1 vs 23.1 with the "
     "corrections, 28.3 without).")
para("<b>Repeats and rearrangements.</b> AIIPS is unaffected by shuffled and inverted blocks (AF error ≤ 0.6) because "
     "walks start from any matched anchor in either orientation; skani loses about 5 points of AF.")
figure(f_groups, 5.3, f"Mean absolute ANI and AF error by condition (from Table {t_groups}).")
figure(f_sweep, 6.6, "Substitution sweep (true AF = 1 throughout). Left: ANI estimate minus truth; right: estimated AF "
       "of A. skani reports no value at 80% ANI. The previous commit's curve coincides with this commit's default.")
h2("6.3 Changes since the previous commit")
para("The sketches of the two commits are identical, so differences come only from the walk changes of Section 4.5. "
     "Over the 32 pairs the default ANI and AF errors are unchanged (0.34 and 5.0 points). Per scenario, the indel "
     "scenarios improve slightly (e.g. −1.23 → −1.21 points at 90% ANI, one record), while the 90%-shared partial-overlap "
     "pairs worsen by up to 0.08 points of ANI (default: 0.23 → 0.26 MAE); these are side effects of the different walk "
     "order and are well within the differences between tools.")
h2("6.4 The likelihood")
figure(f_like, 6.6, f"Log-likelihood (Equation {e_ll}) for the 10 Mbp pair at true ANI 0.9499 (5% substitutions), "
       f"with the pooled counts n<sub>1</sub>/N<sub>1</sub> = 13,410/39,640 anchors and n<sub>2</sub>/N<sub>2</sub> = "
       f"83,194/186,093 probes. Left: anchors only, probes only and joint, each shifted to its maximum; right: the joint "
       f"curve near its maximum q̂ = {mle:.4f}. The shaded likelihood-ratio interval [{ci[0]:.4f}, {ci[1]:.4f}] "
       "excludes the truth, which illustrates both the overconfidence of the independence assumption and the downward "
       "bias of the default estimate that the symmetric-probe correction removes (0.9505 with the correction).")
para("Probes dominate the joint likelihood (they are about five times as numerous); anchors, being longer, are more "
     "sensitive per k-mer. Their separate maxima agree to within 0.001, which supports pooling them.")
h2("6.5 Choice of the low-island tolerance G")
gap_rows = [["G", "ANI MAE", "ANI bias", "AF MAE (all)", "AF MAE (substitutions)"]]
for label, tool in [("0 (stop at first low island)", "aiips:g0"), ("2", "aiips"), ("5", "aiips:g5"), ("10", "aiips:g10"),
                    ("20", "aiips:g20"), ("50 (default)", "aiips:g50")]:
    x = summary([r for r in GAP if r["tool"] == tool])
    xs = summary([r for r in GAP if r["tool"] == tool and r["category"] == "sub"])
    gap_rows.append([label, fmt(x["ani"]), sgn(x['bias']), fmt(x["af"], 1), fmt(xs["af"], 1)])
table(gap_rows, [1.9, 0.9, 0.9, 1.1, 1.5], "Effect of G during development (32 pairs of 5 Mbp genomes, an intermediate "
      "version without trimming and with per-island claims; absolute values differ from the final version, the trend is "
      "what motivated the default).", align_center_from=1)
para("Stopping at the first low island (G = 0) loses most of the alignment below 95% ANI and inflates ANI at low "
     "identity (selection of conserved islands). Larger G only extends true alignments: because probes are compared "
     "only between paired islands, a spurious shared 15-mer between two random islands of about five probes each has "
     "probability of order 25/4<super>15</super> ≈ 2 × 10<super>−8</super>, and trailing low islands are trimmed.")
h2("6.6 Earlier results on 5 Mbp genomes")
rows = [["Variant (previous commit)", "ANI MAE", "ANI bias", "AF MAE", "Wall (ms)", "Memory (MiB)"]]
for tool, lab in [("aiips", "aiips"), ("aiips:sym", "<font face='Mono'>--symmetric-probes</font>"), ("aiips:bidir", "<font face='Mono'>--bidirectional</font>"),
                  ("aiips:sym+bidir", "both flags"), ("skani", "skani 0.4.0 ¹")]:
    x = summary([r for r in R5 if r["tool"] == tool])
    rows.append([lab, fmt(x["ani"]), sgn(x['bias']), fmt(x["af"], 1), f"{x['wall_s']:.0f}", f"{x['rss']:.1f}"])
table(rows, [2.0, 0.8, 0.8, 0.8, 0.9, 1.0], "The same benchmark on 5 Mbp genomes with the previous commit "
      "(benchmarks/simulated-5mb.tsv; wall time at 10 ms resolution, unpinned). It isolates the two corrections: "
      "<font face='Mono'>--symmetric-probes</font> reduces ANI error, <font face='Mono'>--bidirectional</font> reduces "
      "AF error, and they combine. ¹ one pair not reported.",
      align_center_from=1)


h2("6.7 k-mer sizes")
para("To see how the k-mer lengths shape the trade-off, the 10 Mbp benchmark was repeated with six other pairs "
     "(k<sub>1</sub>, k<sub>2</sub>), keeping the default sampling rates (1/500 anchors, 1/100 probes) and the default "
     "comparison options (no flags); skani was run again in the same session as a reference. Because the sampling rates "
     "fix how many k-mers are hashed and stored, time and memory hardly change with k; accuracy does (Tables "
     f"{counters['tab'] + 1} and {counters['tab'] + 2}).")
KTOOLS = [("aiips", "21/15 (default)"), ("aiips:k17_13", "17/13"), ("aiips:k19_13", "19/13"), ("aiips:k21_13", "21/13"),
          ("aiips:k21_17", "21/17"), ("aiips:k25_17", "25/17"), ("aiips:k31_21", "31/21"), ("skani", "skani 0.4.0")]
rows = [["k<sub>1</sub>/k<sub>2</sub>", "ANI MAE", "ANI bias", "AF MAE", "Wall (ms)", "User (ms)", "Memory (MiB)"]]
for t, lab in KTOOLS:
    x = summary([r for r in KSWEEP if r["tool"] == t])
    rows.append([lab + (" ¹" if x["missing"] else ""), fmt(x["ani"]), sgn(x["bias"]), fmt(x["af"], 1),
                 f"{x['wall_s']:.0f}", f"{x['user_s']:.0f}", f"{x['rss']:.1f}"])
table(rows, [1.5, 0.75, 0.75, 0.75, 0.8, 0.8, 0.95], "k-mer sizes: all 32 pairs of 10 Mbp genomes (percentage "
      "points; mean cost per pair, one pinned core, fastest of 3). Sampling rates and options at their defaults. "
      "¹ skani reported no value for the 80% ANI pair.", align_center_from=1)
fams = [("sub", "Substitutions"), ("partial", "Partial overlap"), ("draft", "Indels + drafts"), ("rearr", "Repeats + rearr.")]
rows = [["k<sub>1</sub>/k<sub>2</sub>"] + [f for _, f in fams]]
for t, lab in KTOOLS:
    cells = [lab]
    for c, _ in fams:
        x = summary([r for r in KSWEEP if r["tool"] == t and r["category"] == c])
        cells.append(f"{x['ani']:.2f} / {x['af']:.1f}")
    rows.append(cells)
table(rows, [1.5, 1.25, 1.25, 1.25, 1.25], "k-mer sizes by scenario family: ANI MAE / AF MAE (percentage points).",
      align_center_from=1)
low = sorted({r["scenario"] for r in KSWEEP if r["category"] == "sub" and float(r["true_ani"]) < 0.95},
             key=lambda s: -float(next(r for r in KSWEEP if r["scenario"] == s)["true_ani"]))
sel = [("aiips", "21/15"), ("aiips:k17_13", "17/13"), ("aiips:k21_17", "21/17"), ("aiips:k25_17", "25/17"),
       ("aiips:k31_21", "31/21"), ("skani", "skani")]
rows = [["true ANI"] + [l for _, l in sel]]
for sc in low:
    cells = []
    for t, _ in sel:
        r = next(r for r in KSWEEP if r["scenario"] == sc and r["tool"] == t)
        cells.append("no result" if r["ani"] == "NA" else
                     f"{sgn(100 * (float(r['ani']) - float(r['true_ani'])))} / {float(r['af_a']):.2f}")
    rows.append([f"{float(r['true_ani']):.3f}"] + cells)
table(rows, [0.75, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95], "Low-identity substitution pairs: ANI error (points) / AF of "
      "A (true AF = 1).", align_center_from=1)
para("<b>Shorter k-mers</b> (17/13) survive mutation more often, so walks keep extending at low identity: they give the "
     "lowest AF error overall (3.4 vs 5.0 points), most of all on fragmented drafts, and still align the whole genome at "
     "80% ANI. <b>Longer probes</b> (21/17) give the lowest ANI error but lose AF below 85% ANI. <b>Long k-mers</b> "
     "(25/17, 31/21) leave almost no matching anchors below about 85% ANI, and the alignment collapses (31/21 aligns 3% of "
     "the genome at 80%). The simulated genomes are uniform random sequence, which hides the main cost of short k-mers: "
     "in real genomes 13- and 17-mers recur far more often, producing spurious seeds and probe matches, so a change of "
     "default needs validation on real genome pairs first.")

# ================================================================ 7 Discussion
h1("7. Limitations and next steps")
h2("7.1 Limitations")
bullets([
    "<b>Simulated data only.</b> Genomes are uniform random sequence without the compositional bias, low-complexity "
    "regions and gene-level structure of real genomes; the indel and repeat models are simple. skani's learned ANI "
    "regression was trained on real assemblies and may be disadvantaged here.",
    "<b>Indels.</b> Like other k-mer methods, AIIPS counts a k-mer spanning an indel as a mismatch, so ANI over gap-free "
    "aligned columns is underestimated when indels are frequent.",
    "<b>Fragmented assemblies at low identity.</b> Contigs without any matched anchor are missed, which limits AF on 3 kbp "
    "contigs at 90% ANI (shared with skani).",
    "<b>Uncertainty.</b> The likelihood assumes independent k-mer matches; its curvature overstates precision, so no "
    "confidence interval is reported yet.",
    "<b>Measurements.</b> Most scenario conditions average three pairs; timings on the shared server vary by about ±10% "
    "between runs. skani was run with default parameters only.",
    "<b>Cost floor.</b> Two 64-bit hashes per base dominate run time; further single-thread speed-ups would require a "
    "cheaper hash, which would change the sketch definition.",
])
h2("7.2 Next steps")
bullets([
    "Decide whether <font face='Mono'>--symmetric-probes</font> and <font face='Mono'>--bidirectional</font> become the "
    "defaults; they lower both errors at about 5% extra time.",
    "Validate on real genome pairs with alignment-based reference values (e.g. ANIm or dnadiff identities and aligned "
    "fractions), across a range of taxa, GC content and assembly quality.",
    "Compare against skani with <font face='Mono'>--no-learned-ani</font> and its <font face='Mono'>--medium</font> / "
    "<font face='Mono'>--slow</font> presets (recommended for fragmented assemblies), and against FastANI.",
    "Calibrate an uncertainty estimate (e.g. a block bootstrap over islands) instead of the likelihood curvature.",
    "Validate shorter k-mers (17/13, Section 6.7) on real genomes, where repeated short k-mers may cost accuracy.",
    "Many-versus-many comparisons with an index over anchors, and optional multi-threading.",
])

# ================================================================ References
h1("References")
refs = [
    "Broder AZ. On the resemblance and containment of documents. <i>Proceedings of Compression and Complexity of "
    "Sequences</i>, 1997, 21–29.",
    "Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. Mash: fast genome and metagenome "
    "distance estimation using MinHash. <i>Genome Biology</i> 17:132, 2016.",
    "Irber L, Brooks PT, Reiter T, Pierce-Ward NT, Hera MR, Koslicki D, Brown CT. Lightweight compositional analysis of "
    "metagenomes with FracMinHash and minimum metagenome covers. <i>bioRxiv</i>, 2022.",
    "Hera MR, Pierce-Ward NT, Koslicki D. Deriving confidence intervals for mutation rates across a wide range of "
    "evolutionary distances using FracMinHash. <i>Genome Research</i> 33:1061–1068, 2023.",
    "Blanca A, Harris RS, Koslicki D, Medvedev P. The statistics of k-mers from a sequence undergoing a simple mutation "
    "process without spurious matches. <i>Journal of Computational Biology</i> 29(2):155–168, 2022.",
    "Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. <i>Nature "
    "Methods</i> 20:1661–1665, 2023.",
    "Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic "
    "genomes reveals clear species boundaries. <i>Nature Communications</i> 9:5114, 2018.",
    "Koslicki D, Zabeti H. Improving MinHash via the containment index with applications to metagenomic analysis. "
    "<i>Applied Mathematics and Computation</i> 354:206–215, 2019.",
    "Steele GL, Lea D, Flood CH. Fast splittable pseudorandom number generators. <i>Proceedings of OOPSLA</i>, 2014.",
    "Appleby A. MurmurHash3 (64-bit finalizer fmix64). https://github.com/aappleby/smhasher, 2011.",
]
for i, r in enumerate(refs, 1):
    story.append(Paragraph(f"[{i}]&nbsp;&nbsp;{r}", styles["ref"]))

# ================================================================ Appendix
h1("Appendix: reproducing the results")
code("# build (Rust toolchain from environment.yml)\n"
     "cargo build --locked --release\n"
     "cargo test --locked --release\n\n"
     "# simulate 32 pairs of 10 Mbp genomes with truth (about 1 minute)\n"
     "python3 scripts/simulate.py sim10m --length 10000000\n\n"
     "# benchmark: aiips (default and both flags), an earlier build, skani\n"
     "python3 scripts/benchmark.py sim10m/truth.tsv results.tsv \\\n"
     "    --aiips target/release/aiips --skani \"$(command -v skani)\" \\\n"
     "    --repeat 3 --cpu 63 --variant both=\"--symmetric-probes --bidirectional\" \\\n"
     "    --other-aiips committed=/path/to/previous/aiips")
code("# k-mer size sweep (Section 6.7)\n"
     "python3 scripts/benchmark.py sim10m/truth.tsv ksweep.tsv --aiips target/release/aiips \\\n"
     "    --skani \"$(command -v skani)\" --repeat 3 --cpu 63 \\\n"
     "    --variant k17_13=\"--anchor-k 17 --probe-k 13\" --variant k21_17=\"--anchor-k 21 --probe-k 17\"  # etc.\n\n"
     "# this report\n"
     "python3 reports/build_report.py . reports/aiips_progress_report.pdf")
para("All tables and figures are computed from the committed files in <font face='Mono'>benchmarks/</font>: "
     "<font face='Mono'>simulated-10mb.tsv</font> with <font face='Mono'>simulated-10mb-truth.tsv</font>, "
     "<font face='Mono'>simulated-10mb-ksweep.tsv</font>, <font face='Mono'>simulated-5mb.tsv</font> and "
     "<font face='Mono'>simulated-5mb-gap-sweep.tsv</font>.")


# ---------------------------------------------------------------- build
class Doc(SimpleDocTemplate):
    def afterFlowable(self, flowable):
        if isinstance(flowable, Heading):
            level, text = flowable._toc
            self.notify("TOCEntry", (level, text, self.page))


def decorate(canvas, doc):
    canvas.saveState()
    canvas.setFont("Sans", 7.5)
    canvas.setFillColor(colors.HexColor("#666666"))
    if doc.page > 1:
        canvas.drawString(0.8 * inch, 0.5 * inch, "AIIPS progress report — 17 November 2026")
        canvas.drawRightString(letter[0] - 0.8 * inch, 0.5 * inch, f"{doc.page}")
    canvas.restoreState()


doc = Doc(OUT, pagesize=letter, leftMargin=0.85 * inch, rightMargin=0.85 * inch, topMargin=0.8 * inch,
          bottomMargin=0.8 * inch, title="AIIPS progress report", author="aiips project",
          subject="Anchor-Indexed Island-Probe Sketching: theory, implementation and benchmark")
doc.multiBuild(story, onFirstPage=decorate, onLaterPages=decorate)
print("wrote", OUT)
