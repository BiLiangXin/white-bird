"""绘制后续研究计划的甘特图（docs/figures/D1_research_timeline.png）。

用法：
  python -I src/plan/plot_timeline.py --out docs/figures
"""
from __future__ import annotations

import argparse
import sys
from datetime import date, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "audit"))
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from gwfg_independent_audit import DAY, INK, INK2, MUTED, NIGHT  # noqa: E402

START = date(2026, 11, 16)  # 第1周周一

# (工作包, 起始周, 结束周, 所属论文)
TASKS = [
    ("WP0 分析预案与数据授权", 1, 1, 1),
    ("WP1 外部数据提取与质控", 1, 3, 1),
    ("WP2 两驻留期昼夜格局", 3, 4, 1),
    ("WP3 月光天然实验", 4, 6, 1),
    ("WP4 机制判别", 6, 8, 1),
    ("WP5 月夜农田利用", 7, 9, 1),
    ("WP6 稳健性与功效", 8, 10, 1),
    ("WP7 图件与写作", 10, 14, 1),
    ("WP8 内审、修改与投稿", 15, 16, 1),
    ("WP9 取回原始GPS与晨昏轨迹", 1, 6, 2),
    ("WP10 全日轨迹重建出发事件", 10, 14, 2),
    ("WP11 出发风选择模型", 14, 20, 2),
    ("WP12 论文2写作", 20, 26, 2),
]
MILESTONES = [
    (1, "M1 预案定稿"),
    (6, "M2 验证集复现？"),
    (10, "M3 结果冻结"),
    (16, "M4 论文1投稿"),
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    n = len(TASKS)
    for i, (name, s, e, paper) in enumerate(TASKS):
        y = n - i
        col = DAY if paper == 1 else NIGHT
        ax.barh(y, e - s + 1, left=s - 0.5, height=0.62, color=col, edgecolor="white", linewidth=1.5)
        ax.text(0.3, y, name, ha="right", va="center", fontsize=9.5, color=INK)
    for wk, label in MILESTONES:
        ax.axvline(wk + 0.5, color=MUTED, linewidth=1, linestyle=":")
        ax.plot(wk + 0.5, n + 0.95, marker="D", color=INK, markersize=7)
        ax.text(wk + 0.5, n + 1.45, label, ha="center", fontsize=8.5, color=INK2)
    ax.axvspan(11.5, 13.5, color=MUTED, alpha=0.15, linewidth=0)
    ax.text(12.5, 0.15, "春节假期", ha="center", fontsize=8.5, color=INK2)

    weeks = list(range(1, 27, 2))
    ax.set_xticks(weeks, [f"第{w}周\n{(START + timedelta(weeks=w - 1)).strftime('%m/%d')}" for w in weeks],
                  fontsize=8.5)
    ax.set_yticks([])
    ax.set_xlim(0.4, 26.6)
    ax.set_ylim(-0.3, n + 2.0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="y", visible=False)
    handles = [Patch(color=DAY, label="论文1：昼夜活动与月光（主线）"),
               Patch(color=NIGHT, label="论文2：迁飞出发与风选择（视数据而定）")]
    ax.legend(handles=handles, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.86), fontsize=9)
    fig.suptitle("图D1  后续研究计划（2026年11月—2027年5月）", x=0.01, ha="left", color=INK, fontsize=12)
    fig.subplots_adjust(left=0.22, right=0.98, top=0.9, bottom=0.12)
    fig.savefig(args.out / "D1_research_timeline.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
