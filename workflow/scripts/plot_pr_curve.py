# Copyright 2025 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score


def _concat_df(files, cols, label):
    frames = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        sub = df[cols].copy()
        sub["y"] = label
        frames.append(sub)
    if not frames:
        return pd.DataFrame(columns=cols + ["y"])
    return pd.concat(frames, ignore_index=True)


cols = ["Danc", "Dplus", "df", "fd", "DD", "U", "Q"]
stats = {
    "Danc": "$D_{anc}$",
    "Dplus": "$D^+$",
    "df": "$d_f$",
    "fd": "$f_d$",
    "DD": "$D_D$",
    "U": "$U50$",
    "Q": "$Q95$",
}

wo_all = list(snakemake.input.wo)
neu_all = list(snakemake.input.neu)
sel_all = list(snakemake.input.sel)

R = 1000
recall_grid = np.linspace(0, 1, 200)
mean_curves = {col: [] for col in cols}
ap_logs = {col: [] for col in cols}

for r in range(R):
    random.seed(2487 + r)
    wo_files = random.sample(wo_all, 980)
    neu_files = random.sample(neu_all, 18)
    sel_files = random.sample(sel_all, 2)

    df_wo = _concat_df(wo_files, cols=cols, label=0)
    df_neu = _concat_df(neu_files, cols=cols, label=0)
    df_sel = _concat_df(sel_files, cols=cols, label=1)
    df = pd.concat([df_wo, df_neu, df_sel], ignore_index=True)
    y_true = df["y"].to_numpy()

    for col in cols:
        y_score = pd.to_numeric(df[col], errors="coerce").to_numpy()
        mask = np.isfinite(y_score)
        if mask.sum() == 0:
            continue
        precision, recall, _ = precision_recall_curve(y_true[mask], y_score[mask])
        ap = average_precision_score(y_true[mask], y_score[mask])
        ap_logs[col].append(ap)

        # reverse recall from 1->0 to 0->1
        # because xp must be increasing in np.interp()
        prec_interp = np.interp(recall_grid, recall[::-1], precision[::-1])
        mean_curves[col].append(prec_interp)

plt.figure(figsize=(5, 5), dpi=300)

pos_rate = 0.002
plt.hlines(pos_rate, 0, 1, linestyles="--", linewidth=1, label=f"baseline = {pos_rate:.3f}")

for col in cols:
    if not ap_logs[col]:
        continue
    arr_ap = np.array(ap_logs[col])
    mean_ap = arr_ap.mean()
    mean_prec = np.mean(mean_curves[col], axis=0)
    plt.plot(recall_grid, mean_prec, lw=2, label=f"{stats[col]} (mean AP = {mean_ap:.2f})")

plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title(f"$s$ = {snakemake.wildcards.s}")
plt.legend(loc="lower left", fontsize=8)
plt.tight_layout()
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
