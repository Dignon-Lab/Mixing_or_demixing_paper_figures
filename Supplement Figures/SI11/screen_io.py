#!/usr/bin/env python3
"""Shared data loader for Figure 6 / SI Figure 11: parses analysis_37377.out
(make-whole corrected run, 18 partners x LAF-1 V2, 5 compositions, 298 K)."""
import os, re
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_FILE = os.path.join(HERE, "analysis_37377.out")
F1_FILE = os.path.join(HERE, "figure1_data.csv")

def _partner_to_label(p):
    if p == "LAF_1_IDR_WT": return "WT"
    if p == "LAF_1_V1": return "V1"
    if p == "LAF_1_V3": return "V3"
    if p.startswith("protein"): return "P" + p.replace("protein", "")
    return p

_SIM = re.compile(r"^LAF_1_V2_(.+?)_(\d+)_(\d+)_298$")
_REC = re.compile(r"^Record:\s*LAF_1_V2_(.+?)_(\d+)_(\d+)_298$")
_CONC = re.compile(r"^(\S+):\s*dense=([\d.]+)\S*\s*mM\s+dilute=([\d.]+)\S*\s*µM")
_DEMIX = re.compile(r"^D \(demixing index\):\s*([\d.]+)\s*±\s*([\d.]+)")

def load(out_file=OUT_FILE, f1_file=F1_FILE):
    lines = open(out_file, encoding="utf-8").read().splitlines()
    rows = {}
    cur = None
    for ln in lines:
        s = ln.strip()
        m = _SIM.match(s)
        if m:
            cur = (_partner_to_label(m.group(1)), int(m.group(2)))
            rows.setdefault(cur, {}); continue
        if cur:
            c = _CONC.match(s)
            if c:
                who, dense, dil = c.group(1), float(c.group(2)), float(c.group(3))
                if who == "V2":
                    rows[cur]["V2_dense"] = dense; rows[cur]["V2_dil"] = dil
                else:
                    rows[cur]["B_dense"] = dense; rows[cur]["B_dil"] = dil
    cur = None
    for ln in lines:
        s = ln.strip()
        m = _REC.match(s)
        if m:
            cur = (_partner_to_label(m.group(1)), int(m.group(2))); continue
        if cur:
            d = _DEMIX.match(s)
            if d:
                rows.setdefault(cur, {})["D"] = float(d.group(1))
                rows[cur]["D_sem"] = float(d.group(2))
                cur = None
    df = pd.DataFrame([{"label": k[0], "pct": k[1], **v} for k, v in rows.items()])
    df["K_B"] = df["B_dense"] / (df["B_dil"] / 1000.0)
    f1 = pd.read_csv(f1_file)[["label", "SCD"]]
    df = df.merge(f1, on="label", how="left")
    return df.sort_values(["SCD", "pct"]).reset_index(drop=True)
