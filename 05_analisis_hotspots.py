#!/usr/bin/env python3
"""
05_analisis_hotspots.py
Parseo de descomposicion MM-GBSA por residuo
Sistema: Semaglutida (PROA/A) + GLP-1R (PROB/B)

Input:  FINAL_DECOMP_MMPBSA.dat de cada condicion
Output: CSVs de hotspots + tablas resumen

Hotspot: residuo con TOTAL_avg < -1.5 kcal/mol

Uso: python3 05_analisis_hotspots.py
     (ejecutar desde la raiz del proyecto)
"""

import os
import re
import pandas as pd
import numpy as np

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Analisis_Global")
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Resultados_Hotspots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

REPLICAS = ["Replica_1", "Replica_2", "Replica_3"]
PHS      = ["ph50", "ph70", "ph80"]
PH_LABEL = {"ph50": "pH 5.0", "ph70": "pH 7.4", "ph80": "pH 8.0"}

HOTSPOT_UMBRAL = -1.5  # kcal/mol


def parse_decomp(filepath):
    """
    Lee FINAL_DECOMP_MMPBSA.dat y extrae la seccion
    'Total Energy Decomposition' del bloque Complex.
    """
    rows = []
    in_complex  = False
    in_total    = False
    header_done = False

    with open(filepath, "r") as f:
        for line in f:
            line = line.rstrip("\n")

            if line.strip() == "Complex:":
                in_complex = True
                continue
            if in_complex and line.strip() in ("Receptor:", "Ligand:"):
                break
            if not in_complex:
                continue

            if "Total Energy Decomposition" in line:
                in_total    = True
                header_done = False
                continue
            if not in_total:
                continue

            if not header_done:
                if line.startswith("Residue,"):
                    continue
                if line.startswith(",Avg."):
                    header_done = True
                    continue

            if not re.match(r"^[LR]:[AB]:", line):
                continue

            parts = line.split(",")
            if len(parts) < 18:
                continue

            res_id  = parts[0]
            tokens  = res_id.split(":")
            if len(tokens) < 4:
                continue

            mol_type = tokens[0]
            chain    = tokens[1]
            resname  = tokens[2]
            resnum   = int(tokens[3])

            try:
                rows.append({
                    "residue"     : res_id,
                    "mol_type"    : mol_type,
                    "chain"       : chain,
                    "resname"     : resname,
                    "resnum"      : resnum,
                    "internal_avg": float(parts[1]),
                    "vdw_avg"     : float(parts[4]),
                    "vdw_sd"      : float(parts[5]),
                    "elec_avg"    : float(parts[7]),
                    "elec_sd"     : float(parts[8]),
                    "polar_avg"   : float(parts[10]),
                    "nonpolar_avg": float(parts[13]),
                    "total_avg"   : float(parts[16]),
                    "total_sd"    : float(parts[17]),
                })
            except (ValueError, IndexError):
                continue

    return pd.DataFrame(rows)


# ==========================================================
# LOOP PRINCIPAL
# ==========================================================
all_data = []

print("\n" + "=" * 60)
print("  05_analisis_hotspots.py")
print("  Umbral hotspot: TOTAL < {:.1f} kcal/mol".format(HOTSPOT_UMBRAL))
print("=" * 60)

for replica in REPLICAS:
    for ph in PHS:
        dat_file = os.path.join(BASE, replica, ph, "Decomp", "FINAL_DECOMP_MMPBSA.dat")
        if not os.path.isfile(dat_file):
            print(f"  [SKIP] {replica}/{ph}")
            continue

        df = parse_decomp(dat_file)
        if df.empty:
            print(f"  [WARN] Sin datos: {replica}/{ph}")
            continue

        df["replica"]  = replica
        df["ph"]       = ph
        df["ph_label"] = PH_LABEL[ph]
        df["rep_num"]  = int(replica.split("_")[1])
        all_data.append(df)

        n_hot = (df[df["mol_type"] == "R"]["total_avg"] < HOTSPOT_UMBRAL).sum()
        print(f"  [OK] {replica}/{ph}: {len(df)} residuos, {n_hot} hotspots receptor")

master = pd.concat(all_data, ignore_index=True)
print(f"\n  Total filas: {len(master)}")


# ==========================================================
# TABLA 1: Hotspots del receptor (promedio entre replicas)
# ==========================================================
rec = master[master["mol_type"] == "R"].copy()
rec_mean = (
    rec.groupby(["ph", "ph_label", "resname", "resnum"])
    .agg(
        total_mean=("total_avg", "mean"),
        total_sd  =("total_avg", "std"),
        vdw_mean  =("vdw_avg",   "mean"),
        elec_mean =("elec_avg",  "mean"),
        polar_mean=("polar_avg", "mean"),
        n_replicas=("total_avg", "count"),
    )
    .reset_index()
)
rec_mean["residue_label"] = rec_mean["resname"] + rec_mean["resnum"].astype(str)

hotspots = rec_mean[rec_mean["total_mean"] < HOTSPOT_UMBRAL].sort_values("total_mean")

out_hot = os.path.join(OUTPUT_DIR, "hotspots_receptor_por_pH.csv")
hotspots.to_csv(out_hot, index=False, float_format="%.3f")
print(f"\n  [SAVE] {out_hot}")

# Tabla pivote: residuos vs pH
pivot = hotspots.pivot_table(
    index=["resname", "resnum", "residue_label"],
    columns="ph_label", values="total_mean", aggfunc="first"
).reset_index()
pivot = pivot.sort_values(
    [c for c in ["pH 5.0", "pH 7.4", "pH 8.0"] if c in pivot.columns][:1]
)
out_pivot = os.path.join(OUTPUT_DIR, "hotspots_pivote_pH.csv")
pivot.to_csv(out_pivot, index=False, float_format="%.3f")
print(f"  [SAVE] {out_pivot}")


# ==========================================================
# TABLA 2: Top 15 hotspots por pH
# ==========================================================
print("\n" + "-" * 60)
for ph in PHS:
    label = PH_LABEL[ph]
    sub = hotspots[hotspots["ph"] == ph].nsmallest(15, "total_mean")
    print(f"\n  Top 15 hotspots receptor — {label}:")
    print(f"  {'Residuo':<10} {'Total(kcal/mol)':>16} {'SD':>8} {'VdW':>8} {'Elec':>8}")
    print(f"  {'-' * 52}")
    for _, row in sub.iterrows():
        print(f"  {row['residue_label']:<10} {row['total_mean']:>16.3f} "
              f"{row['total_sd']:>8.3f} {row['vdw_mean']:>8.3f} {row['elec_mean']:>8.3f}")


# ==========================================================
# TABLA 3: Contribucion del ligando
# ==========================================================
lig = master[master["mol_type"] == "L"].copy()
lig_mean = (
    lig.groupby(["ph", "ph_label", "resname", "resnum"])
    .agg(total_mean=("total_avg", "mean"), total_sd=("total_avg", "std"),
         vdw_mean=("vdw_avg", "mean"), elec_mean=("elec_avg", "mean"))
    .reset_index()
)
lig_mean["residue_label"] = lig_mean["resname"] + lig_mean["resnum"].astype(str)
lig_mean = lig_mean.sort_values("total_mean")

out_lig = os.path.join(OUTPUT_DIR, "contribucion_ligando_por_pH.csv")
lig_mean.to_csv(out_lig, index=False, float_format="%.3f")
print(f"\n  [SAVE] {out_lig}")


# ==========================================================
# TABLA 4: Hotspots consistentes en los 3 pH
# ==========================================================
hot_all = set()
for ph in PHS:
    sub = set(hotspots[hotspots["ph"] == ph]["residue_label"].tolist())
    hot_all = sub if not hot_all else hot_all & sub

print(f"\n  Hotspots en los 3 pH ({len(hot_all)} residuos):")
consec = pivot[pivot["residue_label"].isin(hot_all)].copy()
if not consec.empty:
    ph_cols = [c for c in ["pH 5.0", "pH 7.4", "pH 8.0"] if c in consec.columns]
    consec["mean_all_pH"] = consec[ph_cols].mean(axis=1)
    consec = consec.sort_values("mean_all_pH")
    for _, row in consec.iterrows():
        vals = " | ".join(f"{c}: {row[c]:.2f}" for c in ph_cols)
        print(f"    {row['residue_label']:<10}  {vals}")
    out_consec = os.path.join(OUTPUT_DIR, "hotspots_consistentes_3pH.csv")
    consec.to_csv(out_consec, index=False, float_format="%.3f")
    print(f"\n  [SAVE] {out_consec}")


# ==========================================================
# TABLAS MAESTRAS
# ==========================================================
out_master = os.path.join(OUTPUT_DIR, "decomposicion_completa_receptor.csv")
rec_mean.to_csv(out_master, index=False, float_format="%.4f")
print(f"\n  [SAVE] {out_master}")

out_master_all = os.path.join(OUTPUT_DIR, "decomposicion_completa_todos.csv")
master.to_csv(out_master_all, index=False, float_format="%.4f")
print(f"  [SAVE] {out_master_all}")

print("\n" + "=" * 60)
print("  COMPLETADO — archivos en:")
print(f"  {OUTPUT_DIR}")
print("  SIGUIENTE PASO: bash 06_extraer_frames.sh")
print("=" * 60 + "\n")
