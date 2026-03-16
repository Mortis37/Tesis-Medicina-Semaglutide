#!/usr/bin/env python3
# ==========================================================
# 4_analisis_hotspots.py
# Parseo de descomposicion MM-GBSA por residuo
# Sistema: Semaglutida (PROA/A) + GLP-1R (PROB/B)
#
# Input:  FINAL_DECOMP_MMPBSA.dat de cada condicion
# Output: CSVs de hotspots + tablas resumen
#
# Hotspot: residuo con TOTAL_avg < -1.5 kcal/mol
# Columnas parseadas: Internal, VdW, Elec, PolarSolv, NonPolar, TOTAL
# ==========================================================

import os
import re
import pandas as pd
import numpy as np

BASE = os.path.expanduser(
    "~/Documentos/Proyecto_Tesis_Renato_MD/Analisis_Global"
)
OUTPUT_DIR = os.path.expanduser(
    "~/Documentos/Proyecto_Tesis_Renato_MD/Resultados_Hotspots"
)
os.makedirs(OUTPUT_DIR, exist_ok=True)

REPLICAS = ["Replica_1", "Replica_2", "Replica_3"]
PHS      = ["ph50", "ph70", "ph80"]
PH_LABEL = {"ph50": "pH 5.0", "ph70": "pH 7.4", "ph80": "pH 8.0"}

HOTSPOT_UMBRAL = -1.5  # kcal/mol


# ----------------------------------------------------------
# PARSER
# ----------------------------------------------------------
def parse_decomp(filepath):
    """
    Lee FINAL_DECOMP_MMPBSA.dat y extrae la seccion
    'Total Energy Decomposition' del bloque Complex.
    Devuelve DataFrame con columnas:
      residue, chain, resname, resnum,
      internal_avg, vdw_avg, elec_avg,
      polar_avg, nonpolar_avg, total_avg, total_sd
    """
    rows = []
    in_complex  = False
    in_total    = False
    header_done = False

    with open(filepath, "r") as f:
        for line in f:
            line = line.rstrip("\n")

            # Detectar bloque Complex
            if line.strip() == "Complex:":
                in_complex = True
                continue
            # Salir del bloque Complex al llegar a Receptor o Ligand solo
            if in_complex and line.strip() in ("Receptor:", "Ligand:"):
                break

            if not in_complex:
                continue

            # Detectar subseccion Total Energy Decomposition
            if "Total Energy Decomposition" in line:
                in_total    = True
                header_done = False
                continue

            if not in_total:
                continue

            # Saltar las 2 lineas de cabecera
            if not header_done:
                if line.startswith("Residue,"):
                    continue          # primera cabecera
                if line.startswith(",Avg."):
                    header_done = True
                    continue

            # Linea de datos: L:A:HSD:7,... o R:B:THR:29,...
            if not re.match(r"^[LR]:[AB]:", line):
                continue

            parts = line.split(",")
            if len(parts) < 18:
                continue

            res_id  = parts[0]          # L:A:HSD:7
            tokens  = res_id.split(":")
            if len(tokens) < 4:
                continue

            mol_type = tokens[0]        # L=ligando, R=receptor
            chain    = tokens[1]        # A o B
            resname  = tokens[2]        # HSD, GLU, etc.
            resnum   = int(tokens[3])

            try:
                internal_avg  = float(parts[1])
                internal_sd   = float(parts[2])
                vdw_avg       = float(parts[4])
                vdw_sd        = float(parts[5])
                elec_avg      = float(parts[7])
                elec_sd       = float(parts[8])
                polar_avg     = float(parts[10])
                polar_sd      = float(parts[11])
                nonpolar_avg  = float(parts[13])
                nonpolar_sd   = float(parts[14])
                total_avg     = float(parts[16])
                total_sd      = float(parts[17])
            except (ValueError, IndexError):
                continue

            rows.append({
                "residue"     : res_id,
                "mol_type"    : mol_type,
                "chain"       : chain,
                "resname"     : resname,
                "resnum"      : resnum,
                "internal_avg": internal_avg,
                "vdw_avg"     : vdw_avg,
                "vdw_sd"      : vdw_sd,
                "elec_avg"    : elec_avg,
                "elec_sd"     : elec_sd,
                "polar_avg"   : polar_avg,
                "nonpolar_avg": nonpolar_avg,
                "total_avg"   : total_avg,
                "total_sd"    : total_sd,
            })

    return pd.DataFrame(rows)


# ----------------------------------------------------------
# LOOP PRINCIPAL
# ----------------------------------------------------------
all_data = []

print("\n" + "="*60)
print("  4_analisis_hotspots.py — Descomposicion MM-GBSA")
print("  Umbral hotspot: TOTAL < {:.1f} kcal/mol".format(HOTSPOT_UMBRAL))
print("="*60)

for replica in REPLICAS:
    for ph in PHS:
        dat_file = os.path.join(
            BASE, replica, ph, "Decomp", "FINAL_DECOMP_MMPBSA.dat"
        )
        if not os.path.isfile(dat_file):
            print(f"  [SKIP] No encontrado: {replica}/{ph}")
            continue

        df = parse_decomp(dat_file)
        if df.empty:
            print(f"  [WARN] Sin datos: {replica}/{ph}")
            continue

        df["replica"] = replica
        df["ph"]      = ph
        df["ph_label"]= PH_LABEL[ph]
        df["rep_num"] = int(replica.split("_")[1])

        all_data.append(df)
        n_hot = (df[df["mol_type"]=="R"]["total_avg"] < HOTSPOT_UMBRAL).sum()
        print(f"  [OK] {replica}/{ph}: {len(df)} residuos, {n_hot} hotspots receptor")

# Concatenar todo
master = pd.concat(all_data, ignore_index=True)
print(f"\n  Total filas: {len(master)}")


# ----------------------------------------------------------
# TABLA 1: Hotspots del receptor — promedio entre replicas
# Residuos con total_avg < HOTSPOT_UMBRAL en >= 2 replicas
# ----------------------------------------------------------
rec = master[master["mol_type"] == "R"].copy()

# Promedio por condicion de pH (entre las 3 replicas)
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
rec_mean["residue_label"] = (
    rec_mean["resname"] + rec_mean["resnum"].astype(str)
)

# Hotspots: total_mean < umbral
hotspots = rec_mean[rec_mean["total_mean"] < HOTSPOT_UMBRAL].copy()
hotspots = hotspots.sort_values("total_mean")

# Guardar tabla completa de hotspots
out_hot = os.path.join(OUTPUT_DIR, "hotspots_receptor_por_pH.csv")
hotspots.to_csv(out_hot, index=False, float_format="%.3f")
print(f"\n  [SAVE] {out_hot}")

# Tabla pivote: residuos vs pH
pivot = hotspots.pivot_table(
    index=["resname", "resnum", "residue_label"],
    columns="ph_label",
    values="total_mean",
    aggfunc="first"
).reset_index()
pivot = pivot.sort_values(
    [c for c in ["pH 5.0","pH 7.4","pH 8.0"] if c in pivot.columns][:1]
)
out_pivot = os.path.join(OUTPUT_DIR, "hotspots_pivote_pH.csv")
pivot.to_csv(out_pivot, index=False, float_format="%.3f")
print(f"  [SAVE] {out_pivot}")


# ----------------------------------------------------------
# TABLA 2: Top 15 hotspots por pH
# ----------------------------------------------------------
print("\n" + "-"*60)
for ph in PHS:
    label = PH_LABEL[ph]
    sub = hotspots[hotspots["ph"] == ph].nsmallest(15, "total_mean")
    print(f"\n  Top 15 hotspots receptor — {label}:")
    print(f"  {'Residuo':<10} {'Total(kcal/mol)':>16} {'±SD':>8} {'VdW':>8} {'Elec':>8}")
    print(f"  {'-'*52}")
    for _, row in sub.iterrows():
        print(f"  {row['residue_label']:<10} {row['total_mean']:>16.3f} "
              f"{row['total_sd']:>8.3f} {row['vdw_mean']:>8.3f} {row['elec_mean']:>8.3f}")


# ----------------------------------------------------------
# TABLA 3: Contribucion del ligando (semaglutida) por residuo
# ----------------------------------------------------------
lig = master[master["mol_type"] == "L"].copy()
lig_mean = (
    lig.groupby(["ph", "ph_label", "resname", "resnum"])
    .agg(
        total_mean=("total_avg", "mean"),
        total_sd  =("total_avg", "std"),
        vdw_mean  =("vdw_avg",   "mean"),
        elec_mean =("elec_avg",  "mean"),
    )
    .reset_index()
)
lig_mean["residue_label"] = lig_mean["resname"] + lig_mean["resnum"].astype(str)
lig_mean = lig_mean.sort_values("total_mean")

out_lig = os.path.join(OUTPUT_DIR, "contribucion_ligando_por_pH.csv")
lig_mean.to_csv(out_lig, index=False, float_format="%.3f")
print(f"\n  [SAVE] {out_lig}")


# ----------------------------------------------------------
# TABLA 4: Resumen DeltaG promedio por pH (entre replicas)
# ----------------------------------------------------------
print("\n" + "-"*60)
print("\n  Resumen DeltaG global (del script 2):")
deltag = {
    "pH 5.0": [-143.67, -91.97, -81.16],
    "pH 7.4": [-178.43, -120.08, -95.26],
    "pH 8.0": [-181.76, -142.87, -129.41],
}
print(f"  {'pH':<8} {'Media':>10} {'SD':>8} {'SEM':>8}")
print(f"  {'-'*36}")
for label, vals in deltag.items():
    m = np.mean(vals)
    s = np.std(vals, ddof=1)
    sem = s / np.sqrt(len(vals))
    print(f"  {label:<8} {m:>10.2f} {s:>8.2f} {sem:>8.2f}")


# ----------------------------------------------------------
# TABLA 5: Residuos "consistentes" — hotspot en los 3 pH
# ----------------------------------------------------------
hot_all = set()
for ph in PHS:
    sub = set(
        hotspots[hotspots["ph"]==ph]["residue_label"].tolist()
    )
    if not hot_all:
        hot_all = sub
    else:
        hot_all &= sub

print(f"\n  Hotspots presentes en los 3 pH ({len(hot_all)} residuos):")
consec = pivot[pivot["residue_label"].isin(hot_all)].copy()
if not consec.empty:
    ph_cols = [c for c in ["pH 5.0","pH 7.4","pH 8.0"] if c in consec.columns]
    consec["mean_all_pH"] = consec[ph_cols].mean(axis=1)
    consec = consec.sort_values("mean_all_pH")
    for _, row in consec.iterrows():
        vals = " | ".join(
            f"{c}: {row[c]:.2f}" for c in ph_cols if c in row.index
        )
        print(f"    {row['residue_label']:<10}  {vals}")
    out_consec = os.path.join(OUTPUT_DIR, "hotspots_consistentes_3pH.csv")
    consec.to_csv(out_consec, index=False, float_format="%.3f")
    print(f"\n  [SAVE] {out_consec}")


# ----------------------------------------------------------
# TABLA MAESTRA — todos los residuos del receptor
# ----------------------------------------------------------
out_master = os.path.join(OUTPUT_DIR, "decomposicion_completa_receptor.csv")
rec_mean.to_csv(out_master, index=False, float_format="%.4f")
print(f"\n  [SAVE] Tabla maestra receptor: {out_master}")

out_master_all = os.path.join(OUTPUT_DIR, "decomposicion_completa_todos.csv")
master.to_csv(out_master_all, index=False, float_format="%.4f")
print(f"  [SAVE] Tabla maestra completa: {out_master_all}")

print("\n" + "="*60)
print("  COMPLETADO — archivos en:")
print(f"  {OUTPUT_DIR}")
print("="*60 + "\n")
