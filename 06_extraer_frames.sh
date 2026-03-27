#!/bin/bash
# ==========================================================
# 06_extraer_frames.sh
# Extrae frames PDB representativos para visualizacion (ChimeraX)
# Frame: 175 ns (punto medio de la ventana MM-GBSA 150-200 ns)
#
# Output: Analisis_Global/Frames/frame_175ns_{ph}_{replica}.pdb
#
# Uso: bash 06_extraer_frames.sh
# ==========================================================

set -uo pipefail

BASE="$(cd "$(dirname "$0")" && pwd)"
REPLICAS="Replica_1 Replica_2 Replica_3"
PHs="ph50 ph70 ph80"
ANALISIS_DIR="$BASE/Analisis_Global"
FRAMES_DIR="$ANALISIS_DIR/Frames"
DUMP_TIME=175000  # 175 ns en ps

mkdir -p "$FRAMES_DIR"

echo ""
echo "=========================================================="
echo "  06_extraer_frames.sh — Frames para visualizacion"
echo "  Tiempo: 175 ns (punto medio ventana MM-GBSA)"
echo "=========================================================="

for REPLICA in $REPLICAS; do
    REP_NUM=${REPLICA##*_}
    for PH in $PHs; do

        SIM_DIR="$BASE/$REPLICA/$PH"
        OUT_DIR="$ANALISIS_DIR/$REPLICA/$PH"
        NDX="$OUT_DIR/index_analisis.ndx"
        XTC="$OUT_DIR/Visuales/centrado.xtc"
        TPR="$SIM_DIR/step7_1.tpr"

        OUT_PDB="$FRAMES_DIR/frame_175ns_${PH}_R${REP_NUM}.pdb"

        if [ -f "$OUT_PDB" ]; then
            echo "  [OK] $OUT_PDB ya existe."
            continue
        fi

        if [ ! -f "$TPR" ] || [ ! -f "$XTC" ]; then
            echo "  [SKIP] Faltan archivos para $REPLICA/$PH"
            continue
        fi

        echo "  Extrayendo: $REPLICA / $PH ..."
        printf "Protein\nSystem\n" | gmx trjconv \
            -s "$TPR" \
            -f "$XTC" \
            -o "$OUT_PDB" \
            -dump $DUMP_TIME \
            -pbc mol \
            -center 2>/dev/null

        if [ -f "$OUT_PDB" ]; then
            echo "  [OK] $OUT_PDB"
        else
            echo "  [ERROR] No se pudo extraer frame."
        fi

    done
done

echo ""
echo "=========================================================="
echo "  Frames generados:"
ls -lh "$FRAMES_DIR"/*.pdb 2>/dev/null || echo "  (ninguno)"
echo ""
echo "  SIGUIENTE PASO: Rscript 07_graficos.R"
echo "=========================================================="
