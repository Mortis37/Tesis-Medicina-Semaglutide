#!/bin/bash
# ==========================================================
# 03_mmgbsa_energia.sh
# ENERGIA LIBRE DE UNION — MM-GBSA (gmx_MMPBSA)
# Sistema: Semaglutida (PROA) + GLP-1R (PROB) + POPC
#
# Prerequisito: 02_analisis_estructural.sh completado
#   Necesita: centrado.xtc + index_analisis.ndx
#
# FRAME_INICIO=3000 (150 ns) — justificado por convergencia
# igb=5 (GB Onufriev modificado)
# saltcon=0.150 M NaCl (condicion fisiologica)
#
# Nota: memopt no disponible en gmx_MMPBSA 1.6.4 (requiere >= 1.6.7)
# El sitio de union semaglutida-GLP1R es ECD extracelular,
# no dominio TM — impacto de membrana en deltaG es minimo.
# mcenter se mide y documenta pero no se pasa al calculo.
#
# Uso: bash 03_mmgbsa_energia.sh
# ==========================================================

set -uo pipefail
export LC_NUMERIC=C

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA

export OMP_NUM_THREADS=20

FRAME_INICIO=3000
FRAME_FIN=4000
INTERVALO=10
N_FRAMES=$(( (FRAME_FIN - FRAME_INICIO) / INTERVALO ))

BASE="$(cd "$(dirname "$0")" && pwd)"
REPLICAS="Replica_1 Replica_2 Replica_3"
PHs="ph50 ph70 ph80"
ANALISIS_DIR="$BASE/Analisis_Global"
INICIO_GLOBAL=$(date +%s)

echo ""
echo "=========================================================="
echo "  03_mmgbsa_energia.sh — MM-GBSA"
echo "  Frames        : $FRAME_INICIO a $FRAME_FIN (cada $INTERVALO)"
echo "  Frames totales: $N_FRAMES por condicion (= 150-200 ns)"
echo "  igb=5  saltcon=0.150 M"
echo "  Inicio        : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================================="

for REPLICA in $REPLICAS; do
    for PH in $PHs; do

        T_INICIO=$(date +%s)
        SIM_DIR="$BASE/$REPLICA/$PH"
        OUT_DIR="$ANALISIS_DIR/$REPLICA/$PH"
        ENERGIA_DIR="$OUT_DIR/Energia"
        LOG="$ENERGIA_DIR/mmgbsa.log"

        echo ""
        echo "=========================================================="
        echo "  >>> [$REPLICA / $PH]  $(date '+%H:%M:%S')"
        echo "=========================================================="

        NDX="$OUT_DIR/index_analisis.ndx"
        XTC="$OUT_DIR/Visuales/centrado.xtc"
        TPR="$SIM_DIR/step7_1.tpr"
        TOP="$SIM_DIR/topol.top"

        # Verificar archivos
        FALTAN=0
        for f in "$NDX" "$XTC" "$TPR" "$TOP"; do
            [ ! -f "$f" ] && echo "  [SKIP] Falta: $f" && FALTAN=1
        done
        [ $FALTAN -eq 1 ] && continue

        # Si ya existe, mostrar resultado y saltar
        RESULTADO="$ENERGIA_DIR/FINAL_RESULTS_MMGBSA.dat"
        if [ -f "$RESULTADO" ]; then
            echo "  [OK] Ya calculado."
            grep "TOTAL" "$RESULTADO" | tail -1
            continue
        fi

        mkdir -p "$ENERGIA_DIR"

        # Medir centro de membrana (documentacion)
        MCENTER_FILE="$ENERGIA_DIR/membrane_center.txt"
        if [ ! -f "$MCENTER_FILE" ]; then
            MCENTER=$(grep "^ATOM" "$SIM_DIR/step5_input.pdb" | \
                awk '$3=="P" || $3==" P " {
                    z=substr($0,47,8)+0; sum+=z; n++
                } END {
                    if(n>0) printf "%.2f", sum/n;
                    else print "0.00"
                }')
            echo "$MCENTER" > "$MCENTER_FILE"
            echo "  [INFO] mcenter: ${MCENTER} A (documentado)"
        else
            MCENTER=$(cat "$MCENTER_FILE")
            echo "  [INFO] mcenter (cached): ${MCENTER} A"
        fi

        # Input MM-GBSA
        cat > "$ENERGIA_DIR/mmpbsa.in" << EOF
&general
  sys_name="${REPLICA}_${PH}",
  startframe=$FRAME_INICIO,
  endframe=$FRAME_FIN,
  interval=$INTERVALO,
  verbose=2,
/
&gb
  igb=5,
  saltcon=0.150,
/
EOF

        echo "  [INFO] Lanzando gmx_MMPBSA (OMP=$OMP_NUM_THREADS)..."
        cd "$ENERGIA_DIR"

        gmx_MMPBSA -O \
            -i  mmpbsa.in \
            -cs "$TPR" \
            -ci "$NDX" \
            -cg Receptor Ligando \
            -ct "$XTC" \
            -cp "$TOP" \
            -o  FINAL_RESULTS_MMGBSA.dat \
            -eo FINAL_CLEAN.csv \
            -nogui 2>&1 | tee -a "$LOG"

        cd "$BASE"

        if [ ! -f "$RESULTADO" ]; then
            echo "  [ERROR] No se genero FINAL_RESULTS_MMGBSA.dat — ver $LOG"
            continue
        fi

        T_FIN=$(date +%s)
        T_MIN=$(( (T_FIN - T_INICIO) / 60 ))

        echo ""
        echo "  [$REPLICA / $PH] DeltaG:"
        grep "TOTAL" "$RESULTADO" | tail -1
        echo "  Tiempo: ${T_MIN}m"

    done
done

# Resumen final
T_FIN_GLOBAL=$(date +%s)
T_TOTAL=$(( T_FIN_GLOBAL - INICIO_GLOBAL ))

echo ""
echo "=========================================================="
echo "  03_mmgbsa_energia.sh FINALIZADO"
echo "  Tiempo total: $(( T_TOTAL / 60 ))m"
echo ""
printf "  %-14s %-6s  %s\n" "REPLICA" "pH" "DeltaG (kcal/mol)"
echo "  $(printf '%0.s-' {1..55})"
for REPLICA in $REPLICAS; do
    for PH in $PHs; do
        R="$ANALISIS_DIR/$REPLICA/$PH/Energia/FINAL_RESULTS_MMGBSA.dat"
        if [ -f "$R" ]; then
            LINEA=$(grep "TOTAL" "$R" | tail -1 | awk '{printf "%.2f +/- %.2f", $2, $4}')
            printf "  %-14s %-6s  %s kcal/mol\n" "$REPLICA" "$PH" "$LINEA"
        else
            printf "  %-14s %-6s  [NO CALCULADO]\n" "$REPLICA" "$PH"
        fi
    done
    echo ""
done
echo "  SIGUIENTE PASO: bash 04_mmgbsa_decomp.sh"
echo "=========================================================="
