#!/bin/bash
# ==========================================================
# 01_produccion_md.sh
# Corrida de produccion — GROMACS / CHARMM-GUI
# Sistema: Semaglutida + GLP-1R + POPC (CHARMM36m)
#
# Lanza gmx mdrun para las 3 replicas x 3 condiciones de pH
# Los archivos .tpr fueron generados por CHARMM-GUI
# Tiempo de produccion: 200 ns por replica (100M pasos, dt=2fs)
#
# Prerequisitos:
#   - GROMACS 2025+ instalado
#   - Estructura de carpetas:
#       Replica_{1,2,3}/ph{50,70,80}/step7_1.tpr
#   - GPU recomendada (ej: RTX 4090)
#
# Uso: bash 01_produccion_md.sh
# ==========================================================

set -uo pipefail

BASE="$(cd "$(dirname "$0")" && pwd)"
REPLICAS="Replica_1 Replica_2 Replica_3"
PHs="ph50 ph70 ph80"
INICIO_GLOBAL=$(date +%s)

echo ""
echo "=========================================================="
echo "  01_produccion_md.sh — PRODUCCION MD"
echo "  Replicas : $REPLICAS"
echo "  pHs      : $PHs"
echo "  Inicio   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================================="

for REPLICA in $REPLICAS; do
    for PH in $PHs; do

        T_INICIO=$(date +%s)
        SIM_DIR="$BASE/$REPLICA/$PH"

        echo ""
        echo "=========================================================="
        echo "  >>> [$REPLICA / $PH]  $(date '+%H:%M:%S')"
        echo "=========================================================="

        # Verificar TPR
        TPR="$SIM_DIR/step7_1.tpr"
        if [ ! -f "$TPR" ]; then
            echo "  [SKIP] Falta $TPR"
            continue
        fi

        # Si ya existe la trayectoria completa, saltar
        XTC="$SIM_DIR/step7_1.xtc"
        GRO="$SIM_DIR/step7_1.gro"
        if [ -f "$XTC" ] && [ -f "$GRO" ]; then
            echo "  [OK] Trayectoria ya completa — omitiendo."
            continue
        fi

        # Verificar checkpoint para continuar corrida interrumpida
        CPT="$SIM_DIR/step7_1.cpt"
        if [ -f "$CPT" ]; then
            echo "  [INFO] Checkpoint encontrado — continuando corrida..."
            cd "$SIM_DIR"
            gmx mdrun \
                -s step7_1.tpr \
                -deffnm step7_1 \
                -cpi step7_1.cpt \
                -v \
                -ntomp "${OMP_NUM_THREADS:-8}" \
                -pin on
        else
            echo "  [INFO] Lanzando produccion desde cero..."
            cd "$SIM_DIR"
            gmx mdrun \
                -s step7_1.tpr \
                -deffnm step7_1 \
                -v \
                -ntomp "${OMP_NUM_THREADS:-8}" \
                -pin on
        fi

        cd "$BASE"

        T_FIN=$(date +%s)
        T_ELAPSED=$(( T_FIN - T_INICIO ))
        T_HR=$(( T_ELAPSED / 3600 ))
        T_MIN=$(( (T_ELAPSED % 3600) / 60 ))

        echo ""
        echo "  [$REPLICA / $PH] completado en ${T_HR}h ${T_MIN}m"
        echo "----------------------------------------------------------"

    done
done

# Resumen final
T_FIN_GLOBAL=$(date +%s)
T_TOTAL=$(( T_FIN_GLOBAL - INICIO_GLOBAL ))
T_HR_TOTAL=$(( T_TOTAL / 3600 ))
T_MIN_TOTAL=$(( (T_TOTAL % 3600) / 60 ))

echo ""
echo "=========================================================="
echo "  01_produccion_md.sh FINALIZADO"
echo "  Tiempo total: ${T_HR_TOTAL}h ${T_MIN_TOTAL}m"
echo ""
echo "  Estado de trayectorias:"
for REPLICA in $REPLICAS; do
    for PH in $PHs; do
        XTC="$BASE/$REPLICA/$PH/step7_1.xtc"
        if [ -f "$XTC" ]; then
            SIZE=$(du -h "$XTC" | cut -f1)
            printf "    %-12s %-6s  OK (%s)\n" "$REPLICA" "$PH" "$SIZE"
        else
            printf "    %-12s %-6s  [NO COMPLETADO]\n" "$REPLICA" "$PH"
        fi
    done
done
echo ""
echo "  SIGUIENTE PASO: bash 02_analisis_estructural.sh"
echo "=========================================================="
