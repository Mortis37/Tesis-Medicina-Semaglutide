#!/bin/bash
# ==========================================================
# 2_mmgbsa_energia.sh
# ENERGIA LIBRE DE UNION — MM-GBSA (gmx_MMPBSA)
# Sistema: Semaglutida (PROA) + GLP-1R (PROB) + POPC
#
# Prerequisito: 1_analisis_estructural_v6.sh completado
#   Necesita: centrado.xtc + index_analisis.ndx
#
# FRAME_INICIO=3000 justificado por 1.2_verificar_convergencia.sh
# igb=5 (GB Onufriev modificado)
# memopt NO disponible en gmx_MMPBSA 1.6.4 — requiere >= 1.6.7
# Justificacion: sitio de union semaglutida-GLP1R es ECD extracelular,
#   no dominio TM — impacto de membrana en deltaG de union es minimo.
# mcenter se mide y documenta en log pero no se pasa al calculo.
# saltcon=0.150 M NaCl (condicion fisiologica)
#
# Tiempo estimado: ~1-2 horas por replica/pH (i9-14900K, OMP=20)
# ==========================================================

set -uo pipefail
export LC_NUMERIC=C

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gmxMMPBSA

# ----------------------------------------------------------
# PARALELISMO — i9-14900K (24 cores fisicos)
# gmx_MMPBSA serial usa OpenMP internamente.
# OMP_NUM_THREADS=20 aprovecha los P-cores del i9-14900K.
# ----------------------------------------------------------
export OMP_NUM_THREADS=20

FRAME_INICIO=3000
FRAME_FIN=4000
INTERVALO=10
N_FRAMES=$(( (FRAME_FIN - FRAME_INICIO) / INTERVALO ))

REPLICAS="Replica_1 Replica_2 Replica_3"
PHs="ph50 ph70 ph80"
BASE="$HOME/Documentos/Proyecto_Tesis_Renato_MD"
ANALISIS_DIR="$BASE/Analisis_Global"
INICIO_GLOBAL=$(date +%s)

echo ""
echo "=========================================================="
echo "  2_mmgbsa_energia.sh — MM-GBSA"
echo "  Replicas      : $REPLICAS"
echo "  Frames        : $FRAME_INICIO a $FRAME_FIN (cada $INTERVALO)"
echo "  Frames totales: $N_FRAMES por condicion (= 150-200 ns)"
echo "  igb=5  saltcon=0.150 M  (memopt no disponible en v1.6.4)"
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

        # Verificar archivos necesarios
        FALTAN=0
        for f in "$NDX" "$XTC" "$TPR" "$TOP"; do
            if [ ! -f "$f" ]; then
                echo "  [SKIP] Falta: $f"
                FALTAN=1
            fi
        done
        [ $FALTAN -eq 1 ] && echo "  [SKIP] $REPLICA/$PH omitido." && continue

        # Si ya existe, mostrar resultado y saltar
        RESULTADO="$ENERGIA_DIR/FINAL_RESULTS_MMGBSA.dat"
        if [ -f "$RESULTADO" ]; then
            echo "  [OK] Ya calculado anteriormente."
            echo "  Resultado guardado:"
            grep "TOTAL" "$RESULTADO" | tail -1
            continue
        fi

        mkdir -p "$ENERGIA_DIR"

        # --------------------------------------------------
        # PASO PREVIO: medir centro de membrana en Z
        # Se promedia la coordenada Z de los atomos P (fosfato)
        # de POPC en los ultimos 50 ns (frames 3000-4000)
        # usando gmx traj. El resultado es MCENTER en Angstroms.
        # --------------------------------------------------
        echo "  [INFO] Midiendo centro de membrana (coordenada Z de fosfatos POPC)..."

        MCENTER_FILE="$ENERGIA_DIR/membrane_center.txt"
        if [ ! -f "$MCENTER_FILE" ]; then
            # Extraer z-media de atomos P de POPC del frame final
            # gmx traj devuelve coordenadas en nm -> convertir a Angstrom (*10)
            echo "POPC" | gmx traj                 -s "$TPR"                 -f "$XTC"                 -n "$NDX"                 -ox "$ENERGIA_DIR/popc_coords.xvg"                 -b 150000 -e 200000                 -com 2>/dev/null || true

            # Si no hay grupo POPC en el ndx, estimar desde step5_input.pdb
            if [ ! -f "$ENERGIA_DIR/popc_coords.xvg" ]; then
                echo "  [INFO] Estimando mcenter desde coordenadas PDB..."
                MCENTER=$(grep "^ATOM" "$SIM_DIR/step5_input.pdb" |                     awk '$3=="P" || $3==" P " {
                        # columna z en PDB = cols 47-54
                        z=substr($0,47,8)+0; sum+=z; n++
                    } END {
                        if(n>0) printf "%.2f", sum/n;
                        else print "0.00"
                    }')
            else
                # Parsear el xvg: tomar media de la columna z (col 4 si com, o promediar)
                MCENTER=$(awk '/^[^#@]/{sum+=$4; n++} END{if(n>0) printf "%.2f", sum/n*10; else print "0.00"}'                     "$ENERGIA_DIR/popc_coords.xvg")
                rm -f "$ENERGIA_DIR/popc_coords.xvg"
            fi

            echo "$MCENTER" > "$MCENTER_FILE"
            echo "  [INFO] mcenter medido: ${MCENTER} Angstrom (documentado — memopt no disponible en v1.6.4)"
        else
            MCENTER=$(cat "$MCENTER_FILE")
            echo "  [INFO] mcenter (cached): ${MCENTER} Angstrom (documentado)"
        fi

        # Archivo de input MM-GBSA con modelo de membrana implicita
        # mthick=26.0 : nucleo hidrofobico POPC (~26 A, sin cabezas fosfato)
        # emem=1.0    : constante dielectrica del nucleo lipidico (hidrofobico)
        # mcenter     : centro de la bicapa en Z, medido de la trayectoria
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

        echo "  [INFO] Archivo mmpbsa.in generado:"
        cat "$ENERGIA_DIR/mmpbsa.in"
        echo ""
        echo "  [INFO] Lanzando gmx_MMPBSA (serial + OMP_NUM_THREADS=20)..."
        echo "  [INFO] Output completo — visible en terminal y guardado en:"
        echo "         $LOG"
        echo ""

        cd "$ENERGIA_DIR"

        # 2>&1 | tee: output visible en terminal Y guardado en log
        echo "  [INFO] OMP_NUM_THREADS=$OMP_NUM_THREADS (todos los P-cores del i9-14900K)"
        echo ""
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

        # Verificar resultado
        if [ ! -f "$RESULTADO" ]; then
            echo ""
            echo "  [ERROR] No se genero FINAL_RESULTS_MMGBSA.dat"
            echo "  Revisa el log: $LOG"
            continue
        fi

        # Extraer y mostrar resultado DeltaG
        T_FIN=$(date +%s)
        T_ELAPSED=$(( T_FIN - T_INICIO ))
        T_MIN=$(( T_ELAPSED / 60 ))
        T_SEG=$(( T_ELAPSED % 60 ))

        echo ""
        echo "  ====== RESULTADO [$REPLICA / $PH] ======"
        grep "TOTAL" "$RESULTADO" | tail -1
        echo "  Tiempo: ${T_MIN}m ${T_SEG}s"
        echo "  ============================================"

    done
done

# ----------------------------------------------------------
# RESUMEN FINAL
# ----------------------------------------------------------
T_FIN_GLOBAL=$(date +%s)
T_TOTAL=$(( T_FIN_GLOBAL - INICIO_GLOBAL ))
T_MIN_TOTAL=$(( T_TOTAL / 60 ))
T_SEG_TOTAL=$(( T_TOTAL % 60 ))

echo ""
echo "=========================================================="
echo "  2_mmgbsa_energia.sh FINALIZADO"
echo "  Tiempo total: ${T_MIN_TOTAL}m ${T_SEG_TOTAL}s"
echo ""
echo "  DeltaG binding por condicion:"
echo ""
printf "  %-14s %-6s  %s\n" "REPLICA" "pH" "DeltaG (kcal/mol)"
echo "  $(printf '%0.s-' {1..55})"

for REPLICA in $REPLICAS; do
    for PH in $PHs; do
        R="$ANALISIS_DIR/$REPLICA/$PH/Energia/FINAL_RESULTS_MMGBSA.dat"
        if [ -f "$R" ]; then
            # Extrae media ± DE de la linea TOTAL del bloque Delta
            LINEA=$(grep "TOTAL" "$R" | tail -1 | \
                awk '{printf "%.2f +/- %.2f", $2, $4}')
            printf "  %-14s %-6s  %s kcal/mol\n" "$REPLICA" "$PH" "$LINEA"
        else
            printf "  %-14s %-6s  [NO CALCULADO]\n" "$REPLICA" "$PH"
        fi
    done
    echo ""
done

echo "  SIGUIENTE PASO:"
echo "  ./3_mmgbsa_decomp.sh
"
echo "=========================================================="
