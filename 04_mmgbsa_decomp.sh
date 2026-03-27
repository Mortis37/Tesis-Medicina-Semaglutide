#!/bin/bash
# ==========================================================
# 04_mmgbsa_decomp.sh
# DESCOMPOSICION POR RESIDUO — MM-GBSA (gmx_MMPBSA)
# Sistema: Semaglutida (PROA) + GLP-1R (PROB) + POPC
#
# Prerequisito: 03_mmgbsa_energia.sh completado
#
# idecomp=2 : descomposicion por pares (residuo-residuo)
# Identifica hotspots farmacologicos (DeltaG < -1.5 kcal/mol)
#
# Tiempo estimado: ~3-5 horas por condicion
#
# Uso: bash 04_mmgbsa_decomp.sh
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
echo "  04_mmgbsa_decomp.sh — DESCOMPOSICION POR RESIDUO"
echo "  Frames : $FRAME_INICIO-$FRAME_FIN (cada $INTERVALO = $N_FRAMES frames)"
echo "  idecomp=2  igb=5  saltcon=0.150 M"
echo "  Inicio : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================================="

for REPLICA in $REPLICAS; do
    for PH in $PHs; do

        T_INICIO=$(date +%s)
        SIM_DIR="$BASE/$REPLICA/$PH"
        OUT_DIR="$ANALISIS_DIR/$REPLICA/$PH"
        DECOMP_DIR="$OUT_DIR/Decomp"
        LOG="$DECOMP_DIR/mmgbsa_decomp.log"

        echo ""
        echo "  >>> [$REPLICA / $PH]  $(date '+%H:%M:%S')"

        NDX="$OUT_DIR/index_analisis.ndx"
        XTC="$OUT_DIR/Visuales/centrado.xtc"
        TPR="$SIM_DIR/step7_1.tpr"
        TOP="$SIM_DIR/topol.top"

        FALTAN=0
        for f in "$NDX" "$XTC" "$TPR" "$TOP"; do
            [ ! -f "$f" ] && echo "  [SKIP] Falta: $f" && FALTAN=1
        done
        [ $FALTAN -eq 1 ] && continue

        RESULTADO="$DECOMP_DIR/FINAL_RESULTS_DECOMP.dat"
        if [ -f "$RESULTADO" ]; then
            echo "  [OK] Ya calculado."
            continue
        fi

        mkdir -p "$DECOMP_DIR"

        # Reutilizar mcenter de script 03
        MCENTER_FILE="$OUT_DIR/Energia/membrane_center.txt"
        if [ -f "$MCENTER_FILE" ]; then
            MCENTER=$(cat "$MCENTER_FILE")
        else
            MCENTER=$(grep "^ATOM" "$SIM_DIR/step5_input.pdb" | \
                awk '$3=="P"||$3==" P "{z=substr($0,47,8)+0;sum+=z;n++}
                     END{if(n>0)printf "%.2f",sum/n;else print "0.00"}')
        fi
        echo "  [INFO] mcenter: ${MCENTER} A (documentado)"

        # Identificar residuos en la interfaz (<= 6 A)
        RESLIST_FILE="$DECOMP_DIR/interface_residues.txt"
        if [ ! -f "$RESLIST_FILE" ]; then
            echo "  [INFO] Identificando interfaz Receptor-Ligando (<=6 A)..."

            FRAME_PDB="$DECOMP_DIR/frame_ref.pdb"
            echo "System" | gmx trjconv \
                -s "$TPR" -f "$XTC" \
                -o "$FRAME_PDB" -dump 175000 2>/dev/null || true

            PRINT_RES=""
            if [ -f "$FRAME_PDB" ]; then
                gmx select \
                    -s "$FRAME_PDB" -n "$NDX" \
                    -select "same residue as (group Receptor and within 0.6 of group Ligando)" \
                    -on "$DECOMP_DIR/interface.ndx" \
                    -os "$DECOMP_DIR/interface_sizes.xvg" 2>/dev/null || true

                if [ -f "$DECOMP_DIR/interface.ndx" ]; then
                    ATOM_LIST=$(grep -A 9999 "\[ Receptor_interface \]\|\[ selection_0 \]" \
                        "$DECOMP_DIR/interface.ndx" 2>/dev/null | \
                        grep -v "^\[" | tr ' ' '\n' | grep -v "^$" | head -500)

                    if [ -n "$ATOM_LIST" ]; then
                        PRINT_RES=$(awk -v atoms="$ATOM_LIST" '
                            BEGIN{split(atoms,arr,"\n");for(i in arr)atomset[arr[i]]=1}
                            /^ATOM/{atomnum=substr($0,7,5)+0;resnum=substr($0,23,4)+0
                                    if(atomnum in atomset)residues[resnum]=1}
                            END{out="";for(r in residues)out=(out=="")?r:out","r;print out}
                        ' "$FRAME_PDB")
                    fi
                fi
                rm -f "$FRAME_PDB" "$DECOMP_DIR/interface.ndx" \
                      "$DECOMP_DIR/interface_sizes.xvg" 2>/dev/null || true
            fi

            if [ -z "${PRINT_RES:-}" ]; then
                echo "  [WARN] Fallback: usando todos los residuos del receptor"
                PRINT_RES=$(grep "^ATOM" "$SIM_DIR/step5_input.pdb" | \
                    awk '{seg=substr($0,73,4);gsub(/ /,"",seg);
                          if(seg=="PROB")res=substr($0,23,4)+0}
                         END{print "1-"res}')
            fi

            echo "$PRINT_RES" > "$RESLIST_FILE"
            echo "  [INFO] Residuos interfaz: $PRINT_RES"
        else
            PRINT_RES=$(cat "$RESLIST_FILE")
        fi

        # Input descomposicion
        cat > "$DECOMP_DIR/mmpbsa_decomp.in" << EOF
&general
  sys_name="${REPLICA}_${PH}_decomp",
  startframe=$FRAME_INICIO,
  endframe=$FRAME_FIN,
  interval=$INTERVALO,
  verbose=2,
/
&gb
  igb=5,
  saltcon=0.150,
/
&decomp
  idecomp=2,
  csv_format=1,
  dec_verbose=3,
  print_res="A/7-36 B/${PRINT_RES}",
/
EOF

        echo "  [INFO] Lanzando gmx_MMPBSA decomp (OMP=$OMP_NUM_THREADS)..."
        cd "$DECOMP_DIR"

        gmx_MMPBSA -O \
            -i  mmpbsa_decomp.in \
            -cs "$TPR" \
            -ci "$NDX" \
            -cg Receptor Ligando \
            -ct "$XTC" \
            -cp "$TOP" \
            -o  FINAL_RESULTS_DECOMP.dat \
            -eo FINAL_DECOMP.csv \
            -nogui 2>&1 | tee -a "$LOG"

        cd "$BASE"

        if [ ! -f "$RESULTADO" ]; then
            echo "  [ERROR] No se genero resultado — ver $LOG"
            continue
        fi

        T_MIN=$(( ($(date +%s) - T_INICIO) / 60 ))
        echo "  [$REPLICA / $PH] completado en ${T_MIN}m"

    done
done

# Resumen
echo ""
echo "=========================================================="
echo "  04_mmgbsa_decomp.sh FINALIZADO"
echo "  Tiempo total: $(( ($(date +%s) - INICIO_GLOBAL) / 60 ))m"
echo ""
printf "  %-14s %-6s  %-6s  %s\n" "REPLICA" "pH" "DAT" "CSV"
echo "  $(printf '%0.s-' {1..45})"
for REPLICA in $REPLICAS; do
    for PH in $PHs; do
        D="$ANALISIS_DIR/$REPLICA/$PH/Decomp"
        DAT="N"; CSV="N"
        [ -f "$D/FINAL_RESULTS_DECOMP.dat" ] && DAT="OK"
        [ -f "$D/FINAL_DECOMP.csv"         ] && CSV="OK"
        printf "  %-14s %-6s  %-6s  %s\n" "$REPLICA" "$PH" "$DAT" "$CSV"
    done
done
echo ""
echo "  SIGUIENTE PASO: python3 05_analisis_hotspots.py"
echo "=========================================================="
