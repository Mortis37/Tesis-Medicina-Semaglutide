#!/bin/bash
# ==========================================================
# 3_mmgbsa_decomp.sh
# DESCOMPOSICION POR RESIDUO — MM-GBSA (gmx_MMPBSA)
# Sistema: Semaglutida (PROA) + GLP-1R (PROB) + POPC
#
# Prerequisito: 2_mmgbsa_energia.sh completado
#   (necesita: centrado.xtc + index_analisis.ndx + topol.top)
#
# idecomp=2 : descomposicion por pares (residuo-residuo)
#   Calcula la contribucion energetica de cada residuo del
#   receptor a la union con semaglutida.
#   Permite identificar hotspots farmacologicos (DeltaG < -1.5 kcal/mol)
#
# ADVERTENCIA: idecomp=2 es 10-20x mas lento que energia sola.
#   Tiempo estimado: ~3-5 horas por condicion (RTX 4090)
#   Total para 9 condiciones: ~1-2 dias
#
# FRAME_INICIO=3000 (150-200 ns) — justificado por 1.2_verificar_convergencia.sh
# igb=5 (GB Onufriev) — identico a 2_mmgbsa_energia.sh
# memopt no disponible en gmx_MMPBSA 1.6.4
# mcenter reutilizado del script 2 (documentado en log)
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
echo "  3_mmgbsa_decomp.sh — DESCOMPOSICION POR RESIDUO"
echo "  Replicas      : $REPLICAS"
echo "  Frames        : $FRAME_INICIO a $FRAME_FIN (cada $INTERVALO)"
echo "  Frames totales: $N_FRAMES por condicion (= 150-200 ns)"
echo "  idecomp=2  igb=5  saltcon=0.150 M"
echo "  Inicio        : $(date '+%Y-%m-%d %H:%M:%S')"
echo "  ADVERTENCIA   : Este calculo tarda varias horas por condicion."
echo "=========================================================="

for REPLICA in $REPLICAS; do
    for PH in $PHs; do

        T_INICIO=$(date +%s)
        SIM_DIR="$BASE/$REPLICA/$PH"
        OUT_DIR="$ANALISIS_DIR/$REPLICA/$PH"
        DECOMP_DIR="$OUT_DIR/Decomp"
        LOG="$DECOMP_DIR/mmgbsa_decomp.log"

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

        # Si ya existe resultado, mostrar resumen y saltar
        RESULTADO="$DECOMP_DIR/FINAL_RESULTS_DECOMP.dat"
        if [ -f "$RESULTADO" ]; then
            echo "  [OK] Ya calculado anteriormente."
            echo "  Top 5 residuos por contribucion energetica:"
            grep -A 999 "Residue Energy Decomposition" "$RESULTADO" 2>/dev/null | \
                grep "^[A-Z]" | sort -k2 -n | head -5 || \
                echo "  (usa Python para parsear el .dat completo)"
            continue
        fi

        mkdir -p "$DECOMP_DIR"

        # --------------------------------------------------
        # Reutilizar mcenter calculado por 2_mmgbsa_energia.sh
        # Si no existe, calcularlo igual que en el script 2
        # --------------------------------------------------
        MCENTER_FILE="$OUT_DIR/Energia/membrane_center.txt"
        if [ -f "$MCENTER_FILE" ]; then
            MCENTER=$(cat "$MCENTER_FILE")
            echo "  [INFO] mcenter (desde script 2): ${MCENTER} Angstrom"
        else
            echo "  [INFO] Midiendo mcenter desde PDB..."
            MCENTER=$(grep "^ATOM" "$SIM_DIR/step5_input.pdb" |                 awk '$3=="P" || $3==" P " {
                    z=substr($0,47,8)+0; sum+=z; n++
                } END {
                    if(n>0) printf "%.2f", sum/n;
                    else print "0.00"
                }')
            echo "  [INFO] mcenter medido: ${MCENTER} Angstrom (documentado — memopt no disponible en v1.6.4)"
        fi

        # --------------------------------------------------
        # PASO PREVIO: identificar residuos en la interfaz
        # Criterio: residuos del receptor a <= 6 Å del ligando
        # Se usa gmx select sobre el frame promedio (t=175 ns)
        # El resultado se guarda como lista para print_res
        # --------------------------------------------------
        RESLIST_FILE="$DECOMP_DIR/interface_residues.txt"
        if [ ! -f "$RESLIST_FILE" ]; then
            echo "  [INFO] Identificando residuos en interfaz Receptor-Ligando (<=6 A)..."

            # Extraer un frame representativo (~175 ns = frame 3500)
            FRAME_PDB="$DECOMP_DIR/frame_ref.pdb"
            echo "System" | gmx trjconv                 -s "$SIM_DIR/step7_1.tpr"                 -f "$OUT_DIR/Visuales/centrado.xtc"                 -o "$FRAME_PDB"                 -dump 175000 2>/dev/null || true

            if [ -f "$FRAME_PDB" ]; then
                # gmx select: residuos del receptor dentro de 6A del ligando
                # Salida: lista de numeros de residuo (resid)
                gmx select                     -s "$FRAME_PDB"                     -n "$NDX"                     -select "same residue as (group Receptor and within 0.6 of group Ligando)"                     -on "$DECOMP_DIR/interface.ndx"                     -os "$DECOMP_DIR/interface_sizes.xvg" 2>/dev/null || true

                if [ -f "$DECOMP_DIR/interface.ndx" ]; then
                    # Extraer numeros de residuo desde el PDB
                    # Cruzar atomos del ndx con resid del PDB
                    ATOM_LIST=$(grep -A 9999 "\[ Receptor_interface \]\|\[ selection_0 \]"                         "$DECOMP_DIR/interface.ndx" 2>/dev/null |                         grep -v "^\[" | tr ' ' '
' | grep -v "^$" | head -500)

                    if [ -n "$ATOM_LIST" ]; then
                        # Convertir numeros de atomo a residuos unicos del PDB
                        PRINT_RES=$(awk -v atoms="$ATOM_LIST" '
                            BEGIN {
                                split(atoms, arr, "
")
                                for(i in arr) atomset[arr[i]] = 1
                            }
                            /^ATOM/ {
                                atomnum = substr($0,7,5)+0
                                resnum  = substr($0,23,4)+0
                                if(atomnum in atomset) residues[resnum] = 1
                            }
                            END {
                                out = ""
                                for(r in residues) {
                                    out = (out == "") ? r : out "," r
                                }
                                print out
                            }
                        ' "$FRAME_PDB")
                    fi
                fi
            fi

            # Fallback: si falla gmx select, usar todos los residuos del receptor
            if [ -z "${PRINT_RES:-}" ]; then
                echo "  [WARN] No se pudo calcular interfaz — usando todos los residuos del receptor"
                # Rango de residuos PROB desde el PDB
                PRINT_RES=$(grep "^ATOM" "$SIM_DIR/step5_input.pdb" |                     awk '{seg=substr($0,73,4); gsub(/ /,"",seg);
                          if(seg=="PROB") res=substr($0,23,4)+0}
                         END{print "1-" res}')
            fi

            echo "$PRINT_RES" > "$RESLIST_FILE"
            echo "  [INFO] Residuos interfaz: $PRINT_RES"
            rm -f "$DECOMP_DIR/interface.ndx" "$DECOMP_DIR/interface_sizes.xvg"                   "$DECOMP_DIR/frame_ref.pdb" 2>/dev/null || true
        else
            PRINT_RES=$(cat "$RESLIST_FILE")
            echo "  [INFO] Residuos interfaz (cached): $PRINT_RES"
        fi

        # --------------------------------------------------
        # Input MM-GBSA con descomposicion por pares + membrana
        # idecomp=2 en &general (requerido por AmberTools)
        # print_res: solo residuos en la interfaz (<= 6 A)
        #   reduce tiempo de calculo 10-20x vs todos los residuos
        # --------------------------------------------------
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

        echo "  [INFO] Archivo mmpbsa_decomp.in generado:"
        cat "$DECOMP_DIR/mmpbsa_decomp.in"
        echo ""
        echo "  [INFO] Lanzando gmx_MMPBSA decomp (OMP_NUM_THREADS=20)..."
        echo "  [INFO] Output completo guardado en:"
        echo "         $LOG"
        echo ""

        cd "$DECOMP_DIR"

        echo "  [INFO] OMP_NUM_THREADS=$OMP_NUM_THREADS"
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

        # Verificar resultado
        if [ ! -f "$RESULTADO" ]; then
            echo ""
            echo "  [ERROR] No se genero FINAL_RESULTS_DECOMP.dat"
            echo "  Revisa el log: $LOG"
            continue
        fi

        T_FIN=$(date +%s)
        T_ELAPSED=$(( T_FIN - T_INICIO ))
        T_MIN=$(( T_ELAPSED / 60 ))
        T_SEG=$(( T_ELAPSED % 60 ))

        echo ""
        echo "  ====== COMPLETADO [$REPLICA / $PH] ======"
        echo "  Resultado: $RESULTADO"
        echo "  CSV:       $DECOMP_DIR/FINAL_DECOMP.csv"
        echo "  Tiempo:    ${T_MIN}m ${T_SEG}s"
        echo "  ============================================"

    done
done

# ----------------------------------------------------------
# RESUMEN FINAL — listar los CSV generados
# ----------------------------------------------------------
T_FIN_GLOBAL=$(date +%s)
T_TOTAL=$(( T_FIN_GLOBAL - INICIO_GLOBAL ))
T_MIN_TOTAL=$(( T_TOTAL / 60 ))
T_SEG_TOTAL=$(( T_TOTAL % 60 ))

echo ""
echo "=========================================================="
echo "  3_mmgbsa_decomp.sh FINALIZADO"
echo "  Tiempo total: ${T_MIN_TOTAL}m ${T_SEG_TOTAL}s"
echo ""
echo "  Archivos generados por condicion:"
echo ""
printf "  %-14s %-6s  %-6s  %s\n" "REPLICA" "pH" "DAT" "CSV"
echo "  $(printf '%0.s-' {1..65})"

for REPLICA in $REPLICAS; do
    for PH in $PHs; do
        D="$ANALISIS_DIR/$REPLICA/$PH/Decomp"
        DAT="N"; CSV="N"
        [ -f "$D/FINAL_RESULTS_DECOMP.dat" ] && DAT="OK"
        [ -f "$D/FINAL_DECOMP.csv"         ] && CSV="OK"
        printf "  %-14s %-6s  %-6s  %s\n" "$REPLICA" "$PH" "$DAT" "$CSV"
    done
    echo ""
done

echo "  SIGUIENTE PASO:"
echo "  Usar Python/pandas para parsear los CSV y extraer hotspots:"
echo "  (residuos con DeltaG < -1.5 kcal/mol = contribucion favorable)"
echo "=========================================================="
