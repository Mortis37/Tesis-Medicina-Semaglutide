#!/bin/bash
# ==========================================================
# 02_analisis_estructural.sh
# Analisis estructural completo — CHARMM/GROMACS
# Sistema: Semaglutida (PROA) + GLP-1R (PROB) + membrana POPC
# PROA: atomos 1-464   (Ligando / Semaglutida)
# PROB: atomos 465-6802 (Receptor / GLP-1R)
# 4001 frames = 200 ns por replica
#
# Pasos:
#   1. Indice por segmento CHARMM (Ligando / Receptor)
#   2. Centrar trayectoria (-pbc mol)
#   3. RMSD (receptor, ligando vs receptor, backbone)
#   4. RMSF (receptor, ligando)
#   5. Puentes de hidrogeno receptor-ligando
#   6. Radio de giro (receptor, complejo)
#   7. Verificacion de convergencia (ventanas A/B/C)
#   8. RMSD vs estructura cristalografica (PDB 7KI0)
#
# Output: Analisis_Global/Replica_X/phYY/{RMSD,RMSF,HBond,Rg,Visuales}/
# Ningun archivo se genera dentro de las carpetas de replica originales.
#
# Uso: bash 02_analisis_estructural.sh
# ==========================================================

set -uo pipefail
export LC_NUMERIC=C

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true

BASE="$(cd "$(dirname "$0")" && pwd)"
REPLICAS="Replica_1 Replica_2 Replica_3"
PHs="ph50 ph70 ph80"
ANALISIS_DIR="$BASE/Analisis_Global"
GRAFICOS_DIR="$ANALISIS_DIR/Graficos"
INICIO_GLOBAL=$(date +%s)

mkdir -p "$GRAFICOS_DIR"

echo ""
echo "=========================================================="
echo "  02_analisis_estructural.sh"
echo "  Replicas : $REPLICAS"
echo "  pHs      : $PHs"
echo "  Output   : $ANALISIS_DIR"
echo "  Inicio   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================================="

# ----------------------------------------------------------
# FUNCION AUXILIAR: ejecuta un comando GROMACS mostrando
# su output en pantalla Y guardandolo en el log
# ----------------------------------------------------------
run_gmx() {
    local LOG="$1"; shift
    echo "  [CMD] $*" | tee -a "$LOG"
    "$@" 2>&1 | tee -a "$LOG"
    return "${PIPESTATUS[0]}"
}

# ==========================================================
# PARTE A: ANALISIS POR REPLICA/pH
# ==========================================================
for REPLICA in $REPLICAS; do
    for PH in $PHs; do

        T_INICIO=$(date +%s)
        SIM_DIR="$BASE/$REPLICA/$PH"
        OUT_DIR="$ANALISIS_DIR/$REPLICA/$PH"

        echo ""
        echo "=========================================================="
        echo "  >>> [$REPLICA / $PH]  $(date '+%H:%M:%S')"
        echo "=========================================================="

        # Verificar archivos minimos
        FALTAN=0
        for archivo in step5_input.pdb step7_1.tpr step7_1.xtc; do
            if [ ! -f "$SIM_DIR/$archivo" ]; then
                echo "  [SKIP] Falta $archivo en $SIM_DIR"
                FALTAN=1
            fi
        done
        [ $FALTAN -eq 1 ] && echo "  [SKIP] $REPLICA/$PH omitido." && continue

        mkdir -p "$OUT_DIR"/{RMSD,RMSF,HBond,Visuales,Rg}
        LOG="$OUT_DIR/pipeline.log"
        echo "" >> "$LOG"
        echo "=== $REPLICA/$PH  $(date '+%Y-%m-%d %H:%M:%S') ===" >> "$LOG"

        # --------------------------------------------------
        # PASO 1: INDICE POR SEGMENTO CHARMM
        # --------------------------------------------------
        NDX="$OUT_DIR/index_analisis.ndx"
        if [ ! -f "$NDX" ]; then
            echo ""
            echo "  [1/6] Creando indice por segmento CHARMM..."

            LAST_PROA=$(grep "^ATOM" "$SIM_DIR/step5_input.pdb" | \
                awk '{seg=substr($0,73,4); gsub(/ /,"",seg);
                      if(seg=="PROA") last=substr($0,7,5)+0}
                     END{print last}')

            FIRST_PROB=$(grep "^ATOM" "$SIM_DIR/step5_input.pdb" | \
                awk '{seg=substr($0,73,4); gsub(/ /,"",seg);
                      if(seg=="PROB"){n=substr($0,7,5)+0;
                      if(first==""||n<first) first=n}}
                     END{print first}')

            LAST_PROB=$(grep "^ATOM" "$SIM_DIR/step5_input.pdb" | \
                awk '{seg=substr($0,73,4); gsub(/ /,"",seg);
                      if(seg=="PROB") last=substr($0,7,5)+0}
                     END{print last}')

            echo "    PROA: atomos 1 a $LAST_PROA  (Ligando / Semaglutida)"
            echo "    PROB: atomos $FIRST_PROB a $LAST_PROB  (Receptor / GLP-1R)"

            run_gmx "$LOG" gmx make_ndx \
                -f "$SIM_DIR/step5_input.pdb" \
                -o "$NDX" << EOF
a 1-${LAST_PROA}
a ${FIRST_PROB}-${LAST_PROB}
q
EOF

            sed -i "s/\[ a_1-${LAST_PROA} \]/\[ Ligando \]/g"             "$NDX"
            sed -i "s/\[ a_${FIRST_PROB}-${LAST_PROB} \]/\[ Receptor \]/g" "$NDX"

            if ! grep -q "\[ Ligando \]" "$NDX" || ! grep -q "\[ Receptor \]" "$NDX"; then
                echo "  [ERROR] Indice sin grupos Ligando/Receptor."
                grep "^\[" "$NDX"
                continue
            fi
            echo "  [1/6] Indice OK."
        else
            echo ""
            echo "  [1/6] Indice existente — omitiendo."
        fi

        # --------------------------------------------------
        # PASO 2: CENTRAR TRAYECTORIA
        # --------------------------------------------------
        XTC_CENTRADO="$OUT_DIR/Visuales/centrado.xtc"
        if [ ! -f "$XTC_CENTRADO" ]; then
            echo ""
            echo "  [2/6] Centrando trayectoria..."
            echo "1
0" | run_gmx "$LOG" gmx trjconv \
                -s "$SIM_DIR/step7_1.tpr" \
                -f "$SIM_DIR/step7_1.xtc" \
                -o "$XTC_CENTRADO" \
                -pbc mol -center
            echo "  [2/6] Centrado OK."
        else
            echo ""
            echo "  [2/6] Trayectoria centrada existente — omitiendo."
        fi

        # --------------------------------------------------
        # PASO 3: RMSD
        # --------------------------------------------------
        echo ""
        echo "  [3/6] Calculando RMSD..."

        if [ ! -f "$OUT_DIR/RMSD/rmsd_receptor.xvg" ]; then
            echo "Receptor
Receptor" | run_gmx "$LOG" gmx rms \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSD/rmsd_receptor.xvg" -tu ns
        fi

        if [ ! -f "$OUT_DIR/RMSD/rmsd_ligando_vs_receptor.xvg" ]; then
            echo "Receptor
Ligando" | run_gmx "$LOG" gmx rms \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSD/rmsd_ligando_vs_receptor.xvg" -tu ns
        fi

        if [ ! -f "$OUT_DIR/RMSD/rmsd_backbone.xvg" ]; then
            echo "Backbone
Backbone" | run_gmx "$LOG" gmx rms \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSD/rmsd_backbone.xvg" -tu ns
        fi
        echo "  [3/6] RMSD OK."

        # --------------------------------------------------
        # PASO 4: RMSF (-res promedia por residuo)
        # --------------------------------------------------
        echo ""
        echo "  [4/6] Calculando RMSF..."

        if [ ! -f "$OUT_DIR/RMSF/rmsf_receptor.xvg" ]; then
            echo "Receptor" | run_gmx "$LOG" gmx rmsf \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSF/rmsf_receptor.xvg" -res
        fi

        if [ ! -f "$OUT_DIR/RMSF/rmsf_ligando.xvg" ]; then
            echo "Ligando" | run_gmx "$LOG" gmx rmsf \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSF/rmsf_ligando.xvg" -res
        fi
        echo "  [4/6] RMSF OK."

        # --------------------------------------------------
        # PASO 5: PUENTES DE HIDROGENO
        # --------------------------------------------------
        echo ""
        echo "  [5/6] Calculando puentes de hidrogeno..."

        if [ ! -f "$OUT_DIR/HBond/hbonds_Rec_Lig.xvg" ]; then
            echo "Receptor
Ligando" | run_gmx "$LOG" gmx hbond \
                -s "$SIM_DIR/step7_1.tpr" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -num "$OUT_DIR/HBond/hbonds_Rec_Lig.xvg"
            rm -f \#hbond* hbond.ndx* 2>/dev/null || true
        fi
        echo "  [5/6] HBond OK."

        # --------------------------------------------------
        # PASO 6: RADIO DE GIRO
        # --------------------------------------------------
        echo ""
        echo "  [6/6] Calculando radio de giro..."

        if [ ! -f "$OUT_DIR/Rg/rg_receptor.xvg" ]; then
            echo "Receptor" | run_gmx "$LOG" gmx gyrate \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/Rg/rg_receptor.xvg"
        fi

        if [ ! -f "$OUT_DIR/Rg/rg_complejo.xvg" ]; then
            echo "Protein" | run_gmx "$LOG" gmx gyrate \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/Rg/rg_complejo.xvg"
        fi
        echo "  [6/6] Rg OK."

        # Resumen de tiempo
        T_FIN=$(date +%s)
        T_ELAPSED=$(( T_FIN - T_INICIO ))
        T_MIN=$(( T_ELAPSED / 60 ))
        T_SEG=$(( T_ELAPSED % 60 ))
        XVG_COUNT=$(find "$OUT_DIR" -name "*.xvg" 2>/dev/null | wc -l)

        echo ""
        echo "  [$REPLICA / $PH] completado en ${T_MIN}m ${T_SEG}s — $XVG_COUNT archivos .xvg"
        echo "----------------------------------------------------------"

    done
done


# ==========================================================
# PARTE B: VERIFICACION DE CONVERGENCIA
# ==========================================================
echo ""
echo "=========================================================="
echo "  VERIFICACION DE CONVERGENCIA — RMSD RECEPTOR"
echo "=========================================================="

# Ventanas en nanosegundos
V_A_INI=0;   V_A_FIN=100
V_B_INI=100; V_B_FIN=200
V_C_INI=150; V_C_FIN=200

calcular_stats() {
    local FILE="$1" T_INI="$2" T_FIN="$3"
    awk -v ini="$T_INI" -v fin="$T_FIN" '
        /^[#@]/ { next }
        $1 >= ini && $1 <= fin { suma+=$2; suma2+=$2*$2; n++ }
        END {
            if (n==0) { printf "N/A" }
            else {
                media=suma/n; var=(suma2/n)-(media*media)
                de=(var>0)?sqrt(var):0
                printf "%.4f +/- %.4f nm (n=%d)", media, de, n
            }
        }
    ' "$FILE"
}

CONVERGE_OK=0
CONVERGE_WARN=0

printf "\n  %-14s %-6s  %-28s %-28s %-28s\n" \
    "REPLICA" "pH" "A: 0-100 ns" "B: 100-200 ns" "C: 150-200 ns"
echo "  $(printf '%0.s-' {1..104})"

for REPLICA in $REPLICAS; do
    for PH in $PHs; do
        XVG="$ANALISIS_DIR/$REPLICA/$PH/RMSD/rmsd_receptor.xvg"
        [ ! -f "$XVG" ] && continue

        STATS_A=$(calcular_stats "$XVG" $V_A_INI $V_A_FIN)
        STATS_B=$(calcular_stats "$XVG" $V_B_INI $V_B_FIN)
        STATS_C=$(calcular_stats "$XVG" $V_C_INI $V_C_FIN)

        printf "  %-14s %-6s  %-28s %-28s %-28s\n" \
            "$REPLICA" "$PH" "$STATS_A" "$STATS_B" "$STATS_C"
    done
    echo ""
done

# Criterio GPCR: DE_C < 0.10 nm, delta B->C < 0.05 nm
echo "  ANALISIS (criterio GPCR: DE_C < 0.10 nm, delta < 0.05 nm):"
echo ""

for REPLICA in $REPLICAS; do
    for PH in $PHs; do
        XVG="$ANALISIS_DIR/$REPLICA/$PH/RMSD/rmsd_receptor.xvg"
        [ ! -f "$XVG" ] && continue

        MEDIA_B=$(awk -v ini=$V_B_INI -v fin=$V_B_FIN '
            /^[#@]/{next} $1>=ini&&$1<=fin{s+=$2;n++}
            END{if(n>0) printf "%.4f",s/n; else print "0"}' "$XVG")

        MEDIA_C=$(awk -v ini=$V_C_INI -v fin=$V_C_FIN '
            /^[#@]/{next} $1>=ini&&$1<=fin{s+=$2;n++}
            END{if(n>0) printf "%.4f",s/n; else print "0"}' "$XVG")

        DE_C=$(awk -v ini=$V_C_INI -v fin=$V_C_FIN '
            /^[#@]/{next} $1>=ini&&$1<=fin{s+=$2;s2+=$2*$2;n++}
            END{if(n>0){m=s/n;v=(s2/n)-(m*m);printf "%.4f",(v>0)?sqrt(v):0}else print "0"}' "$XVG")

        DELTA=$(awk -v b="$MEDIA_B" -v c="$MEDIA_C" \
            'BEGIN{d=b-c;if(d<0)d=-d;printf "%.4f",d}')

        CONVERGE=$(awk -v de="$DE_C" -v d="$DELTA" \
            'BEGIN{if(de+0<0.10&&d+0<0.05)print "OK";else print "REVISAR"}')

        if [ "$CONVERGE" = "OK" ]; then
            ICONO="[OK]"; (( CONVERGE_OK++ )) || true
        else
            ICONO="[!!]"; (( CONVERGE_WARN++ )) || true
        fi

        printf "  %s %-14s %-6s  media_B=%.4f  media_C=%.4f  DE_C=%.4f  delta=%.4f -> %s\n" \
            "$ICONO" "$REPLICA" "$PH" "$MEDIA_B" "$MEDIA_C" "$DE_C" "$DELTA" "$CONVERGE"
    done
done

echo ""
echo "  Convergidos: $CONVERGE_OK  |  Revisar: $CONVERGE_WARN"
if [ "$CONVERGE_WARN" -eq 0 ]; then
    echo "  FRAME_INICIO=3000 (150 ns) confirmado para MM-GBSA."
else
    echo "  Revisar graficos antes de lanzar MM-GBSA."
fi


# ==========================================================
# PARTE C: RMSD vs ESTRUCTURA CRISTALOGRAFICA (PDB 7KI0)
# Valida que la pose de semaglutida se conserve tras equilibracion
# ==========================================================
echo ""
echo "=========================================================="
echo "  RMSD vs CRISTAL (PDB 7KI0) — Validacion de pose"
echo "=========================================================="

CSV_CRISTAL="$GRAFICOS_DIR/rmsd_cristal_data.csv"
echo "replica,ph_label,rmsd_angstrom" > "$CSV_CRISTAL"

for REPLICA in $REPLICAS; do
    REP_NUM=${REPLICA##*_}
    for PH in $PHs; do
        SIM_DIR="$BASE/$REPLICA/$PH"
        OUT_DIR="$ANALISIS_DIR/$REPLICA/$PH"
        NDX="$OUT_DIR/index_analisis.ndx"

        [ ! -f "$SIM_DIR/step7_1.tpr" ] && continue
        [ ! -f "$OUT_DIR/Visuales/centrado.xtc" ] && continue
        [ ! -f "$NDX" ] && continue

        case $PH in
            ph50) PH_LABEL="5.0" ;;
            ph70) PH_LABEL="7.4" ;;
            ph80) PH_LABEL="8.0" ;;
        esac

        TMP_FRAME="$OUT_DIR/Visuales/frame_t0.pdb"
        TMP_RMSD="$OUT_DIR/RMSD/rmsd_vs_cristal.xvg"

        if [ ! -f "$TMP_RMSD" ]; then
            # Extraer primer frame
            echo "System" | gmx trjconv \
                -s "$SIM_DIR/step7_1.tpr" \
                -f "$OUT_DIR/Visuales/centrado.xtc" \
                -o "$TMP_FRAME" \
                -dump 0 2>/dev/null || true

            if [ -f "$TMP_FRAME" ]; then
                # RMSD del backbone del ligando vs estructura de referencia
                echo "Ligando
Ligando" | gmx rms \
                    -s "$SIM_DIR/step5_input.pdb" \
                    -f "$TMP_FRAME" \
                    -n "$NDX" \
                    -o "$TMP_RMSD" 2>/dev/null || true
                rm -f "$TMP_FRAME"
            fi
        fi

        if [ -f "$TMP_RMSD" ]; then
            # Extraer RMSD en nm, convertir a Angstrom
            RMSD_NM=$(awk '/^[^#@]/{print $2; exit}' "$TMP_RMSD")
            RMSD_A=$(awk -v r="$RMSD_NM" 'BEGIN{printf "%.4f", r*10}')
            echo "$REP_NUM,$PH_LABEL,$RMSD_A" >> "$CSV_CRISTAL"
            printf "  %-12s %-6s  RMSD = %s A\n" "$REPLICA" "$PH" "$RMSD_A"
        fi
    done
done

echo ""
echo "  CSV generado: $CSV_CRISTAL"


# ==========================================================
# RESUMEN FINAL
# ==========================================================
T_FIN_GLOBAL=$(date +%s)
T_TOTAL=$(( T_FIN_GLOBAL - INICIO_GLOBAL ))
T_MIN_TOTAL=$(( T_TOTAL / 60 ))
T_SEG_TOTAL=$(( T_TOTAL % 60 ))

echo ""
echo "=========================================================="
echo "  02_analisis_estructural.sh FINALIZADO"
echo "  Tiempo total: ${T_MIN_TOTAL}m ${T_SEG_TOTAL}s"
echo ""
echo "  Archivos .xvg por replica/pH:"
for REPLICA in $REPLICAS; do
    for PH in $PHs; do
        DIR="$ANALISIS_DIR/$REPLICA/$PH"
        XVG_COUNT=$(find "$DIR" -name "*.xvg" 2>/dev/null | wc -l)
        if [ "$XVG_COUNT" -ge 8 ]; then ESTADO="OK ($XVG_COUNT xvg)"
        elif [ "$XVG_COUNT" -gt 0 ]; then ESTADO="INCOMPLETO ($XVG_COUNT xvg)"
        else ESTADO="SIN DATOS"; fi
        echo "    $REPLICA/$PH  —  $ESTADO"
    done
done
echo ""
echo "  SIGUIENTE PASO: bash 03_mmgbsa_energia.sh"
echo "=========================================================="
