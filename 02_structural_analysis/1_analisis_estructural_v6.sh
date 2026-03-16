#!/bin/bash
# ==========================================================
# EXTRACTOR ESTRUCTURAL V6 — CHARMM/GROMACS 2025
# Sistema: Semaglutida (PROA) + GLP-1R (PROB) + membrana POPC
# PROA: atomos 1-464   (Ligando / Semaglutida)
# PROB: atomos 465-6802 (Receptor / GLP-1R)
# 4001 frames = 200 ns por replica
# Replicas: 1, 2 y 3
# Pasos: indice, centrar, RMSD, RMSF, HBond, Rg
# ---------------------------------------------------------
# CAMBIOS vs V5:
#   - Replica_3 incluida
#   - Todo el output de GROMACS visible en terminal Y en log
#   - Radio de giro (Rg) agregado como Paso 6
#   - Temporizador por replica/pH
#   - Limpieza de archivos temporales en el directorio correcto
#   - set -e removido: si un paso falla, el pH se salta pero
#     el script sigue con el siguiente — no se cuelga
# ==========================================================

set -uo pipefail

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true

REPLICAS="Replica_1 Replica_2 Replica_3"
PHs="ph50 ph70 ph80"
BASE="$HOME/Documentos/Proyecto_Tesis_Renato_MD"
ANALISIS_DIR="$BASE/Analisis_Global"
INICIO_GLOBAL=$(date +%s)

echo ""
echo "=========================================================="
echo "  EXTRACTOR ESTRUCTURAL V6 — CHARMM/GROMACS 2025"
echo "  Replicas : $REPLICAS"
echo "  pHs      : $PHs"
echo "  Base     : $BASE"
echo "  Output   : $ANALISIS_DIR"
echo "  Inicio   : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=========================================================="

# ----------------------------------------------------------
# FUNCION AUXILIAR: ejecuta un comando GROMACS mostrando
# su output en pantalla Y guardandolo en el log
# Uso: run_gmx LOGFILE <comando gmx con argumentos>
# ----------------------------------------------------------
run_gmx() {
    local LOG="$1"
    shift
    echo "  [CMD] $*" | tee -a "$LOG"
    "$@" 2>&1 | tee -a "$LOG"
    return "${PIPESTATUS[0]}"
}

# ----------------------------------------------------------
# BUCLE PRINCIPAL
# ----------------------------------------------------------
for REPLICA in $REPLICAS; do
    for PH in $PHs; do

        T_INICIO=$(date +%s)
        SIM_DIR="$BASE/$REPLICA/$PH"
        OUT_DIR="$ANALISIS_DIR/$REPLICA/$PH"

        echo ""
        echo "=========================================================="
        echo "  >>> [$REPLICA / $PH]  $(date '+%H:%M:%S')"
        echo "=========================================================="

        # ---------- verificar archivos minimos ----------
        FALTAN=0
        for archivo in step5_input.pdb step7_1.tpr step7_1.xtc; do
            if [ ! -f "$SIM_DIR/$archivo" ]; then
                echo "  [SKIP] Falta $archivo en $SIM_DIR"
                FALTAN=1
            fi
        done
        [ $FALTAN -eq 1 ] && echo "  [SKIP] $REPLICA/$PH omitido." && continue

        mkdir -p "$OUT_DIR"/{RMSD,RMSF,HBond,Visuales,Rg}
        LOG="$OUT_DIR/pipeline_v6.log"
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
            echo "    PROA: 1-$LAST_PROA  PROB: $FIRST_PROB-$LAST_PROB" >> "$LOG"

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
                echo "  [ERROR] Indice sin grupos Ligando/Receptor. Grupos encontrados:"
                grep "^\[" "$NDX"
                echo "  Revisa: $NDX"
                echo "  Revisa: $LOG"
                continue
            fi
            echo "  [1/6] Indice OK."
        else
            echo ""
            echo "  [1/6] Indice existente — omitiendo."
            echo "    Grupos: $(grep '^\[' "$NDX" | tr '\n' ' ')"
        fi

        # --------------------------------------------------
        # PASO 2: CENTRAR TRAYECTORIA
        # -pbc mol : moleculas enteras (sin artefactos PBC)
        # Centro   : Protein (grupo 1 del ndx por defecto)
        # Salida   : todo el sistema (grupo 0)
        # --------------------------------------------------
        XTC_CENTRADO="$OUT_DIR/Visuales/centrado.xtc"
        if [ ! -f "$XTC_CENTRADO" ]; then
            echo ""
            echo "  [2/6] Centrando trayectoria (puede tardar varios minutos)..."
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
        # rmsd_receptor          : estabilidad GLP-1R vs GLP-1R
        # rmsd_ligando_vs_receptor: movimiento semaglutida relativo al receptor
        # rmsd_backbone          : backbone complejo
        # Referencia: t=0 (step5_input.pdb)
        # --------------------------------------------------
        echo ""
        echo "  [3/6] Calculando RMSD (3 archivos)..."

        if [ ! -f "$OUT_DIR/RMSD/rmsd_receptor.xvg" ]; then
            echo "    -> rmsd_receptor..."
            echo "Receptor
Receptor" | run_gmx "$LOG" gmx rms \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSD/rmsd_receptor.xvg" \
                -tu ns
        else
            echo "    -> rmsd_receptor.xvg existente."
        fi

        if [ ! -f "$OUT_DIR/RMSD/rmsd_ligando_vs_receptor.xvg" ]; then
            echo "    -> rmsd_ligando_vs_receptor..."
            echo "Receptor
Ligando" | run_gmx "$LOG" gmx rms \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSD/rmsd_ligando_vs_receptor.xvg" \
                -tu ns
        else
            echo "    -> rmsd_ligando_vs_receptor.xvg existente."
        fi

        if [ ! -f "$OUT_DIR/RMSD/rmsd_backbone.xvg" ]; then
            echo "    -> rmsd_backbone..."
            echo "Backbone
Backbone" | run_gmx "$LOG" gmx rms \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSD/rmsd_backbone.xvg" \
                -tu ns
        else
            echo "    -> rmsd_backbone.xvg existente."
        fi

        echo "  [3/6] RMSD OK."

        # --------------------------------------------------
        # PASO 4: RMSF
        # -res : promedia por residuo (no por atomo)
        # --------------------------------------------------
        echo ""
        echo "  [4/6] Calculando RMSF (2 archivos)..."

        if [ ! -f "$OUT_DIR/RMSF/rmsf_receptor.xvg" ]; then
            echo "    -> rmsf_receptor..."
            echo "Receptor" | run_gmx "$LOG" gmx rmsf \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSF/rmsf_receptor.xvg" \
                -res
        else
            echo "    -> rmsf_receptor.xvg existente."
        fi

        if [ ! -f "$OUT_DIR/RMSF/rmsf_ligando.xvg" ]; then
            echo "    -> rmsf_ligando..."
            echo "Ligando" | run_gmx "$LOG" gmx rmsf \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/RMSF/rmsf_ligando.xvg" \
                -res
        else
            echo "    -> rmsf_ligando.xvg existente."
        fi

        echo "  [4/6] RMSF OK."

        # --------------------------------------------------
        # PASO 5: PUENTES DE HIDROGENO
        # Usa TPR (no PDB) para topologia completa con cargas
        # -num : numero de hbonds por frame (grafico vs tiempo)
        # El hbond.ndx que genera gmx se borra al terminar
        # --------------------------------------------------
        echo ""
        echo "  [5/6] Calculando puentes de hidrogeno (puede tardar)..."

        if [ ! -f "$OUT_DIR/HBond/hbonds_Rec_Lig.xvg" ]; then
            echo "Receptor
Ligando" | run_gmx "$LOG" gmx hbond \
                -s "$SIM_DIR/step7_1.tpr" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -num "$OUT_DIR/HBond/hbonds_Rec_Lig.xvg"
            # limpiar backups que gromacs genera en cwd
            rm -f \#hbond* hbond.ndx* 2>/dev/null || true
        else
            echo "    -> hbonds_Rec_Lig.xvg existente."
        fi

        echo "  [5/6] HBond OK."

        # --------------------------------------------------
        # PASO 6: RADIO DE GIRO (Rg)
        # Mide compacidad global del complejo
        # -o : Rg total + componentes x, y, z
        # --------------------------------------------------
        echo ""
        echo "  [6/6] Calculando radio de giro (Rg)..."

        if [ ! -f "$OUT_DIR/Rg/rg_receptor.xvg" ]; then
            echo "    -> rg_receptor..."
            echo "Receptor" | run_gmx "$LOG" gmx gyrate \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/Rg/rg_receptor.xvg"
        else
            echo "    -> rg_receptor.xvg existente."
        fi

        if [ ! -f "$OUT_DIR/Rg/rg_complejo.xvg" ]; then
            echo "    -> rg_complejo (Protein = Ligando+Receptor)..."
            # Grupo 1 "Protein" = PROA(semaglutida) + PROB(GLP-1R) = complejo completo
            # Verificado: 463 + 6330 = 6793 atomos = Protein
            echo "Protein" | run_gmx "$LOG" gmx gyrate \
                -s "$SIM_DIR/step5_input.pdb" \
                -f "$XTC_CENTRADO" \
                -n "$NDX" \
                -o "$OUT_DIR/Rg/rg_complejo.xvg"
        else
            echo "    -> rg_complejo.xvg existente."
        fi

        echo "  [6/6] Rg OK."

        # ---------- resumen de tiempo ----------
        T_FIN=$(date +%s)
        T_ELAPSED=$(( T_FIN - T_INICIO ))
        T_MIN=$(( T_ELAPSED / 60 ))
        T_SEG=$(( T_ELAPSED % 60 ))

        XVG_COUNT=$(find "$OUT_DIR" -name "*.xvg" 2>/dev/null | wc -l)

        echo ""
        echo "  ✓ [$REPLICA / $PH] completado en ${T_MIN}m ${T_SEG}s"
        echo "    Archivos .xvg generados: $XVG_COUNT"
        echo "    Log: $LOG"
        echo "----------------------------------------------------------"

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
echo "  PIPELINE V6 FINALIZADO"
echo "  Tiempo total: ${T_MIN_TOTAL}m ${T_SEG_TOTAL}s"
echo ""
echo "  Archivos .xvg por replica/pH:"
for REPLICA in $REPLICAS; do
    for PH in $PHs; do
        DIR="$ANALISIS_DIR/$REPLICA/$PH"
        XVG_COUNT=$(find "$DIR" -name "*.xvg" 2>/dev/null | wc -l)
        if [ "$XVG_COUNT" -ge 8 ]; then
            ESTADO="OK ($XVG_COUNT xvg)"
        elif [ "$XVG_COUNT" -gt 0 ]; then
            ESTADO="INCOMPLETO ($XVG_COUNT xvg)"
        else
            ESTADO="SIN DATOS"
        fi
        echo "    $REPLICA/$PH  —  $ESTADO"
    done
done
echo ""
echo "  SIGUIENTE PASO:"
echo "  ./3_mmgbsa_energia.sh"
echo "=========================================================="
