#!/bin/bash
set -euo pipefail

# ==========================================================
# PIPELINE MD - REPLICA_3
# Hardware: Intel i9-14900K (32 hilos) + NVIDIA RTX 4090
# GROMACS optimizado para GPU + CPU mixto
#
# Optimizaciones vs script anterior:
#   -ntmpi 1        : 1 proceso MPI (obligatorio con 1 GPU)
#   -ntomp 16       : 16 hilos OpenMP para calculos CPU
#                     (deja 16 hilos libres para el sistema)
#   -gpu_id 0       : RTX 4090 exclusiva
#   -nb gpu         : non-bonded en GPU (mayor aceleracion)
#   -pme gpu        : PME electrostatics en GPU
#   -bonded gpu     : bonded interactions en GPU
#   -update gpu     : integracion de coordenadas en GPU
#                     (elimina transferencia CPU-GPU cada paso)
#   -maxwarn 1      : permite 1 warning (necesario con CHARMM-GUI)
#
# Con -update gpu la RTX 4090 maneja el loop completo de MD
# sin cuellos de botella CPU-GPU. En sistemas de ~150k atomos
# con membrana esto da ~150-200 ns/dia vs ~80-100 ns/dia sin -update gpu
# ==========================================================

# Detectar automaticamente en que replica estamos
REPLICA=$(basename $(pwd))
echo "=========================================================="
echo "  PIPELINE MD - $REPLICA"
echo "  Hardware: i9-14900K + RTX 4090"
echo "  Objetivo: 200 ns produccion x 3 pH"
echo "=========================================================="

for PH in ph50 ph70 ph80; do

    echo ""
    echo ">>> [$REPLICA / $PH] <<<"
    cd $PH

    # ----------------------------------------------------------
    # PASO 1: MINIMIZACION (Step 6.0)
    # La minimizacion corre en CPU - es rapida y no vale GPU
    # ----------------------------------------------------------
    if [ ! -f step6.0_minimization.gro ]; then
        echo "  [1/3] Minimizacion de energia..."
        gmx grompp \
            -f step6.0_minimization.mdp \
            -o step6.0_minimization.tpr \
            -c step5_input.gro \
            -r step5_input.gro \
            -p topol.top \
            -n index.ndx \
            -maxwarn 1

        gmx mdrun \
            -v \
            -deffnm step6.0_minimization \
            -ntmpi 1 \
            -ntomp 16
    else
        echo "  [1/3] Minimizacion ya completada. Saltando."
    fi

    # ----------------------------------------------------------
    # PASO 2: EQUILIBRADO (Steps 6.1 a 6.6)
    # NVT (6.1-6.3) y NPT (6.4-6.6) con restricciones gradualmente
    # reducidas sobre la proteina y la membrana
    # ----------------------------------------------------------
    echo "  [2/3] Equilibrado (6 fases)..."
    for cnt in 1 2 3 4 5 6; do
        pcnt=$((cnt - 1))
        ISTEP="step6.${cnt}_equilibration"
        if [ $cnt -eq 1 ]; then
            PSTEP="step6.0_minimization"
        else
            PSTEP="step6.${pcnt}_equilibration"
        fi

        if [ ! -f ${ISTEP}.gro ]; then
            echo "    -> Fase $cnt/6..."
            gmx grompp \
                -f ${ISTEP}.mdp \
                -o ${ISTEP}.tpr \
                -c ${PSTEP}.gro \
                -r step5_input.gro \
                -p topol.top \
                -n index.ndx \
                -maxwarn 1

            gmx mdrun \
                -v \
                -deffnm ${ISTEP} \
                -ntmpi 1 \
                -ntomp 16 \
                -gpu_id 0 \
                -nb gpu \
                -pme gpu \
                -bonded gpu
        else
            echo "    -> Fase $cnt/6 ya completada. Saltando."
        fi
    done

    # ----------------------------------------------------------
    # PASO 3: PRODUCCION 200 ns (Step 7)
    # -update gpu es la optimizacion mas importante:
    # mantiene posiciones y velocidades en la GPU entre pasos
    # eliminando la transferencia CPU<->GPU en cada frame
    # ----------------------------------------------------------
    if [ ! -f step7_1.gro ]; then
        echo "  [3/3] Produccion 200 ns (RTX 4090 full power)..."
        gmx grompp \
            -f step7_production.mdp \
            -o step7_1.tpr \
            -c step6.6_equilibration.gro \
            -t step6.6_equilibration.cpt \
            -p topol.top \
            -n index.ndx \
            -maxwarn 1

        gmx mdrun \
            -v \
            -deffnm step7_1 \
            -ntmpi 1 \
            -ntomp 16 \
            -gpu_id 0 \
            -nb gpu \
            -pme gpu \
            -bonded gpu \
            -update gpu
    else
        echo "  [3/3] Produccion ya completada. Saltando."
    fi

    cd ..
    echo ">>> [$REPLICA / $PH] COMPLETADO <<<"

done

echo ""
echo "=========================================================="
echo "  $REPLICA FINALIZADA"
echo "  Archivos de produccion generados:"
for PH in ph50 ph70 ph80; do
    if [ -f $PH/step7_1.gro ]; then
        FRAMES=$(gmx check -f $PH/step7_1.xtc 2>&1 | grep "frames" | awk '{print $NF}')
        echo "  [OK]  $PH/step7_1.xtc ($FRAMES frames)"
    else
        echo "  [!]   $PH/step7_1.xtc no encontrado"
    fi
done
echo ""
echo "  SIGUIENTE PASO:"
echo "  cd ~/Documentos/Proyecto_Tesis_Renato_MD"
echo "  ./1_analisis_estructural_v5.sh  (para analizar Replica_3)"
echo "=========================================================="
