#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════
# run_chimerax.sh  –  Genera las 5 figuras ChimeraX automáticamente
# Corre con:  bash run_chimerax.sh
# ═══════════════════════════════════════════════════════════════════════

cd ~/Documents/2TesisRenato/Proyecto_Tesis_Renato_MD

# Verificar que existen los PDB
if [ ! -f frame_175ns_ph74.pdb ] || [ ! -f frame_175ns_ph50.pdb ]; then
    echo "ERROR: No se encuentran los PDB. Asegúrate de tener:"
    echo "  frame_175ns_ph74.pdb"
    echo "  frame_175ns_ph50.pdb"
    exit 1
fi

echo "Generando figuras ChimeraX..."
echo "Esto abrirá ChimeraX, correrá el script y cerrará automáticamente."
echo "(puede tardar 1-2 minutos)"

# Copiar el script al directorio de trabajo
cp ~/Downloads/chimerax_complete.cxc . 2>/dev/null || true

# Correr ChimeraX con el script
# --offscreen usa rendering sin display (más rápido, misma calidad)
chimerax --offscreen --script chimerax_complete.cxc

echo ""
echo "Figuras generadas:"
ls -lh fig1_sistema_general.png fig2_histidinas_pH74.png \
        fig3_histidinas_pH50.png fig4_GLU387_closeup.png \
        fig5_nucleo_hidrofobico.png 2>/dev/null

