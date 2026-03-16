#!/bin/bash
# gromacs_frames.sh ─ Extrae frames para ChimeraX
# Corre con: bash gromacs_frames.sh

cd ~/Documents/2TesisRenato/Proyecto_Tesis_Renato_MD

# Verificar rutas primero
echo "=== Verificando archivos ==="
ls Replica_1/ph70/*.tpr 2>/dev/null || echo "No hay .tpr en ph70"
ls Replica_1/ph70/*.xtc 2>/dev/null || echo "No hay .xtc en ph70"
echo "---"

# Frame 175 ns, pH 7.4, Replica 1
# El -center necesita que des 2 grupos: primero para centrar, luego para output
# Grupo 1 = proteína, Grupo 0 = sistema completo
printf "Protein\nSystem\n" | gmx trjconv \
  -s Replica_1/ph70/step7_1.tpr \
  -f Replica_1/ph70/step7_1.xtc \
  -o frame_175ns_ph74.pdb \
  -dump 175000 \
  -pbc mol \
  -center \
  -n Replica_1/ph70/index.ndx 2>&1 | tail -10

echo "=== pH 7.4 done ==="

# Frame 175 ns, pH 5.0, Replica 1
printf "Protein\nSystem\n" | gmx trjconv \
  -s Replica_1/ph50/step7_1.tpr \
  -f Replica_1/ph50/step7_1.xtc \
  -o frame_175ns_ph50.pdb \
  -dump 175000 \
  -pbc mol \
  -center \
  -n Replica_1/ph50/index.ndx 2>&1 | tail -10

echo "=== pH 5.0 done ==="
echo "Archivos generados: frame_175ns_ph74.pdb  frame_175ns_ph50.pdb"
