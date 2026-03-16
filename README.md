# Semaglutida–GLP-1R Molecular Dynamics Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![GROMACS](https://img.shields.io/badge/GROMACS-2025.3-orange)](https://www.gromacs.org/)
[![Python](https://img.shields.io/badge/Python-3.10+-green)](https://www.python.org/)

**Efecto del pH sobre la afinidad de unión del complejo semaglutida–GLP-1R: estudio in silico de la plausibilidad de señalización endosomal y su relación con la tolerancia clínica**

Vásquez Laime, Renato Junior — ORCID: [0009-0007-1494-4685](https://orcid.org/0009-0007-1494-4685)  
Asesora: Dra. Muñoz del Carpio Toia, Águeda — ORCID: [0000-0003-0501-7314](https://orcid.org/0000-0003-0501-7314)  
CBCRG — Universidad Católica de Santa María, Arequipa, Perú — 2026

---

## Resumen del estudio

Se evaluó in silico el efecto de tres condiciones de pH (5.0 endosomal, 7.4 sistémico, 8.0 inducido por NAC/SNAC) sobre la energía libre de unión y la estabilidad estructural del complejo semaglutida–GLP-1R en membrana POPC explícita (1,800 ns totales · 9 trayectorias).

| Condición | pH | ΔG binding (kcal·mol⁻¹) | Puentes H (media ± DE) |
|---|---|---|---|
| Endosomal tardío | 5.0 | −105.60 ± 33.41 | 9.7 ± 3.0 |
| Sistémico / plasmático | 7.4 | −131.26 ± 42.70 | 15.2 ± 5.5 |
| Gástrico / NAC-SNAC | 8.0 | −151.35 ± 27.18 | 18.1 ± 4.4 |

---

## Estructura del repositorio

```
semaglutida-glp1r-md/
├── README.md
├── LICENSE
│
├── 00_charmm_base/                   # Archivos de entrada CHARMM-GUI (~98 MB)
│   ├── ph50/                         # pH 5.0 — His protonadas (HSP +1)
│   │   ├── step5_input.gro           #   Coordenadas iniciales del sistema
│   │   ├── step5_input.pdb           #   Estructura PDB
│   │   ├── topol.top                 #   Topología CHARMM36m
│   │   ├── index.ndx                 #   Índice de grupos GROMACS
│   │   ├── step6.0_minimization.mdp  #   Parámetros minimización
│   │   ├── step6.1_equilibration.mdp #   NVT fase 1
│   │   ├── ...
│   │   ├── step6.6_equilibration.mdp #   NPT fase 6 (sin restricciones)
│   │   └── step7_production.mdp      #   ← Producción 200 ns
│   ├── ph70/                         # pH 7.4 — His neutras (HSD)
│   └── ph80/                         # pH 8.0 — His neutras (HSD)
│
├── 01_simulation/                    # Pipeline de simulación MD
│   └── correr_replica.sh
│
├── 02_structural_analysis/           # RMSD, RMSF, H-bonds, radio de giro
│   └── 1_analisis_estructural_v6.sh
│
├── 03_mmgbsa_energy/                 # Energía libre de unión MM-GBSA
│   └── 2_mmgbsa_energia.sh
│
├── 04_mmgbsa_decomp/                 # Descomposición por residuo (idecomp=2)
│   └── 3_mmgbsa_decomp.sh
│
├── 05_hotspots/                      # Análisis de hotspots farmacológicos
│   └── 4_analisis_hotspots.py
│
├── 06_visualization/                 # Visualización molecular
│   ├── chimerax_v12.cxc              #   Script ChimeraX (figuras tesis)
│   ├── gromacs_frames.sh             #   Extracción de frames PDB
│   └── run_chimerax.sh
│
└── 07_figures/                       # Cuadros de la tesis
    └── cuadros_tesis.R               #   Cuadros 1–4 como PNG (ggplot2)
```

> **Reproducibilidad:** los archivos de `00_charmm_base/` son suficientes para reproducir todas las simulaciones desde cero. Las trayectorias (`.xtc`, ~50 GB/réplica) no están incluidas por tamaño — ver sección [Datos de trayectoria](#datos-de-trayectoria).

---

## Requisitos

| Herramienta | Versión | Uso |
|---|---|---|
| [GROMACS](https://www.gromacs.org/) | 2025.3 | Simulación MD |
| [CHARMM-GUI](https://charmm-gui.org/) | — | Sistema ya construido en `00_charmm_base/` |
| [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/) | 1.6.4 | Energía libre MM-GBSA |
| [PropKa](https://github.com/jensengroup/propka) | 3.0 | Estados de protonación |
| [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) | 1.8 | Visualización molecular |
| Python | ≥ 3.10 | Análisis hotspots (`pandas`, `numpy`) |
| R | ≥ 4.3 | Figuras (`ggplot2`, `dplyr`, `patchwork`) |

**Hardware utilizado:** Intel Core i9-14900K + NVIDIA RTX 4090 24 GB (CBCRG-UCSM)

---

## Uso — Pipeline completo

### Paso 0: Preparar directorios de réplicas

Los archivos de `00_charmm_base/` contienen los sistemas ya construidos con CHARMM-GUI Membrane Builder a partir de PDB 7KI0 (crioEM, 2.21 Å). Copiarlos a cada réplica:

```bash
for rep in Replica_1 Replica_2 Replica_3; do
    mkdir -p $rep
    cp -r 00_charmm_base/ph50 $rep/
    cp -r 00_charmm_base/ph70 $rep/
    cp -r 00_charmm_base/ph80 $rep/
done
```

Cada réplica usa los mismos archivos de topología pero diferente semilla de velocidades en el equilibrado (configurada en `step6.1_equilibration.mdp`).

### Paso 1: Simulación MD (200 ns por réplica/pH)

```bash
cd Replica_1
bash ../01_simulation/correr_replica.sh
```

Ejecuta para cada pH: minimización → 6 fases de equilibrado → producción 200 ns.

**Flags clave para RTX 4090** (en `gmx mdrun`):
```
-ntmpi 1 -ntomp 16 -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu
```
Con `-update gpu` el loop completo de MD corre en GPU sin transferencias CPU↔GPU (~150–200 ns/día para ~90,000 átomos).

### Paso 2: Análisis estructural

```bash
bash 02_structural_analysis/1_analisis_estructural_v6.sh
```

Genera por réplica/pH: RMSD del receptor y ligando, RMSF por residuo, puentes de hidrógeno interfaz, radio de giro. Output en `.xvg` para cada condición.

### Paso 3: Energía libre de unión MM-GBSA

```bash
conda activate gmxMMPBSA
bash 03_mmgbsa_energy/2_mmgbsa_energia.sh
```

Parámetros: `igb=5`, `saltcon=0.150 M`, frames `3000–4000` (150–200 ns), intervalo 10. ~1–2 horas por condición.

### Paso 4: Descomposición por residuo

```bash
bash 04_mmgbsa_decomp/3_mmgbsa_decomp.sh
```

`idecomp=2` (pairwise per-residue). Identifica hotspots farmacológicos. ~3–5 horas por condición.

### Paso 5: Análisis de hotspots

```bash
python3 05_hotspots/4_analisis_hotspots.py
```

Umbral: `TOTAL_avg < −1.5 kcal·mol⁻¹`. Genera tablas CSV de hotspots por pH, tabla pivote y residuos consistentes en los 3 pH.

### Paso 6: Visualización molecular

```bash
bash 06_visualization/gromacs_frames.sh   # Extrae frames PDB (t = 175 ns)
bash 06_visualization/run_chimerax.sh     # Genera figuras de tesis
```

### Paso 7: Cuadros de la tesis

```bash
Rscript 07_figures/cuadros_tesis.R
```

---

## Sistema simulado

| Componente | Detalle |
|---|---|
| Ligando (cadena A / PROA) | Semaglutida, residuos 7–36 |
| Receptor (cadena B / PROB) | GLP-1R, residuos 29–423 |
| Membrana | POPC, 140 lípidos |
| Agua | TIP3P |
| Iones | NaCl 150 mM |
| Campo de fuerzas | CHARMM36m |
| Tamaño del sistema | ~90,000 átomos |
| Estructura base | PDB 7KI0 (crioEM, 2.21 Å) |

### Histidinas pH-sensibles (PropKa 3.0)

| Residuo | Localización | pH 5.0 | pH 7.4 / 8.0 | pKa estimado |
|---|---|---|---|---|
| His99 | N-terminal extracelular | HSP (+1) | HSD (neutro) | 5.0–7.4 |
| His171 | N-terminal extracelular | HSP (+1) | HSD (neutro) | 5.0–7.4 |
| His173 | N-terminal extracelular | HSP (+1) | HSD (neutro) | 5.0–7.4 |
| His180 | ECL1 – contacto ligando | HSP (+1) | HSD (neutro) | 5.0–7.4 |
| His212 | ECL2 | HSP (+1) | HSD (neutro) | 5.0–7.4 |
| His363 | TM6 | HSP (+1) | HSD (neutro) | 5.0–7.4 |
| His374 | TM6 / TM7 | HSP (+1) | HSD (neutro) | 5.0–7.4 |

---

## Datos de trayectoria

Las trayectorias completas (`.xtc`, ~50 GB por réplica) no están incluidas en este repositorio. Los archivos de `00_charmm_base/` son suficientes para reproducir el trabajo desde cero ejecutando el pipeline en orden.

Los archivos `.xvg` de análisis y frames representativos en PDB están disponibles bajo solicitud al autor.

---

## Cita

```bibtex
@thesis{vasquez2026semaglutida,
  author  = {Vásquez Laime, Renato Junior},
  title   = {Efecto del pH sobre la afinidad de unión del complejo semaglutida--GLP-1R:
             estudio in silico de la plausibilidad de señalización endosomal
             y su relación con la tolerancia clínica},
  school  = {Universidad Católica de Santa María},
  year    = {2026},
  address = {Arequipa, Perú},
  url     = {https://github.com/renato-vasquez/semaglutida-glp1r-md}
}
```

---

## Licencia

MIT License — ver [LICENSE](LICENSE) para detalles.

---

## Contacto

**Renato Junior Vásquez Laime**  
ORCID: [0009-0007-1494-4685](https://orcid.org/0009-0007-1494-4685)  
CBCRG — Universidad Católica de Santa María, Arequipa, Perú
