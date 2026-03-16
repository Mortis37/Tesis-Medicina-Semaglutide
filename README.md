# Semaglutida–GLP-1R Molecular Dynamics Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![GROMACS](https://img.shields.io/badge/GROMACS-2025.3-orange)](https://www.gromacs.org/)
[![Python](https://img.shields.io/badge/Python-3.10+-green)](https://www.python.org/)

**Efecto del pH sobre la afinidad de unión del complejo semaglutida–GLP-1R: estudio in silico de la plausibilidad de señalización endosomal y su relación con la tolerancia clínica**

Vásquez Laime, Renato Junior — ORCID: [0009-0007-1494-4685](https://orcid.org/0009-0007-1494-4685)   
CBCRG — Universidad Católica de Santa María, Arequipa, Perú — 2026

---

## Resumen del estudio

Se evaluó in silico el efecto de tres condiciones de pH biológicamente relevantes (5.0 endosomal, 7.4 sistémico, 8.0 inducido por NAC/SNAC) sobre la energía libre de unión y la estabilidad estructural del complejo semaglutida–GLP-1R en membrana POPC explícita (1,800 ns totales · 9 trayectorias).

**Resultados principales:**
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
├── 01_simulation/              # Pipeline de simulación MD
│   └── correr_replica.sh       # Minimización → equilibrado → producción 200 ns
│
├── 02_structural_analysis/     # Análisis estructural
│   └── 1_analisis_estructural_v6.sh   # RMSD, RMSF, H-bonds, Rg
│
├── 03_mmgbsa_energy/           # Energía libre de unión
│   └── 2_mmgbsa_energia.sh     # MM-GBSA global (igb=5, NaCl 150 mM)
│
├── 04_mmgbsa_decomp/           # Descomposición por residuo
│   └── 3_mmgbsa_decomp.sh      # idecomp=2 (pairwise per-residue)
│
├── 05_hotspots/                # Identificación de hotspots farmacológicos
│   └── 4_analisis_hotspots.py  # Parseo de CSV + tablas resumen
│
├── 06_visualization/           # Visualización molecular
│   ├── chimerax_v12.cxc        # Script ChimeraX (figuras de tesis)
│   ├── gromacs_frames.sh       # Extracción de frames PDB para ChimeraX
│   └── run_chimerax.sh         # Lanzador de figuras automático
│
└── 07_figures/                 # Figuras y cuadros de la tesis
    └── cuadros_tesis.R         # Cuadros 1–4 como PNG (ggplot2)
```

---

## Requisitos

### Software
| Herramienta | Versión | Uso |
|---|---|---|
| [GROMACS](https://www.gromacs.org/) | 2025.3 | Simulación MD |
| [CHARMM-GUI](https://charmm-gui.org/) | — | Construcción del sistema |
| [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/) | 1.6.4 | Energía libre MM-GBSA |
| [PropKa](https://github.com/jensengroup/propka) | 3.0 | Estados de protonación |
| [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) | 1.8 | Visualización molecular |
| Python | ≥ 3.10 | Análisis hotspots |
| R | ≥ 4.3 | Figuras y cuadros |

### Hardware utilizado
- **CPU:** Intel Core i9-14900K (32 hilos)
- **GPU:** NVIDIA RTX 4090 24 GB
- **Institución:** CBCRG-UCSM, Arequipa, Perú

### Dependencias Python
```bash
pip install pandas numpy
```

### Dependencias R
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork"))
```

---

## Uso — Pipeline completo

### Paso 0: Preparar el sistema
El sistema molecular (semaglutida-GLP-1R en membrana POPC) fue construido con **CHARMM-GUI Membrane Builder** a partir de la estructura crioEM **PDB: 7KI0** (resolución 2.21 Å).  
Los estados de protonación de las histidinas se asignaron con **PropKa 3.0** antes de cada condición de pH.

Estructura del directorio de trabajo esperada:
```
Proyecto_Tesis_Renato_MD/
├── Replica_1/
│   ├── ph50/   (step5_input.gro, topol.top, step6.*.mdp, step7_production.mdp)
│   ├── ph70/
│   └── ph80/
├── Replica_2/
│   └── ...
└── Replica_3/
    └── ...
```

### Paso 1: Simulación MD
```bash
cd Replica_3
bash 01_simulation/correr_replica.sh
```
Ejecuta: minimización (Step 6.0) → equilibrado en 6 fases (Steps 6.1–6.6) → producción 200 ns (Step 7).

**Flags clave de `gmx mdrun` para RTX 4090:**
```
-ntmpi 1 -ntomp 16 -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu
```
Con `-update gpu` la GPU maneja el loop completo de MD sin transferencias CPU↔GPU por frame (~150–200 ns/día para sistemas de ~90,000 átomos con membrana).

### Paso 2: Análisis estructural
```bash
bash 02_structural_analysis/1_analisis_estructural_v6.sh
```
Genera por cada réplica/pH:
- `RMSD/rmsd_receptor.xvg` — RMSD backbone GLP-1R
- `RMSD/rmsd_ligando_vs_receptor.xvg` — RMSD semaglutida vs receptor
- `RMSF/rmsf_receptor.xvg` / `rmsf_ligando.xvg`
- `HBond/hbonds_Rec_Lig.xvg`
- `Rg/rg_receptor.xvg` / `rg_complejo.xvg`

### Paso 3: Energía libre de unión MM-GBSA
```bash
conda activate gmxMMPBSA
bash 03_mmgbsa_energy/2_mmgbsa_energia.sh
```
Parámetros: `igb=5`, `saltcon=0.150 M`, frames `3000–4000` (150–200 ns), intervalo 10.

### Paso 4: Descomposición por residuo
```bash
bash 04_mmgbsa_decomp/3_mmgbsa_decomp.sh
```
Usa `idecomp=2` (pairwise per-residue). Calcula la contribución energética de cada residuo del receptor a la unión con semaglutida. **Advertencia:** 3–5 horas por condición.

### Paso 5: Análisis de hotspots
```bash
python3 05_hotspots/4_analisis_hotspots.py
```
Umbral hotspot: `TOTAL_avg < −1.5 kcal·mol⁻¹`. Genera CSVs con hotspots por pH, tabla pivote y residuos consistentes en los 3 pH.

### Paso 6: Extracción de frames y visualización
```bash
# Extraer frames PDB (t = 175 ns)
bash 06_visualization/gromacs_frames.sh

# Generar figuras de tesis con ChimeraX
bash 06_visualization/run_chimerax.sh
```

### Paso 7: Cuadros de la tesis (R)
```bash
Rscript 07_figures/cuadros_tesis.R
```
Genera `cuadro_1_protocolo.png`, `cuadro_2_mmgbsa.png`, `cuadro_3_histidinas.png`, `cuadro_4_hotspots.png`.

---

## Detalles del sistema simulado

| Componente | Detalle |
|---|---|
| Ligando (cadena A) | Semaglutida, residuos 7–36 |
| Receptor (cadena B) | GLP-1R, residuos 29–423 |
| Membrana | POPC, 140 lípidos |
| Agua | TIP3P |
| Iones | NaCl 150 mM |
| Campo de fuerzas | CHARMM36m |
| Tamaño del sistema | ~90,000 átomos |
| Estructura base | PDB 7KI0 (crioEM, 2.21 Å) |

### Histidinas pH-sensibles (PropKa 3.0)
| Residuo | Localización | pH 5.0 | pH 7.4 / 8.0 |
|---|---|---|---|
| His99 | N-terminal extracelular | HSP (+1) | HSD (neutro) |
| His171 | N-terminal extracelular | HSP (+1) | HSD (neutro) |
| His173 | N-terminal extracelular | HSP (+1) | HSD (neutro) |
| His180 | ECL1 – contacto ligando | HSP (+1) | HSD (neutro) |
| His212 | ECL2 | HSP (+1) | HSD (neutro) |
| His363 | TM6 | HSP (+1) | HSD (neutro) |
| His374 | TM6 / TM7 | HSP (+1) | HSD (neutro) |

---

## Reproducibilidad

> Este repositorio contiene los scripts exactos utilizados para generar los resultados publicados en la tesis.  
> Los archivos de trayectoria (`.xtc`, ~50 GB por réplica) no están incluidos por limitaciones de tamaño.  
> Los frames representativos (PDB) y los archivos `.xvg` de análisis están disponibles bajo solicitud.

Para reproducir exactamente los resultados:
1. Descargar PDB 7KI0 de [RCSB](https://www.rcsb.org/structure/7KI0)
2. Construir el sistema con CHARMM-GUI Membrane Builder (POPC, TIP3P, NaCl 150 mM)
3. Asignar estados de protonación con PropKa 3.0 para cada condición de pH
4. Ejecutar el pipeline en orden (Steps 1–7 arriba)

---

## Cita

Si usas este código en tu trabajo, por favor cita:

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
CBCRG — Universidad Católica de Santa María  
Arequipa, Perú
