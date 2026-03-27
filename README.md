# Semaglutida–GLP-1R Molecular Dynamics Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![GROMACS](https://img.shields.io/badge/GROMACS-2025.3-orange)](https://www.gromacs.org/)
[![Python](https://img.shields.io/badge/Python-3.10+-green)](https://www.python.org/)

**Efecto del pH sobre la afinidad de unión del complejo semaglutida–GLP-1R: estudio in silico de la plausibilidad de señalización endosomal y su relación con la tolerancia clínica**

Vásquez Laime, Renato Junior - ORCID: [0009-0007-1494-4685](https://orcid.org/0009-0007-1494-4685)  
CBCRG - Universidad Católica de Santa María, Arequipa, Perú - 2026

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
Tesis-Medicina-Semaglutide/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── Archivos_Base_CHARMM/               # Archivos de entrada CHARMM-GUI (~90 MB)
│   ├── ph50/                           # pH 5.0 — His protonadas (HSP +1)
│   │   ├── step5_input.gro             #   Coordenadas iniciales
│   │   ├── step5_input.pdb             #   Estructura PDB
│   │   ├── step5_input.psf             #   Topología PSF
│   │   ├── topol.top                   #   Topología GROMACS
│   │   ├── index.ndx                   #   Índice de grupos
│   │   ├── step6.0_minimization.mdp    #   Minimización
│   │   ├── step6.{1-6}_equilibration.mdp  # Equilibrado NVT/NPT (6 fases)
│   │   ├── step7_production.mdp        #   Producción 200 ns
│   │   └── toppar/                     #   Parámetros de campo de fuerza
│   │       ├── forcefield.itp          #     CHARMM36m
│   │       ├── PROA.itp               #     Semaglutida
│   │       ├── PROB.itp               #     GLP-1R
│   │       ├── POPC.itp               #     Membrana
│   │       ├── TIP3.itp               #     Agua
│   │       ├── SOD.itp                #     Na+
│   │       └── CLA.itp               #     Cl-
│   ├── ph70/                           # pH 7.4 — His neutras (HSD)
│   └── ph80/                           # pH 8.0 — His neutras (HSD)
│
├── 01_produccion_md.sh                 # Corrida MD (3 réplicas × 3 pH)
├── 02_analisis_estructural.sh          # RMSD, RMSF, HBond, Rg, convergencia
├── 03_mmgbsa_energia.sh                # Energía libre MM-GBSA (igb=5)
├── 04_mmgbsa_decomp.sh                 # Descomposición por residuo (idecomp=2)
├── 05_analisis_hotspots.py             # Hotspots farmacológicos
├── 06_extraer_frames.sh                # Frames PDB para visualización
├── 07_graficos.R                       # Todas las figuras de la tesis
└── 08_tablas.py                        # Todas las tablas Word de la tesis
```

> **Reproducibilidad:** los archivos de `Archivos_Base_CHARMM/` son suficientes para reproducir todas las simulaciones desde cero. Las trayectorias (`.xtc`, ~50 GB/réplica) no están incluidas por tamaño — ver sección [Datos de trayectoria](#datos-de-trayectoria).

---

## Requisitos

| Herramienta | Versión | Uso |
|---|---|---|
| [GROMACS](https://www.gromacs.org/) | 2025.3 | Simulación MD |
| [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/) | 1.6.4 | Energía libre MM-GBSA |
| [PropKa](https://github.com/jensengroup/propka) | 3.0 | Estados de protonación |
| Python | ≥ 3.10 | Análisis hotspots (`pandas`, `numpy`, `python-docx`) |
| R | ≥ 4.3 | Figuras (`ggplot2`, `dplyr`, `tidyr`) |

**Hardware utilizado:** Intel Core i9-14900K + NVIDIA RTX 4090 24 GB (CBCRG-UCSM)

---

## Uso — Pipeline completo

### Paso 0: Preparar directorios de réplicas

Los archivos de `Archivos_Base_CHARMM/` contienen los sistemas construidos con CHARMM-GUI Membrane Builder a partir de PDB 7KI0 (crioEM, 2.21 Å). Copiarlos a cada réplica:

```bash
for rep in Replica_1 Replica_2 Replica_3; do
    mkdir -p $rep
    cp -r Archivos_Base_CHARMM/ph50 $rep/
    cp -r Archivos_Base_CHARMM/ph70 $rep/
    cp -r Archivos_Base_CHARMM/ph80 $rep/
done
```

Cada réplica usa los mismos archivos de topología pero diferente semilla de velocidades en el equilibrado.

### Paso 1: Simulación MD (200 ns por réplica/pH)

```bash
bash 01_produccion_md.sh
```

Ejecuta `gmx mdrun` para las 3 réplicas × 3 pH. Detecta checkpoints para continuar corridas interrumpidas.

### Paso 2: Análisis estructural

```bash
bash 02_analisis_estructural.sh
```

Genera por réplica/pH: índice CHARMM por segmento, trayectoria centrada, RMSD (receptor, ligando, backbone), RMSF por residuo, puentes de hidrógeno interfaz, radio de giro, verificación de convergencia por ventanas, y validación RMSD vs estructura cristalográfica.

### Paso 3: Energía libre de unión MM-GBSA

```bash
conda activate gmxMMPBSA
bash 03_mmgbsa_energia.sh
```

Parámetros: `igb=5`, `saltcon=0.150 M`, frames 3000–4000 (150–200 ns), intervalo 10. ~1–2 h por condición.

### Paso 4: Descomposición por residuo

```bash
bash 04_mmgbsa_decomp.sh
```

`idecomp=2` (pairwise per-residue). Detecta automáticamente residuos en la interfaz ≤ 6 Å. ~3–5 h por condición.

### Paso 5: Análisis de hotspots

```bash
python3 05_analisis_hotspots.py
```

Umbral: `TOTAL_avg < −1.5 kcal·mol⁻¹`. Genera tablas CSV de hotspots por pH, tabla pivote y residuos consistentes en los 3 pH.

### Paso 6: Frames para visualización

```bash
bash 06_extraer_frames.sh
```

Extrae frames PDB a t = 175 ns (punto medio de la ventana MM-GBSA) para las 9 condiciones.

### Paso 7: Figuras de la tesis

```bash
Rscript 07_graficos.R
```

Genera: ΔG barras, RMSD vs tiempo, puentes H vs tiempo, RMSF, heatmap VdW hotspots, lollipop VdW, validación RMSD vs cristal, y cuadros 1–4 como PNG.

### Paso 8: Tablas Word de la tesis

```bash
python3 08_tablas.py
```

Genera: Cuadros 1–4 (protocolo, energía, histidinas, hotspots), ANOVA 5A/5B/5C, y Anexo 11 (pKa PropKa 3.0).

---

## Principios de organización

1. **Las carpetas de réplica NUNCA se contaminan** con outputs de análisis
2. **Todos los outputs van a `Analisis_Global/`** separados por réplica, pH y tipo
3. **Los archivos CHARMM-GUI originales** se preservan en `Archivos_Base_CHARMM/`
4. **Cada script detecta resultados previos** y no recalcula lo que ya existe
5. **Todas las rutas son relativas** al directorio del proyecto

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
| His99 | N-terminal extracelular | HSP (+1) | HSD (neutro) | 6.16 |
| His171 | N-terminal / TM1 | HSP (+1) | HSD (neutro) | 7.08 |
| His173 | N-terminal extracelular | HSP (+1) | HSD (neutro) | 5.97 |
| His180 | ECL1 – contacto ligando | HSP (+1) | HSD (neutro) | 6.35 |
| His212 | ECL2 | HSP (+1) | HSD (neutro) | 6.37 |
| His363 | TM6 | HSP (+1) | HSD (neutro) | 7.95 |
| His374 | TM6 / TM7 | HSP (+1) | HSD (neutro) | 6.35 |

---

## Datos de trayectoria

Las trayectorias completas (`.xtc`, ~50 GB por réplica) no están incluidas en este repositorio. Los archivos de `Archivos_Base_CHARMM/` son suficientes para reproducir el trabajo desde cero ejecutando el pipeline en orden.

---

## Cita

```bibtex
@thesis{vasquez2026semaglutida,
  author  = {Vásquez Laime, Renato Junior},
  title   = {Efecto del pH sobre la afinidad de unión del complejo
             semaglutida--GLP-1R: estudio in silico de la plausibilidad
             de señalización endosomal y su relación con la tolerancia clínica},
  school  = {Universidad Católica de Santa María},
  year    = {2026},
  address = {Arequipa, Perú},
  url     = {https://github.com/Mortis37/Tesis-Medicina-Semaglutide}
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
