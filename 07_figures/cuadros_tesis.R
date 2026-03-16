#!/usr/bin/env Rscript
# cuadros_tesis.R  ─  Genera tablas formateadas como PNG para la tesis
# Requiere: ggplot2, dplyr, tidyr
# Corre con: Rscript cuadros_tesis.R
# Output: 4 archivos PNG de tablas

library(ggplot2)
library(dplyr)
library(tidyr)

# ══════════════════════════════════════════════════════════════════════════
# Helper: dibujar una tabla como figura con ggplot2 (sin paquetes extra)
# ══════════════════════════════════════════════════════════════════════════
draw_table <- function(df, title, subtitle=NULL, header_fill="#2E74B5",
                       col_widths=NULL, outfile, w=8, h=NULL) {

  nc <- ncol(df); nr <- nrow(df)
  if (is.null(col_widths)) col_widths <- rep(1/nc, nc)
  col_widths <- col_widths / sum(col_widths)  # normalizar

  # Calcular x centres para cada columna
  x_centers <- cumsum(col_widths) - col_widths/2

  # Altura de fila (proporcional al nro de filas)
  row_h <- 1/(nr+2)  # +1 header +1 margenes
  if (is.null(h)) h <- max(3.5, (nr+2)*0.45 + 1.5)

  # Construir data frame de celdas
  cells <- data.frame()

  # Header
  for (j in 1:nc) {
    cells <- rbind(cells, data.frame(
      x=x_centers[j], y=1-row_h*0.5, label=names(df)[j],
      bold=TRUE, color="white", fill=header_fill,
      row_type="header", col=j, row=0
    ))
  }

  # Filas de datos
  for (i in 1:nr) {
    fill_color <- if (i%%2==0) "#EBF3FA" else "white"
    for (j in 1:nc) {
      cells <- rbind(cells, data.frame(
        x=x_centers[j], y=1-row_h*(i+0.5), label=as.character(df[i,j]),
        bold=FALSE, color="gray10", fill=fill_color,
        row_type="data", col=j, row=i
      ))
    }
  }

  # Rectángulos (fondo de celdas)
  rects <- data.frame()
  for (i in 0:nr) {
    fill_color <- if (i==0) header_fill else if (i%%2==0) "#EBF3FA" else "white"
    for (j in 1:nc) {
      rects <- rbind(rects, data.frame(
        xmin=ifelse(j==1,0,cumsum(col_widths)[j-1]),
        xmax=cumsum(col_widths)[j],
        ymin=1-row_h*(i+1), ymax=1-row_h*i,
        fill=fill_color
      ))
    }
  }

  p <- ggplot() +
    # Fondos
    geom_rect(data=rects,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill),
              color="white", linewidth=0.5) +
    scale_fill_identity() +
    # Borde exterior
    annotate("rect", xmin=0, xmax=1,
             ymin=1-row_h*(nr+1), ymax=1,
             fill=NA, color="gray60", linewidth=0.6) +
    # Texto de celdas
    geom_text(data=filter(cells, row_type=="header"),
              aes(x=x, y=y, label=label),
              color="white", fontface="bold", size=3.4,
              hjust=0.5, vjust=0.5, lineheight=1.1) +
    geom_text(data=filter(cells, row_type=="data"),
              aes(x=x, y=y, label=label),
              color="gray10", size=3.2,
              hjust=0.5, vjust=0.5, lineheight=1.1) +
    # Título
    annotate("text", x=0, y=1+0.04, label=title,
             hjust=0, vjust=0, size=4.2, fontface="bold", color="gray10") +
    {if (!is.null(subtitle)) annotate("text", x=0, y=1+0.01, label=subtitle,
             hjust=0, vjust=0, size=3.0, color="gray50") } +
    coord_cartesian(xlim=c(-0.01,1.01),
                    ylim=c(1-row_h*(nr+1)-0.02, 1.08),
                    clip="off") +
    theme_void() +
    theme(plot.margin=margin(12,10,10,10))

  ggsave(outfile, p, width=w, height=h, dpi=300, bg="white")
  cat("Guardado:", outfile, "\n")
}

# ══════════════════════════════════════════════════════════════════════════
# CUADRO 1: Parámetros del protocolo de simulación
# ══════════════════════════════════════════════════════════════════════════
d1 <- data.frame(
  "Parámetro"   = c("dt (paso integración)","nsteps (pasos totales)","Duración / réplica",
                    "Frames totales","Termostato","Temperatura referencia",
                    "Barostato","Presión referencia","Electrostática largo alcance",
                    "Radio de corte (rcoulomb / rvdw)","Campo de fuerzas",
                    "Membrana (lípido)","Agua / iones","Réplicas por pH","Total simulado"),
  "Valor"       = c("0.002 ps (2 fs)","100,000,000","200 ns",
                    "4,000 (1 frame/50 ps)","v-rescale","310.15 K (37 °C)",
                    "C-rescale semi-isotrópico","1.0 bar","PME (Particle Mesh Ewald)",
                    "1.2 nm","CHARMM36m","POPC (140 lípidos)","TIP3P / NaCl 150 mM",
                    "3 (independientes)","1800 ns (3 réplicas × 3 pH × 200 ns)"),
  "Descripción" = c("Integración numérica","200 ns de producción","Por réplica",
                    "Resolución temporal 50 ps","3 grupos independientes","37 °C fisiológico",
                    "Correcto para bicapa lipídica","Isobárico","Alta precisión electrostática",
                    "Estándar CHARMM36m","Optimizado proteína-membrana","1-palmitoil-2-oleoil-PC",
                    "Concentración fisiológica","Semillas distintas","9 sistemas en total"),
  check.names=FALSE
)

draw_table(d1,
  title    = "Cuadro 1. Parámetros del protocolo de simulación (GROMACS 2025.3)",
  subtitle = "Confirmados de los archivos step7_production.mdp generados por CHARMM-GUI",
  outfile  = "cuadro_1_protocolo.png",
  col_widths = c(0.28, 0.30, 0.42), w=10, h=7.5)

# ══════════════════════════════════════════════════════════════════════════
# CUADRO 2: Energía libre de unión MM-GBSA
# ══════════════════════════════════════════════════════════════════════════
d2 <- data.frame(
  "Condición"           = c("Endosomal tardío","Sistémico / plasmático","Gástrico / NAC-SNAC"),
  "pH"                  = c("5.0","7.4","8.0"),
  "ΔG binding (kcal·mol⁻¹)" = c("−105.60","−131.26","−151.35"),
  "DE"                  = c("33.41","42.70","27.18"),
  "SEM"                 = c("19.29","24.65","15.70"),
  "Puentes H (media ± DE)" = c("9.7 ± 3.0","15.2 ± 5.5","18.1 ± 4.4"),
  check.names=FALSE
)

draw_table(d2,
  title    = "Cuadro 2. Energía libre de unión MM-GBSA del complejo semaglutida-GLP-1R",
  subtitle = "igb=5, NaCl 150 mM, frames 150-200 ns, n = 3 réplicas por condición de pH",
  outfile  = "cuadro_2_mmgbsa.png",
  col_widths = c(0.22, 0.09, 0.22, 0.10, 0.10, 0.27), w=11, h=3.5)

# ══════════════════════════════════════════════════════════════════════════
# CUADRO 3: Estado de protonación histidinas pH-sensibles
# ══════════════════════════════════════════════════════════════════════════
d3 <- data.frame(
  "Residuo"         = c("His99","His171","His173","His180","His212","His363","His374"),
  "Localización"    = c("N-terminal extracelular","N-terminal extracelular",
                        "N-terminal extracelular","ECL1 – contacto ligando",
                        "ECL2","TM6 (transmembrana)","TM6 / TM7"),
  "pH 5.0"          = c("HSP (+1)","HSP (+1)","HSP (+1)","HSP (+1)","HSP (+1)","HSP (+1)","HSP (+1)"),
  "pH 7.4"          = c("HSD (neutro)","HSD (neutro)","HSD (neutro)","HSD (neutro)",
                        "HSD (neutro)","HSD (neutro)","HSD (neutro)"),
  "pH 8.0"          = c("HSD (neutro)","HSD (neutro)","HSD (neutro)","HSD (neutro)",
                        "HSD (neutro)","HSD (neutro)","HSD (neutro)"),
  "pKa estimado"    = c("5.0–7.4","5.0–7.4","5.0–7.4","5.0–7.4","5.0–7.4","5.0–7.4","5.0–7.4"),
  check.names=FALSE
)

draw_table(d3,
  title    = "Cuadro 3. Estado de protonación de las histidinas pH-sensibles del GLP-1R (PropKa 3.0)",
  subtitle = "HSP = His doblemente protonada (carga +1)  |  HSD = His neutra  |  ECL = loop extracelular  |  TM = hélice transmembrana",
  outfile  = "cuadro_3_histidinas.png",
  col_widths = c(0.12, 0.28, 0.14, 0.14, 0.14, 0.18), w=11, h=4.5)

# ══════════════════════════════════════════════════════════════════════════
# CUADRO 4: Hotspots farmacológicos – residuos clave del paper
# ══════════════════════════════════════════════════════════════════════════
d4 <- data.frame(
  "Residuo"         = c("MET204","LEU244","LEU183","HIS180*","LEU166","MET233",
                        "TYR205","TRP297","GLU387*","HIS363*","TRP306","LEU359*"),
  "Localización"    = c("ECL1-TM2","TM3","TM2","ECL1","TM2","TM3",
                        "ECL1-TM2","ECL2","TM6 (profundo)","TM6","ECL2-TM4","TM5-TM6"),
  "ΔG pH 5.0"       = c("−6.76","−5.74","−5.87","0.00","−5.68","−6.05",
                        "−6.70","−3.25","−20.96","0.00","−2.14","−4.19"),
  "ΔG pH 7.4"       = c("−6.67","−6.11","−5.98","−6.04","−5.89","−6.18",
                        "−8.02","−3.91","−51.94","−6.06","−1.86","−5.22"),
  "ΔG pH 8.0"       = c("−7.08","−6.11","−6.00","−6.30","−6.21","−6.80",
                        "−7.84","−3.71","−53.98","−5.93","−2.43","−5.30"),
  "pH-sensible"     = c("No","No","No","Sí","No","No","Sí","No","Sí","Sí","No","Sí"),
  check.names=FALSE
)

draw_table(d4,
  title    = "Cuadro 4. Hotspots farmacológicos del GLP-1R – contribución VdW por residuo (kcal·mol⁻¹)",
  subtitle = "* pH-sensible (|ΔΔG| > 1.0 kcal·mol⁻¹ entre pH 5.0 y 7.4)  |  Valores VdW medios de n=9 trayectorias (3 réplicas × 3 pH)",
  outfile  = "cuadro_4_hotspots.png",
  col_widths = c(0.13, 0.20, 0.14, 0.14, 0.14, 0.15),
  w=11, h=5.8)

cat("\nListo. Cuadros generados:\n")
cat("  cuadro_1_protocolo.png\n")
cat("  cuadro_2_mmgbsa.png\n")
cat("  cuadro_3_histidinas.png\n")
cat("  cuadro_4_hotspots.png\n")
