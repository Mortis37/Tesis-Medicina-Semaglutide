#!/usr/bin/env Rscript
# ==========================================================
# 07_graficos.R
# Genera todas las figuras de la tesis
#
# Figuras:
#   fig_deltaG.png            — DeltaG binding por pH (barras)
#   fig_rmsd.png              — RMSD backbone + ligando vs tiempo
#   fig_hbonds.png            — Puentes H receptor-ligando vs tiempo
#   fig_rmsf.png              — RMSF receptor + ligando
#   fig_hotspots_heatmap.png  — Heatmap VdW residuos clave
#   fig_hotspots_lollipop.png — Lollipop VdW top residuos
#   fig_rmsd_cristal.png      — Validacion RMSD vs PDB 7KI0
#   fig_rmsd_cristal_tabla.png
#   cuadro_1_protocolo.png    — Tabla protocolo simulacion
#   cuadro_2_mmgbsa.png       — Tabla energia libre
#   cuadro_3_histidinas.png   — Tabla histidinas pH-sensibles
#   cuadro_4_hotspots.png     — Tabla hotspots farmacologicos
#
# Input:  Analisis_Global/ (xvg) + Resultados_Hotspots/ (csv)
# Output: Analisis_Global/Graficos/
#
# Uso: Rscript 07_graficos.R
#      (ejecutar desde la raiz del proyecto)
# ==========================================================

library(ggplot2)
library(dplyr)
library(tidyr)

OUTDIR   <- "Analisis_Global/Graficos"
BASE_DIR <- "Analisis_Global"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(OUTDIR, f)

PH_LABS <- c("pH 5.0 (endosomal)", "pH 7.4 (sistemico)", "pH 8.0 (NAC/SNAC)")
PH_COLS <- c("pH 5.0 (endosomal)" = "#E74C3C",
             "pH 7.4 (sistemico)"  = "#2E86AB",
             "pH 8.0 (NAC/SNAC)"   = "#27AE60")
phs_dir <- c("ph50", "ph70", "ph80")


# ==================================================================
# HELPERS
# ==================================================================
read_xvg <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dl <- lines[!grepl("^[#@]", lines) & nchar(trimws(lines)) > 0]
  vals <- strsplit(trimws(dl), "\\s+")
  df <- data.frame(x = as.numeric(sapply(vals, `[`, 1)),
                   y = as.numeric(sapply(vals, `[`, 2)))
  df[complete.cases(df), ]
}

rollmean_k <- function(x, k = 81) {
  n <- length(x); out <- rep(NA_real_, n); half <- floor(k / 2)
  for (i in seq_len(n)) {
    lo <- max(1, i - half); hi <- min(n, i + half)
    out[i] <- mean(x[lo:hi], na.rm = TRUE)
  }
  out
}


# ==================================================================
# FIG: DeltaG binding por pH
# ==================================================================
df_dg <- data.frame(
  ph = factor(c("pH 5.0\n(Endosomal)", "pH 7.4\n(Sistemico)", "pH 8.0\n(NAC/SNAC)"),
              levels = c("pH 5.0\n(Endosomal)", "pH 7.4\n(Sistemico)", "pH 8.0\n(NAC/SNAC)")),
  mean = c(-105.60, -131.26, -151.35),
  sd   = c(33.41, 42.70, 27.18))

cols3 <- c("pH 5.0\n(Endosomal)" = "#E74C3C",
           "pH 7.4\n(Sistemico)"  = "#2E86AB",
           "pH 8.0\n(NAC/SNAC)"   = "#27AE60")

fig_dg <- ggplot(df_dg, aes(x = ph, y = mean, fill = ph)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.18, linewidth = 0.8, color = "gray25") +
  geom_text(aes(label = paste0(round(mean, 1), "\nkcal/mol"), y = mean / 2),
            size = 4, fontface = "bold", color = "white", lineheight = 1.1) +
  annotate("segment", x = 1, xend = 1, y = -8, yend = -16, linewidth = 0.5, color = "gray35") +
  annotate("segment", x = 2, xend = 2, y = -8, yend = -16, linewidth = 0.5, color = "gray35") +
  annotate("segment", x = 1, xend = 2, y = -16, yend = -16, linewidth = 0.5, color = "gray35") +
  annotate("text", x = 1.5, y = -20, label = "DDG ~ 26 kcal/mol",
           size = 3.5, color = "gray25", fontface = "italic", vjust = 1) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "gray55", linetype = "dashed") +
  scale_fill_manual(values = cols3) +
  scale_y_continuous(limits = c(-212, 12), breaks = seq(-200, 0, 50), expand = c(0, 0)) +
  labs(x = NULL,
       y = expression(Delta * "G binding (kcal/mol)"),
       title = "Alta afinidad residual a pH endosomal",
       subtitle = "Media +/- DE, n = 3 replicas") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4),
        plot.margin = margin(10, 20, 10, 10))
ggsave(out("fig_deltaG.png"), fig_dg, width = 6.5, height = 5.5, dpi = 300)
cat("OK: fig_deltaG.png\n")


# ==================================================================
# FIG: RMSD
# ==================================================================
time_grid <- seq(0, 200, by = 0.1)

collect_rmsd <- function(mf) {
  result <- list()
  for (i in seq_along(phs_dir)) {
    ph <- phs_dir[i]
    mat <- matrix(NA, nrow = length(time_grid), ncol = 3)
    for (r in 1:3) {
      path <- file.path(BASE_DIR, paste0("Replica_", r), ph, "RMSD", mf)
      if (file.exists(path)) {
        d <- read_xvg(path)
        mat[, r] <- approx(d$x, d$y, xout = time_grid, rule = 2)$y
      }
    }
    n_rep <- max(sum(!is.na(mat[1, ])), 1)
    result[[i]] <- data.frame(
      time = time_grid,
      mean = rollmean_k(rowMeans(mat, na.rm = TRUE), k = 81),
      sem  = rollmean_k(apply(mat, 1, sd, na.rm = TRUE) / sqrt(n_rep), k = 81),
      ph_label = PH_LABS[i])
  }
  bind_rows(result)
}

df_rmsd <- bind_rows(
  collect_rmsd("rmsd_backbone.xvg") %>% mutate(panel = "Receptor GLP-1R (backbone)"),
  collect_rmsd("rmsd_ligando_vs_receptor.xvg") %>% mutate(panel = "Semaglutida (vs receptor)")) %>%
  mutate(ph_label = factor(ph_label, levels = PH_LABS),
         panel = factor(panel, levels = c("Receptor GLP-1R (backbone)", "Semaglutida (vs receptor)")),
         mean_A = mean * 10, sem_A = sem * 10)

fig_rmsd <- ggplot(df_rmsd, aes(x = time, y = mean_A, color = ph_label, fill = ph_label)) +
  geom_ribbon(aes(ymin = pmax(0, mean_A - sem_A), ymax = mean_A + sem_A), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  annotate("rect", xmin = 150, xmax = 200, ymin = -Inf, ymax = Inf, fill = "gray93", alpha = 0.5) +
  geom_vline(xintercept = 150, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  annotate("text", x = 175, y = Inf, label = "Ventana\nMM-GBSA", vjust = 1.3, size = 2.9, color = "gray45") +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_color_manual(values = PH_COLS) + scale_fill_manual(values = PH_COLS) +
  scale_x_continuous(breaks = seq(0, 200, 50), expand = c(0.01, 0)) +
  labs(x = "Tiempo (ns)", y = "RMSD (A)", color = "pH", fill = "pH",
       title = "Convergencia estructural del complejo semaglutida-GLP-1R",
       subtitle = "Suavizado 4 ns | Banda = +/- SEM (n=3) | Zona gris = MM-GBSA 150-200 ns") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "gray95", color = "gray70"),
        strip.text = element_text(size = 10.5, face = "bold"),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, color = "gray40"),
        panel.grid.major.y = element_line(color = "gray92", linewidth = 0.35)) +
  guides(color = guide_legend(nrow = 1))
ggsave(out("fig_rmsd.png"), fig_rmsd, width = 7, height = 6.5, dpi = 300)
cat("OK: fig_rmsd.png\n")


# ==================================================================
# FIG: Puentes de hidrogeno
# ==================================================================
time_grid_ps <- seq(0, 200000, by = 500)
all_hb <- list()
for (i in seq_along(phs_dir)) {
  ph <- phs_dir[i]
  mat <- matrix(NA, nrow = length(time_grid_ps), ncol = 3)
  for (r in 1:3) {
    path <- file.path(BASE_DIR, paste0("Replica_", r), ph, "HBond", "hbonds_Rec_Lig.xvg")
    if (file.exists(path)) {
      d <- read_xvg(path)
      if (max(d$x, na.rm = TRUE) < 1000) d$x <- d$x * 1000
      mat[, r] <- approx(d$x, d$y, xout = time_grid_ps, rule = 2)$y
    }
  }
  all_hb[[i]] <- data.frame(time_ns = time_grid_ps / 1000,
                             mean = rowMeans(mat, na.rm = TRUE),
                             sd = apply(mat, 1, sd, na.rm = TRUE),
                             ph_label = PH_LABS[i])
}

smf <- function(x, k = 21) {
  n <- length(x); out <- rep(NA_real_, n); half <- floor(k / 2)
  for (ii in seq_len(n)) { lo <- max(1, ii - half); hi <- min(n, ii + half); out[ii] <- mean(x[lo:hi], na.rm = TRUE) }
  out
}

smooth_hb <- bind_rows(all_hb) %>%
  mutate(ph_label = factor(ph_label, levels = PH_LABS)) %>%
  group_by(ph_label) %>% mutate(mean_s = smf(mean), sd_s = smf(sd)) %>% ungroup()

fig_hb <- ggplot(smooth_hb, aes(x = time_ns, y = mean_s, color = ph_label, fill = ph_label)) +
  geom_ribbon(aes(ymin = pmax(0, mean_s - sd_s), ymax = mean_s + sd_s), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.75) +
  annotate("rect", xmin = 150, xmax = 200, ymin = -Inf, ymax = Inf, fill = "gray95", alpha = 0.5) +
  geom_vline(xintercept = 150, linetype = "dashed", color = "gray45", linewidth = 0.5) +
  annotate("label", x = 198, y = c(9.7, 15.2, 18.1),
           label = c("9.7 +/- 3.0", "15.2 +/- 5.5", "18.1 +/- 4.4"), hjust = 1, size = 3,
           fill = c("#FDEAEA", "#E8F4F8", "#E8F6EE"),
           color = c("#C0392B", "#1A6080", "#1E8449"), label.size = 0.3) +
  scale_color_manual(values = PH_COLS) + scale_fill_manual(values = PH_COLS) +
  scale_x_continuous(breaks = seq(0, 200, 50), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, 5)) +
  labs(x = "Tiempo (ns)", y = "Puentes de hidrogeno", color = "pH", fill = "pH",
       title = "Interacciones receptor-ligando durante la simulacion",
       subtitle = "Etiquetas = promedio +/- DE en ventana 150-200 ns | Zona gris = MM-GBSA") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, color = "gray40", lineheight = 1.3),
        panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4)) +
  guides(color = guide_legend(nrow = 1))
ggsave(out("fig_hbonds.png"), fig_hb, width = 7, height = 5, dpi = 300)
cat("OK: fig_hbonds.png\n")


# ==================================================================
# FIG: RMSF
# ==================================================================
collect_rmsf <- function(rf) {
  all_ph <- list()
  for (i in seq_along(phs_dir)) {
    ph <- phs_dir[i]; all_rep <- list()
    for (r in 1:3) {
      path <- file.path(BASE_DIR, paste0("Replica_", r), ph, "RMSF", rf)
      if (file.exists(path)) all_rep[[r]] <- read_xvg(path)
    }
    if (length(all_rep) == 0) next
    res_grid <- seq(min(sapply(all_rep, function(d) min(d$x))),
                    max(sapply(all_rep, function(d) max(d$x))))
    mat <- matrix(NA, nrow = length(res_grid), ncol = 3)
    for (r in 1:3) if (!is.null(all_rep[[r]]))
      mat[, r] <- approx(all_rep[[r]]$x, all_rep[[r]]$y, xout = res_grid, rule = 2)$y
    all_ph[[i]] <- data.frame(resnum = res_grid, mean = rowMeans(mat, na.rm = TRUE),
                              sd = apply(mat, 1, sd, na.rm = TRUE), ph_label = PH_LABS[i])
  }
  bind_rows(all_ph)
}

df_rmsf <- bind_rows(
  collect_rmsf("rmsf_receptor.xvg") %>% mutate(panel = "Receptor GLP-1R"),
  collect_rmsf("rmsf_ligando.xvg") %>% mutate(panel = "Semaglutida")) %>%
  mutate(ph_label = factor(ph_label, levels = PH_LABS),
         panel = factor(panel, levels = c("Receptor GLP-1R", "Semaglutida")),
         mean_A = mean * 10, sd_A = sd * 10)

fig_rmsf <- ggplot(df_rmsf, aes(x = resnum, y = mean_A, color = ph_label, fill = ph_label)) +
  geom_ribbon(aes(ymin = pmax(0, mean_A - sd_A), ymax = mean_A + sd_A), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~panel, ncol = 1, scales = "free") +
  scale_color_manual(values = PH_COLS) + scale_fill_manual(values = PH_COLS) +
  labs(x = "Numero de residuo", y = "RMSF (A)", color = "pH", fill = "pH",
       title = "Flexibilidad local del receptor GLP-1R y semaglutida",
       subtitle = "Media +/- DE, n = 3 replicas") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "gray95", color = "gray70"),
        strip.text = element_text(size = 10.5, face = "bold"),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, color = "gray40"),
        panel.grid.major.y = element_line(color = "gray92", linewidth = 0.4)) +
  guides(color = guide_legend(nrow = 1))
ggsave(out("fig_rmsf.png"), fig_rmsf, width = 7.5, height = 6.5, dpi = 300)
cat("OK: fig_rmsf.png\n")


# ==================================================================
# FIG: Heatmap hotspots VdW
# ==================================================================
decomp_path <- file.path(BASE_DIR, "../Resultados_Hotspots/decomposicion_completa_receptor.csv")
if (file.exists(decomp_path)) {
  decomp <- read.csv(decomp_path)
  resnum_label <- c("204" = "MET204", "244" = "LEU244", "183" = "LEU183", "180" = "HIS180",
                     "166" = "LEU166", "233" = "MET233", "387" = "GLU387", "363" = "HIS363", "359" = "LEU359")
  key_resnums <- as.integer(names(resnum_label))

  key_data <- decomp %>% filter(resnum %in% key_resnums) %>%
    group_by(resnum, ph_label) %>%
    summarise(vdw = mean(vdw_mean, na.rm = TRUE), .groups = "drop") %>%
    mutate(display_label = resnum_label[as.character(resnum)],
           is_ph_sensitive = resnum %in% c(180, 363, 387, 359),
           ph_label = factor(ph_label, levels = c("pH 5.0", "pH 7.4", "pH 8.0"),
                             labels = c("pH 5.0\n(endosomal)", "pH 7.4\n(sistemico)", "pH 8.0\n(NAC/SNAC)")))

  ord <- key_data %>% filter(grepl("sist", as.character(ph_label))) %>%
    arrange(vdw) %>% pull(display_label)
  key_data$display_label <- factor(key_data$display_label, levels = rev(unique(ord)))

  fig_heat <- ggplot(key_data, aes(x = ph_label, y = display_label, fill = vdw)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = round(vdw, 1)), size = 3.2, color = "white", fontface = "bold") +
    geom_tile(data = filter(key_data, is_ph_sensitive), color = "#E74C3C", fill = NA, linewidth = 1.3) +
    scale_fill_gradient2(low = "#1A5276", mid = "#2E86C1", high = "#AED6F1", midpoint = -6,
                         name = "VdW\n(kcal/mol)") +
    labs(x = "pH", y = NULL,
         title = "Residuos clave del sitio de union semaglutida-GLP-1R",
         subtitle = "Contribucion VdW | Borde rojo = pH-sensible") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8.5, color = "gray40"),
          legend.position = "right")
  ggsave(out("fig_hotspots_heatmap.png"), fig_heat, width = 7, height = 5.5, dpi = 300)
  cat("OK: fig_hotspots_heatmap.png\n")

  # FIG: Lollipop VdW
  his_resnums <- c(99, 171, 173, 180, 212, 363, 374)
  piv_vdw <- decomp %>%
    mutate(display_label = case_when(resnum == 180 ~ "HIS180", resnum == 363 ~ "HIS363", TRUE ~ residue_label),
           resname_clean = case_when(resnum %in% his_resnums ~ "HIS", TRUE ~ resname)) %>%
    group_by(display_label, resname_clean, resnum, ph_label) %>%
    summarise(vdw_mean = mean(vdw_mean, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = ph_label, values_from = vdw_mean, values_fn = mean) %>%
    rename(pH50 = `pH 5.0`, pH74 = `pH 7.4`, pH80 = `pH 8.0`)

  vdw_hot <- piv_vdw %>%
    filter(!resname_clean %in% c("ARG", "ASP", "GLU", "LYS", "TYR")) %>%
    filter(pH50 <= -2 | pH74 <= -2 | pH80 <= -2) %>%
    mutate(delta_pH = abs(coalesce(pH74, 0) - coalesce(pH50, 0)),
           ph_sensible = delta_pH > 1.0 | resnum %in% his_resnums,
           mean_vdw = rowMeans(cbind(pH50, pH74, pH80), na.rm = TRUE)) %>%
    arrange(mean_vdw) %>%
    distinct(display_label, .keep_all = TRUE) %>%
    head(22) %>%
    mutate(label = factor(display_label, levels = display_label))

  fig_lol <- ggplot(vdw_hot, aes(x = label, y = pH74, color = ph_sensible)) +
    geom_segment(aes(xend = label, y = 0, yend = pH74), linewidth = 0.7) +
    geom_point(size = 4) +
    geom_point(aes(y = pH50), shape = 1, size = 2.8, alpha = 0.7) +
    geom_point(aes(y = pH80), shape = 2, size = 2.8, alpha = 0.7) +
    geom_hline(yintercept = 0, color = "gray55", linewidth = 0.4) +
    geom_hline(yintercept = -2, color = "gray80", linewidth = 0.3, linetype = "dashed") +
    coord_flip() +
    scale_color_manual(values = c("FALSE" = "#5B8DB8", "TRUE" = "#E74C3C"),
                       labels = c("FALSE" = "Estable", "TRUE" = "pH-sensible")) +
    labs(x = NULL, y = "VdW a pH 7.4 (kcal/mol)", color = "Clasificacion",
         title = "Nucleo hidrofobico de union del GLP-1R",
         subtitle = "Relleno = pH 7.4 | Circulo = pH 5.0 | Triangulo = pH 8.0") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9.5, color = "gray40"),
          panel.grid.major.x = element_line(color = "gray92", linewidth = 0.4))
  ggsave(out("fig_hotspots_lollipop.png"), fig_lol, width = 7.5, height = 6.5, dpi = 300)
  cat("OK: fig_hotspots_lollipop.png\n")
} else {
  cat("[SKIP] Hotspot CSV no encontrado — corre 05_analisis_hotspots.py primero\n")
}


# ==================================================================
# FIG: RMSD vs cristal (PDB 7KI0)
# ==================================================================
csv_cristal <- file.path(OUTDIR, "rmsd_cristal_data.csv")
if (file.exists(csv_cristal)) {
  df_c <- read.csv(csv_cristal, colClasses = c(ph_label = "character"))
  df_c$ph_label <- sprintf("%.1f", as.numeric(df_c$ph_label))
  df_c$ph_label <- factor(df_c$ph_label, levels = c("5.0", "7.4", "8.0"))
  df_c$replica  <- factor(df_c$replica)

  colores_ph <- c("5.0" = "#E41A1C", "7.4" = "#377EB8", "8.0" = "#4DAF4A")
  umbral <- 1.5

  p_barras <- ggplot(df_c, aes(x = ph_label, y = rmsd_angstrom, fill = ph_label)) +
    stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x),
                 fun.max = function(x) mean(x) + sd(x),
                 geom = "crossbar", width = 0.45, color = "grey30", linewidth = 0.5, fatten = 2) +
    geom_jitter(aes(shape = replica), width = 0.12, size = 3.5, color = "grey20", stroke = 0.8) +
    geom_hline(yintercept = umbral, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    annotate("text", x = 0.55, y = umbral + 0.04,
             label = paste0("Umbral = ", umbral, " A"), hjust = 0, size = 3.2, color = "grey40") +
    scale_fill_manual(values = colores_ph, guide = "none") +
    scale_shape_manual(values = c(16, 17, 15), name = "Replica", labels = c("R1", "R2", "R3")) +
    scale_y_continuous(limits = c(0, max(df_c$rmsd_angstrom) * 1.3), expand = c(0, 0)) +
    labs(title = "Validacion de la pose de union de semaglutida",
         subtitle = "RMSD backbone vs PDB 7KI0 al inicio de produccion",
         x = "pH", y = "RMSD (A)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10, color = "grey40"),
          axis.title = element_text(face = "bold"), legend.position = "right")
  ggsave(out("fig_rmsd_cristal.png"), p_barras, width = 7, height = 5, dpi = 300)

  p_tabla <- ggplot(df_c, aes(x = ph_label, y = replica,
                               label = sprintf("%.2f A", rmsd_angstrom), fill = rmsd_angstrom)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(size = 3.8, fontface = "bold",
              color = ifelse(df_c$rmsd_angstrom > 1.2, "white", "grey20")) +
    scale_fill_gradient(low = "#C6EFCE", high = "#E41A1C", limits = c(0, umbral),
                        name = "RMSD (A)", oob = scales::squish) +
    scale_y_discrete(labels = c("1" = "Replica 1", "2" = "Replica 2", "3" = "Replica 3")) +
    labs(title = "Tabla validacion: RMSD semaglutida vs PDB 7KI0", x = "pH", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank(), legend.position = "right")
  ggsave(out("fig_rmsd_cristal_tabla.png"), p_tabla, width = 6, height = 4, dpi = 300)

  write.csv(df_c %>% select(replica, ph_label, rmsd_angstrom) %>%
              pivot_wider(names_from = ph_label, values_from = rmsd_angstrom, names_prefix = "pH_"),
            file.path(OUTDIR, "Tabla_RMSD_vs_cristal.csv"), row.names = FALSE)
  cat("OK: fig_rmsd_cristal.png + fig_rmsd_cristal_tabla.png\n")
} else {
  cat("[SKIP] rmsd_cristal_data.csv no encontrado — corre 02_analisis_estructural.sh primero\n")
}


# ==================================================================
# CUADROS como PNG (tablas visuales ggplot)
# ==================================================================
draw_table <- function(df, title, subtitle = NULL, header_fill = "#2E74B5",
                       col_widths = NULL, outfile, w = 8, h = NULL, highlight_rows = NULL) {
  nc <- ncol(df); nr <- nrow(df)
  if (is.null(col_widths)) col_widths <- rep(1 / nc, nc)
  col_widths <- col_widths / sum(col_widths)
  x_centers <- cumsum(col_widths) - col_widths / 2
  row_h <- 1 / (nr + 2)
  if (is.null(h)) h <- max(3.5, (nr + 2) * 0.45 + 1.5)

  cells <- data.frame()
  for (j in 1:nc) cells <- rbind(cells, data.frame(x = x_centers[j], y = 1 - row_h * 0.5,
    label = names(df)[j], row_type = "header", col = j, row = 0, stringsAsFactors = FALSE))
  for (i in 1:nr) for (j in 1:nc) cells <- rbind(cells, data.frame(x = x_centers[j],
    y = 1 - row_h * (i + 0.5), label = as.character(df[i, j]), row_type = "data",
    col = j, row = i, stringsAsFactors = FALSE))

  rects <- data.frame()
  for (i in 0:nr) {
    fc <- if (i == 0) header_fill
          else if (!is.null(highlight_rows) && i %in% highlight_rows) "#FFF3CD"
          else if (i %% 2 == 0) "#EBF3FA" else "white"
    for (j in 1:nc) rects <- rbind(rects, data.frame(
      xmin = ifelse(j == 1, 0, cumsum(col_widths)[j - 1]), xmax = cumsum(col_widths)[j],
      ymin = 1 - row_h * (i + 1), ymax = 1 - row_h * i, fill = fc))
  }

  subtitle_lines <- if (!is.null(subtitle)) length(strsplit(subtitle, "\n")[[1]]) else 0
  top_margin <- 0.06 + 0.04 * subtitle_lines

  p <- ggplot() +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              color = "white", linewidth = 0.5) +
    scale_fill_identity() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 1 - row_h * (nr + 1), ymax = 1,
             fill = NA, color = "gray60", linewidth = 0.6) +
    geom_text(data = filter(cells, row_type == "header"), aes(x = x, y = y, label = label),
              color = "white", fontface = "bold", size = 3.4, hjust = 0.5, vjust = 0.5, lineheight = 1.1) +
    geom_text(data = filter(cells, row_type == "data"), aes(x = x, y = y, label = label),
              color = "gray10", size = 3.2, hjust = 0.5, vjust = 0.5, lineheight = 1.1) +
    annotate("text", x = 0, y = 1 + 0.055, label = title, hjust = 0, vjust = 0,
             size = 4.2, fontface = "bold", color = "gray10")

  if (!is.null(subtitle))
    p <- p + annotate("text", x = 0, y = 1 + 0.015, label = subtitle, hjust = 0, vjust = 0,
                       size = 2.9, color = "gray50")

  p <- p + coord_cartesian(xlim = c(-0.01, 1.01),
                            ylim = c(1 - row_h * (nr + 1) - 0.02, 1 + top_margin), clip = "off") +
    theme_void() + theme(plot.margin = margin(12, 10, 10, 10))
  ggsave(outfile, p, width = w, height = h, dpi = 300, bg = "white")
  cat("OK:", outfile, "\n")
}

# Cuadro 1 — Protocolo
d1 <- data.frame(
  "Parametro" = c("dt", "nsteps", "Duracion/replica", "Frames totales", "Termostato",
    "Temperatura", "Barostato", "Presion", "Electrostatica",
    "Radio de corte", "Campo de fuerzas", "Membrana", "Agua/iones",
    "Replicas por pH", "Total simulado"),
  "Valor" = c("0.002 ps (2 fs)", "100,000,000", "200 ns", "4,000 (1 frame/50 ps)", "v-rescale",
    "310.15 K (37 C)", "C-rescale semi-isotropico", "1.0 bar", "PME",
    "1.2 nm", "CHARMM36m", "POPC (140 lipidos)", "TIP3P / NaCl 150 mM",
    "3 (independientes)", "1800 ns (3x3x200 ns)"),
  "Descripcion" = c("Integracion numerica", "200 ns produccion", "Por replica", "Resolucion 50 ps",
    "3 grupos independientes", "Fisiologico", "Bicapa lipidica", "Isobarico",
    "Alta precision", "Estandar CHARMM36m", "Proteina-membrana", "1-palmitoil-2-oleoil-PC",
    "Concentracion fisiologica", "Semillas distintas", "9 sistemas"),
  check.names = FALSE)

draw_table(d1, title = "Cuadro 1. Parametros del protocolo de simulacion (GROMACS 2025)",
  subtitle = "Confirmados de step7_production.mdp (CHARMM-GUI)",
  outfile = out("cuadro_1_protocolo.png"), col_widths = c(0.28, 0.30, 0.42), w = 10, h = 8)

# Cuadro 2 — Energia libre
d2 <- data.frame(
  "Condicion" = c("Endosomal tardio", "Sistemico/plasmatico", "Gastrico/NAC-SNAC"),
  "pH" = c("5.0", "7.4", "8.0"),
  "DG binding (kcal/mol)" = c("-105.60", "-131.26", "-151.35"),
  "DE" = c("33.41", "42.70", "27.18"),
  "SEM" = c("19.29", "24.65", "15.70"),
  "Puentes H (media +/- DE)" = c("9.7 +/- 3.0", "15.2 +/- 5.5", "18.1 +/- 4.4"),
  check.names = FALSE)
draw_table(d2, title = "Cuadro 2. Energia libre de union MM-GBSA semaglutida-GLP-1R",
  subtitle = "igb=5, NaCl 150 mM, frames 150-200 ns, n=3 replicas",
  outfile = out("cuadro_2_mmgbsa.png"), col_widths = c(0.22, 0.09, 0.22, 0.10, 0.10, 0.27), w = 11, h = 4.5)

# Cuadro 3 — Histidinas
d3 <- data.frame(
  "Residuo" = c("His99", "His171", "His173", "His180", "His212", "His363", "His374"),
  "Localizacion" = c("N-term ECD", "N-term ECD", "N-term ECD", "ECL1 contacto lig.", "ECL2", "TM6", "TM6/TM7"),
  "pH 5.0" = rep("HSP (+1)", 7), "pH 7.4" = rep("HSD (neutro)", 7),
  "pH 8.0" = rep("HSD (neutro)", 7), "pKa" = rep("5.0-7.4", 7),
  check.names = FALSE)
draw_table(d3, title = "Cuadro 3. Histidinas pH-sensibles del GLP-1R (PropKa 3.0)",
  subtitle = "HSP = doblemente protonada (+1) | HSD = neutra | ECL = loop extracelular | TM = helice TM",
  outfile = out("cuadro_3_histidinas.png"), col_widths = c(0.12, 0.28, 0.14, 0.14, 0.14, 0.18), w = 11, h = 5)

# Cuadro 4 — Hotspots
d4 <- data.frame(
  "Residuo" = c("GLU387", "HIS180*", "HIS363*", "LEU359*", "MET204", "MET233", "LEU244", "LEU183", "LEU166"),
  "Localizacion" = c("TM6 (profundo)", "ECL1 contacto lig.", "TM6", "TM5-TM6", "ECL1-TM2", "TM3", "TM3", "TM2", "TM2"),
  "DG pH 5.0" = c("-20.96", "-5.27", "-4.08", "-4.19", "-6.76", "-6.05", "-5.74", "-5.87", "-5.68"),
  "DG pH 7.4" = c("-51.94", "-6.04", "-6.06", "-5.22", "-6.67", "-6.18", "-6.11", "-5.98", "-5.89"),
  "DG pH 8.0" = c("-53.98", "-6.30", "-5.93", "-5.30", "-7.08", "-6.80", "-6.11", "-6.00", "-6.21"),
  "pH-sensible" = c("Si", "Si", "Si", "Si", "No", "No", "No", "No", "No"),
  check.names = FALSE)
draw_table(d4, title = "Cuadro 4. Hotspots farmacologicos del GLP-1R (kcal/mol)",
  subtitle = paste0("* pH-sensible (|DDG| > 1.0 entre pH 5.0 y 7.4) | ",
    "GLU387: DG total pairwise (VdW+elec+GB); VdW aislado = -3.85/-3.16/-2.19\n",
    "22 hotspots totales, 11 pH-sensibles | n = 9 trayectorias"),
  outfile = out("cuadro_4_hotspots.png"),
  col_widths = c(0.15, 0.22, 0.13, 0.13, 0.13, 0.13),
  highlight_rows = c(1, 2, 3, 4), w = 12, h = 6.5)

cat("\n=== Figuras guardadas en:", OUTDIR, "===\n")
cat("SIGUIENTE PASO: python3 08_tablas.py\n")
