# =============================================================================
# LiDAR Point Cloud Analysis using lidR v8.1
# =============================================================================
# CHANGES from v8:
#   [v8.1-FIX1] make_biomass_raster() rewritten — rasterize from crown
#               polygons (tree_met) instead of raw point cloud.
#               Eliminates std::bad_alloc on large files (>1GB LAS).
#   [v8.1-FIX2] terra::vect() coercion fixed — explicit sf→SpatVector path.
#   [v8.1-FIX3] filter_duplicates() added after readLAS (Step 1).
#   [v8.1-FIX4] IoU length mismatch guard added in compute_seg_metrics().
#   [v8.1-FIX5] gc() calls added at key steps to release RAM proactively.
#
# Segmentation References:
#   [A] Dalponte & Coomes (2016) Methods Ecol Evol 7:1427-1436
#   [B] Silva et al. (2016) Methods Ecol Evol - silva2016()
# Biomass References:
#   [1] Soraya et al. (2025) ecoeet 26(2):178-192
#   [2] Komiyama et al. (2005) Forest Ecology & Management
#   [3] Chave et al. (2014) Global Change Biology
#   [4] IPCC (2006) AFOLU Guidelines CF=0.47
# =============================================================================

# -- 0. Packages (auto-install if missing) ------------------------------------
required_pkgs <- c(
  "lidR", "sf", "terra", "raster", "ggplot2", "dplyr",
  "viridis", "tidyr", "tools", "future", "patchwork",
  "furrr", "parallel", "doParallel"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat(sprintf("Installing missing packages: %s\n", paste(missing_pkgs, collapse = ", ")))
  install.packages(missing_pkgs, dependencies = TRUE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))
cat("All packages loaded OK.\n\n")

# -- macOS X11 Graphics Fix --
if (Sys.info()["sysname"] == "Darwin") {
  options(device = "quartz")
  cat("[macOS] Using quartz graphics backend (X11 not required)\n\n")
}

# =============================================================================
# PARALLEL PROCESSING SETUP
# =============================================================================
max_available_cores <- parallel::detectCores()
system_reserve      <- 2
recommended_cores   <- max(1, max_available_cores - system_reserve)

cat(sprintf("System cores available: %d\n", max_available_cores))
cat(sprintf("Recommended cores to use: %d\n", recommended_cores))
cat("Note: Adjust n_cores in USER SETTINGS if needed\n\n")

chunk_size_m   <- 50
chunk_buffer_m <- 15
n_cores        <- as.integer(min(8, max_available_cores - 1))

future::plan(future::multisession, workers = n_cores)
cat(sprintf("Future plan set to: multisession with %d workers\n", n_cores))

if (Sys.info()["sysname"] == "Windows") {
  cat("Parallel backend: future::multisession (Windows compatible)\n")
} else {
  cat(sprintf("Parallel backend: future::multisession (%d cores)\n", n_cores))
}

if (!is.na(Sys.getenv("OMP_NUM_THREADS"))) {
  Sys.setenv(OMP_NUM_THREADS = n_cores)
  cat(sprintf("OpenMP threads set to: %d\n", n_cores))
}

cat("\n=== PARALLELISM DIAGNOSTIC ===\n")
cat(sprintf("  CPU cores (logical)  : %d\n", parallel::detectCores()))
cat(sprintf("  n_cores configured   : %d\n", n_cores))
cat(sprintf("  lidR OpenMP threads  : %d\n", lidR::get_lidr_threads()))
cat(sprintf("  future plan          : %s\n", class(future::plan())[1]))
cat(sprintf("  OMP_NUM_THREADS      : %s\n", Sys.getenv("OMP_NUM_THREADS")))
cat("==============================\n\n")

# =============================================================================
# USER SETTINGS
# =============================================================================
input_las    <- "./data/1E_05_v1.las"
output_dir   <- "./output_v8/1E_05_test1"
site_name    <- "1E_05"

target_epsg  <- 32749
wood_density <- 0.57
a_feld       <- 0.5279
b_feld       <- 0.5782
CF           <- 0.47

min_pts_export <- 5
match_dist_m   <- 2.0

cat(sprintf("USER SETTINGS: Using %d cores for parallel processing\n", n_cores))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

safe_write_shp <- function(obj, path, force_2d = TRUE) {
  if (force_2d) obj <- sf::st_zm(obj, drop = TRUE, what = "ZM")
  exts <- c(".shp",".shx",".dbf",".prj",".cpg",".qpj",".sbn",".sbx")
  base <- tools::file_path_sans_ext(path)
  for (e in exts) { f <- paste0(base, e); if (file.exists(f)) file.remove(f) }
  suppressWarnings(sf::st_write(obj, path, delete_dsn = FALSE, quiet = TRUE))
  cat(sprintf("    Saved SHP : %s\n", basename(path)))
}

safe_write_gpkg <- function(obj, path) {
  if (file.exists(path)) file.remove(path)
  sf::st_write(obj, path, delete_dsn = FALSE, quiet = TRUE)
  cat(sprintf("    Saved GPKG: %s\n", basename(path)))
}

assign_crs_if_missing <- function(las, epsg_code) {
  if (is.na(sf::st_crs(las))) {
    cat(sprintf("  WARNING: CRS is NA - assigning EPSG:%d\n", epsg_code))
    lidR::projection(las) <- sf::st_crs(epsg_code)
  } else {
    cat(sprintf("  CRS OK: %s\n", sf::st_crs(las)$input))
  }
  return(las)
}

# ---------------------------------------------------------------------------
# [v8.1-FIX1] make_biomass_raster() — polygon-based rasterization
# Rasterizes from crown polygon centroids (tiny sf, ~n_trees rows) instead
# of raw point cloud (125M rows). Eliminates std::bad_alloc entirely.
# ---------------------------------------------------------------------------
make_biomass_raster <- function(tree_met_sf, chm_template, attribute_name,
                                normalize_to_ha = FALSE) {
  centroids <- suppressWarnings(sf::st_centroid(tree_met_sf))
  attr_vals <- sf::st_drop_geometry(centroids)[[attribute_name]]
  coords    <- sf::st_coordinates(centroids)

  df_pts <- data.frame(X = coords[,1], Y = coords[,2], value = attr_vals)
  df_pts <- df_pts[!is.na(df_pts$value), ]

  pts_vect <- terra::vect(
    as.matrix(df_pts[, c("X","Y")]),
    type = "points",
    atts = data.frame(value = df_pts$value),
    crs  = terra::crs(chm_template))

  # Rasterize: sum all tree values falling in each pixel
  r_sum <- terra::rasterize(pts_vect, chm_template, field = "value", fun = "sum")

  if (normalize_to_ha) {
    # pixel area in m2 → convert sum(kg/pixel) to Mg/ha
    pix_res      <- terra::res(chm_template)          # e.g. c(0.25, 0.25)
    pix_area_m2  <- pix_res[1] * pix_res[2]           # e.g. 0.0625 m2
    # kg/pixel ÷ 1000 → Mg/pixel; ÷ pix_area_m2 × 10000 → Mg/ha
    r_out <- r_sum / 1000 / pix_area_m2 * 10000
  } else {
    r_out <- r_sum
  }

  names(r_out) <- attribute_name
  return(r_out)
}


# ---------------------------------------------------------------------------
# VALIDATION HELPER: compute_seg_metrics()
# [v8.1-FIX4] IoU length mismatch guard added
# Reference = Silva2016; Test = Dalponte2016
# ---------------------------------------------------------------------------
compute_seg_metrics <- function(ref_met, test_met,
                                 match_dist_m = 2.0,
                                 algo_name    = "Test") {
  cat(sprintf("  [Validation] %s vs Reference (Silva2016)...\n", algo_name))

  ref_pts  <- sf::st_centroid(ref_met)
  test_pts <- sf::st_centroid(test_met)

  nn_idx  <- sf::st_nearest_feature(ref_pts, test_pts)
  nn_dist <- as.numeric(sf::st_distance(
    ref_pts, test_pts[nn_idx, ], by_element = TRUE))

  matched  <- nn_dist <= match_dist_m
  n_ref    <- nrow(ref_met)
  n_test   <- nrow(test_met)
  TP       <- sum(matched)
  FP       <- n_test - TP
  FN       <- n_ref  - TP

  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall    <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  f1        <- ifelse((precision + recall) > 0,
                      2 * precision * recall / (precision + recall), 0)

  ref_match  <- ref_met[matched, ]
  test_match <- test_met[nn_idx[matched], ]

  # [v8.1-FIX4] Guard against geometry length mismatch in st_intersection
  iou_vals <- tryCatch({
    n_pairs    <- nrow(ref_match)
    inter_area <- numeric(n_pairs)
    ref_area   <- as.numeric(sf::st_area(ref_match))
    test_area  <- as.numeric(sf::st_area(test_match))
    for (k in seq_len(n_pairs)) {
      inter_k <- tryCatch(
        sf::st_intersection(
          sf::st_geometry(ref_match)[k],
          sf::st_geometry(test_match)[k]),
        error = function(e) NULL)
      inter_area[k] <- if (!is.null(inter_k) && length(inter_k) > 0)
        as.numeric(sf::st_area(inter_k)) else 0
    }
    union_area <- ref_area + test_area - inter_area
    ifelse(union_area > 0, inter_area / union_area, NA_real_)
  }, error = function(e) rep(NA_real_, TP))

  mean_iou <- mean(iou_vals, na.rm = TRUE)

  h_ref   <- sf::st_drop_geometry(ref_match)$z_max
  h_test  <- sf::st_drop_geometry(test_match)$z_max
  rmse_h  <- sqrt(mean((h_ref  - h_test)^2,  na.rm = TRUE))
  mae_h   <- mean(abs(h_ref   - h_test),      na.rm = TRUE)

  ca_ref  <- sf::st_drop_geometry(ref_match)$Crown_area_m2
  ca_test <- sf::st_drop_geometry(test_match)$Crown_area_m2
  rmse_ca <- sqrt(mean((ca_ref - ca_test)^2,  na.rm = TRUE))
  mae_ca  <- mean(abs(ca_ref  - ca_test),      na.rm = TRUE)

  result <- data.frame(
    Algorithm  = algo_name,
    N_ref      = n_ref,   N_detected = n_test,
    TP = TP,   FP = FP,   FN = FN,
    Precision  = round(precision, 4),
    Recall     = round(recall,    4),
    F1_Score   = round(f1,        4),
    Mean_IoU   = round(mean_iou,  4),
    RMSE_H_m   = round(rmse_h,    4),
    MAE_H_m    = round(mae_h,     4),
    RMSE_CA_m2 = round(rmse_ca,   4),
    MAE_CA_m2  = round(mae_ca,    4),
    stringsAsFactors = FALSE
  )

  cat(sprintf("    N_ref=%d | N_det=%d | TP=%d FP=%d FN=%d\n",
              n_ref, n_test, TP, FP, FN))
  cat(sprintf("    Precision=%.3f | Recall=%.3f | F1=%.3f | IoU=%.3f\n",
              precision, recall, f1, mean_iou))
  cat(sprintf("    RMSE_H=%.3fm | MAE_H=%.3fm | RMSE_CA=%.2fm2 | MAE_CA=%.2fm2\n",
              rmse_h, mae_h, rmse_ca, mae_ca))

  list(metrics      = result,
       matched_ref  = ref_match,
       matched_test = test_match,
       iou_vals     = iou_vals,
       nn_dist      = nn_dist,
       matched_flag = matched)
}

# =============================================================================
# SETUP DIRECTORIES
# =============================================================================
dir_rasters   <- file.path(output_dir, "rasters")
dir_plots     <- file.path(output_dir, "plots")
dir_vectors   <- file.path(output_dir, "vectors")
dir_trees     <- file.path(output_dir, "individual_trees")
dir_trees_las <- file.path(dir_trees,  "las")
dir_trees_csv <- file.path(dir_trees,  "attributes")
dir_valid     <- file.path(output_dir, "validation")

for (d in c(dir_rasters, dir_plots, dir_vectors,
            dir_trees, dir_trees_las, dir_trees_csv, dir_valid)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("=== LiDAR Analysis Pipeline v8.1 ===\n")
cat("Input :", input_las,  "\n")
cat("Output:", output_dir, "\n\n")

# =============================================================================
# STEP 1: Read LAS
# =============================================================================
cat("[Step 1] Reading LAS file...\n")

file_size_gb <- file.info(input_las)$size / (1024^3)
cat(sprintf("  File size: %.2f GB\n", file_size_gb))

ctg_peek  <- readLAScatalog(input_las)
z_min_raw <- floor(  min(ctg_peek@data$Min.Z, na.rm = TRUE))
z_max_raw <- ceiling(max(ctg_peek@data$Max.Z, na.rm = TRUE))
cat(sprintf("  Points (header) : %s\n",
            format(sum(ctg_peek@data$Number.of.point.records), big.mark=",")))
cat(sprintf("  Z range (header): %.2f to %.2f m\n", z_min_raw, z_max_raw))
cat(sprintf("  LAS format      : %s\n", ctg_peek@data$Point.Data.Format.ID))

las_filter_str <- paste0(
  "-drop_z_below ", z_min_raw - 5,
  " -drop_z_above ", z_max_raw + 5)

options(lidR.progress = FALSE)
cat("  Reading LAS (progress disabled)...\n")
las_raw <- readLAS(files = input_las, select = "xyzrnic",
                   filter = las_filter_str)
options(lidR.progress = TRUE)

# [v8.1-FIX3] Remove duplicate points flagged by las_check
n_before_dedup <- nrow(las_raw@data)
las_raw        <- filter_duplicates(las_raw)
n_after_dedup  <- nrow(las_raw@data)
n_dupes        <- n_before_dedup - n_after_dedup
if (n_dupes > 0) {
  cat(sprintf("  [Dedup] Removed %s duplicate points (%s remaining)\n",
              format(n_dupes, big.mark=","),
              format(n_after_dedup, big.mark=",")))
}
gc()

print(las_raw)
las_check(las_raw)
cat("\n--- Point Cloud Summary ---\n")
cat("Points  :", format(nrow(las_raw@data), big.mark=","), "\n")
cat("Z range :", round(min(las_raw$Z),2), "to", round(max(las_raw$Z),2), "m\n")
las_raw <- assign_crs_if_missing(las_raw, target_epsg)
cat("\n")

# =============================================================================
# STEP 2: Ground Classification
# =============================================================================
cat("[Step 2] Ground Classification...\n")
has_ground <- any(las_raw$Classification == 2L)
if (!has_ground) {
  cat("  No ground points. Running PMF...\n")
  las_gnd <- classify_ground(las_raw, algorithm = pmf(ws = 5, th = 3))
} else {
  cat("  Ground already classified. Skipping.\n")
  las_gnd <- las_raw
}
rm(las_raw); gc()   # [v8.1-FIX5] free raw LAS from RAM
n_gnd <- sum(las_gnd$Classification == 2L, na.rm = TRUE)
cat(sprintf("  Ground points: %s (%.1f%%)\n",
            format(n_gnd, big.mark=","),
            n_gnd / nrow(las_gnd@data) * 100))

# =============================================================================
# STEP 3: DTM
# =============================================================================
cat("[Step 3] Generating DTM...\n")
dtm <- suppressWarnings(
  rasterize_terrain(las = las_gnd, res = 0.25, algorithm = tin()))
terra::writeRaster(dtm, file.path(dir_rasters,"DTM_0.25m.tif"), overwrite=TRUE)
cat("  Saved: DTM_0.25m.tif\n")

# =============================================================================
# STEP 4: Height Normalization
# =============================================================================
cat("[Step 4] Normalizing heights...\n")
las_norm <- suppressWarnings(
  normalize_height(las = las_gnd, algorithm = tin()))
rm(las_gnd); gc()   # [v8.1-FIX5] free ground LAS
las_norm <- filter_poi(las_norm, Z >= 0 & Z <= 50)
cat(sprintf("  Z after norm: %.2f to %.2f m\n", min(las_norm$Z), max(las_norm$Z)))

# =============================================================================
# STEP 5: CHM
# =============================================================================
cat("\n[Step 5] Generating CHMs...\n")
chm_p2r     <- rasterize_canopy(las = las_norm, res = 0.25,
                                 algorithm = p2r(subcircle = 0.15))
chm_pitfree <- rasterize_canopy(las = las_norm, res = 0.25,
                                 algorithm = pitfree(
                                   thresholds = c(0,2,5,10,15),
                                   max_edge   = c(0,1.5)))
chm_smooth  <- terra::focal(chm_pitfree, w = matrix(1,3,3),
                             fun = median, na.rm = TRUE)

terra::writeRaster(chm_p2r,     file.path(dir_rasters,"CHM_p2r_0.25m.tif"),     overwrite=TRUE)
terra::writeRaster(chm_pitfree, file.path(dir_rasters,"CHM_pitfree_0.25m.tif"), overwrite=TRUE)
terra::writeRaster(chm_smooth,  file.path(dir_rasters,"CHM_smooth_0.25m.tif"),  overwrite=TRUE)
cat("  Saved: CHM p2r / pitfree / smooth\n\n")

# =============================================================================
# STEP 6: Individual Tree Detection
# =============================================================================
cat("[Step 6] Individual Tree Detection...\n")

ttops_fixed_2m <- locate_trees(chm_smooth, algorithm = lmf(ws = 2))
ttops_fixed_3m <- locate_trees(chm_smooth, algorithm = lmf(ws = 3))
f_vwf          <- function(x) { x * 0.1 + 1 }
ttops_variable <- locate_trees(chm_smooth, algorithm = lmf(ws = f_vwf, hmin = 2))

cat(sprintf("  Fixed 2m    : %d trees (raw)\n", nrow(ttops_fixed_2m)))
cat(sprintf("  Fixed 3m    : %d trees (raw)\n", nrow(ttops_fixed_3m)))
cat(sprintf("  Variable WS : %d trees (raw)\n", nrow(ttops_variable)))

# ---------------------------------------------------------------------------
# CROWN RADIUS PRE-FILTER
# Reference: Jucker et al. (2017) Global Ecol Biogeogr 26:1261-1275
# ---------------------------------------------------------------------------
MIN_CROWN_RADIUS_M     <- 1
MIN_CROWN_AREA_M2      <- pi * MIN_CROWN_RADIUS_M^2
MIN_H_FOR_CROWN_RADIUS <- (MIN_CROWN_RADIUS_M / 0.557) ^ (1 / 0.809)

n_before       <- nrow(ttops_variable)
ttops_variable <- ttops_variable %>% dplyr::filter(Z >= MIN_H_FOR_CROWN_RADIUS)
n_after        <- nrow(ttops_variable)

cat(sprintf("\n  [Crown Radius Pre-filter]\n"))
cat(sprintf("    Min crown radius target : %.1f m\n",  MIN_CROWN_RADIUS_M))
cat(sprintf("    Min crown area target   : %.2f m2\n", MIN_CROWN_AREA_M2))
cat(sprintf("    Min apex height (Jucker): %.2f m\n",  MIN_H_FOR_CROWN_RADIUS))
cat(sprintf("    Treetops before filter  : %d\n",      n_before))
cat(sprintf("    Treetops after filter   : %d\n",      n_after))
cat(sprintf("    Removed (too small/low) : %d\n\n",    n_before - n_after))

if (n_after == 0)
  stop("ERROR: No treetops remaining after crown radius pre-filter.")

# =============================================================================
# STEP 7: Individual Tree Segmentation — TWO ALGORITHMS
# =============================================================================
cat("[Step 7] Individual Tree Segmentation - 2 Algorithms...\n")

# --- 7A: Dalponte2016 (reference) ---
cat("  [7A] Dalponte2016 (reference)...\n")
las_seg_dalp <- segment_trees(
  las = las_norm,
  algorithm = dalponte2016(
    chm      = chm_smooth,
    treetops = ttops_variable,
    th_tree  = 2.5,
    th_seed  = 0.45,
    th_cr    = 0.5,
    max_cr   = 15))
n_dalp_raw <- length(unique(na.omit(las_seg_dalp$treeID)))
if (n_dalp_raw == 0) warning("WARNING: Dalponte2016 produced no valid trees")
cat(sprintf("    Trees segmented (raw): %d\n", n_dalp_raw))

# --- 7B: Silva2016 (primary) ---
cat("  [7B] Silva2016 (primary)...\n")
las_seg_silv <- tryCatch({
  segment_trees(
    las = las_norm,
    algorithm = silva2016(
      chm           = chm_smooth,
      treetops      = ttops_variable,
      max_cr_factor = 0.45,
      exclusion     = 0.3))
}, error = function(e) {
  cat(sprintf("    WARNING: silva2016 failed: %s\n  -> Using Dalponte fallback.\n",
              conditionMessage(e)))
  las_seg_dalp
})
n_silv_raw <- length(unique(na.omit(las_seg_silv$treeID)))
if (n_silv_raw == 0) warning("WARNING: Silva2016 produced no valid trees")
cat(sprintf("    Trees segmented (raw): %d\n\n", n_silv_raw))

# ---------------------------------------------------------------------------
# POST-SEGMENTATION FILTER
# ---------------------------------------------------------------------------
MIN_PTS_SEGMENT <- 20L
MIN_PTS_PER_M2  <- 5.0

filter_segments <- function(las_s, label) {
  cat(sprintf("  [Filter: %s]\n", label))

  cm <- crown_metrics(
    las  = las_s,
    func = ~list(z_max = max(Z), n_points = length(Z)),
    geom = "convex") %>%
    dplyr::mutate(
      Crown_area_m2  = as.numeric(sf::st_area(geometry)),
      Crown_radius_m = sqrt(Crown_area_m2 / pi),
      pt_density_m2  = n_points / pmax(Crown_area_m2, 0.01))

  n_total   <- nrow(cm)
  valid_ids <- cm %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(
      Crown_area_m2 >= MIN_CROWN_AREA_M2,
      n_points      >= MIN_PTS_SEGMENT,
      pt_density_m2 >= MIN_PTS_PER_M2) %>%
    dplyr::pull(treeID)

  n_valid   <- length(valid_ids)

  rej_area <- cm %>% sf::st_drop_geometry() %>%
    dplyr::filter(Crown_area_m2 < MIN_CROWN_AREA_M2) %>% nrow()
  rej_pts  <- cm %>% sf::st_drop_geometry() %>%
    dplyr::filter(Crown_area_m2 >= MIN_CROWN_AREA_M2,
                  n_points < MIN_PTS_SEGMENT) %>% nrow()
  rej_dens <- cm %>% sf::st_drop_geometry() %>%
    dplyr::filter(Crown_area_m2 >= MIN_CROWN_AREA_M2,
                  n_points >= MIN_PTS_SEGMENT,
                  pt_density_m2 < MIN_PTS_PER_M2) %>% nrow()

  cat(sprintf("    Segments total          : %d\n",   n_total))
  cat(sprintf("    Rejected — area < %.1fm2: %d  (crown radius < %.1fm)\n",
              MIN_CROWN_AREA_M2, rej_area, MIN_CROWN_RADIUS_M))
  cat(sprintf("    Rejected — pts  < %d    : %d\n",   MIN_PTS_SEGMENT, rej_pts))
  cat(sprintf("    Rejected — dens < %.0f/m2: %d\n",  MIN_PTS_PER_M2,  rej_dens))
  cat(sprintf("    Retained (valid trees)  : %d\n\n", n_valid))

  las_s@data$treeID[!(las_s@data$treeID %in% valid_ids)] <- NA_integer_
  list(las_filtered = las_s, valid_ids = valid_ids, n_valid = n_valid)
}

res_dalp <- filter_segments(las_seg_dalp, "Dalponte2016")
res_silv <- filter_segments(las_seg_silv, "Silva2016")

las_seg_dalp <- res_dalp$las_filtered
las_seg_silv <- res_silv$las_filtered
n_dalp       <- res_dalp$n_valid
n_silv       <- res_silv$n_valid

cat(sprintf("  [Post-filter summary]\n"))
cat(sprintf("    Dalponte2016 : %d -> %d trees (removed %d)\n",
            n_dalp_raw, n_dalp, n_dalp_raw - n_dalp))
cat(sprintf("    Silva2016    : %d -> %d trees (removed %d)\n\n",
            n_silv_raw, n_silv, n_silv_raw - n_silv))

# Primary = Silva2016
las_seg     <- las_seg_silv
n_trees_seg <- n_silv

# =============================================================================
# STEP 8: Tree-Level Crown Metrics (both algorithms)
# =============================================================================
cat("[Step 8] Computing Tree-Level Metrics (2 algorithms)...\n")

compute_crown_metrics <- function(las_s, label) {
  cat(sprintf("  [%s] ", label))
  m <- crown_metrics(
    las  = las_s,
    func = ~list(
      z_max       = max(Z),
      z_mean      = mean(Z),
      z_sd        = sd(Z),
      z_cv        = sd(Z) / mean(Z),
      z_p25       = quantile(Z, 0.25),
      z_p50       = quantile(Z, 0.50),
      z_p75       = quantile(Z, 0.75),
      z_p95       = quantile(Z, 0.95),
      i_mean      = mean(Intensity),
      i_max       = max(Intensity),
      i_sd        = sd(Intensity),
      n_points    = length(Z),
      cover       = sum(Z > 2) / length(Z),
      pct_1st_rtn = sum(ReturnNumber == 1L) / length(ReturnNumber) * 100),
    geom = "convex") %>%
    mutate(
      algorithm     = label,
      DBH_cm        = pmax((z_max / exp(a_feld)) ^ (1 / b_feld), 0.1),
      Crown_area_m2 = as.numeric(sf::st_area(geometry)))
  cat(sprintf("%d trees\n", nrow(m)))
  return(m)
}

tree_met_dalp <- compute_crown_metrics(las_seg_dalp, "Dalponte2016")
tree_met_silv <- compute_crown_metrics(las_seg_silv, "Silva2016")
gc()  # [v8.1-FIX5]

# Primary tree metrics = Silva2016 + biomass
tree_met <- tree_met_silv %>%
  mutate(
    rho          = wood_density,
    AGB_teak_kg  = 0.153  * (DBH_cm ^ 2.382),
    AGB_Kom_kg   = 0.251  * wood_density * (DBH_cm ^ 2.46),
    BGB_Kom_kg   = 0.199  * (wood_density ^ 0.899) * (DBH_cm ^ 2.22),
    TB_Kom_kg    = AGB_Kom_kg + BGB_Kom_kg,
    AGB_Chave_kg = 0.0673 * ((wood_density * (DBH_cm ^ 2) * z_max) ^ 0.976),
    AGC_teak_kg  = AGB_teak_kg  * CF,
    AGC_Kom_kg   = TB_Kom_kg    * CF,
    AGC_Chave_kg = AGB_Chave_kg * CF)

total_area_ha <- as.numeric(
  sf::st_area(sf::st_convex_hull(sf::st_union(tree_met)))) / 10000
tree_met_csv  <- sf::st_drop_geometry(tree_met)
cat(sprintf("\n  Area: %.4f ha | Trees (Silva2016): %d\n\n",
            total_area_ha, nrow(tree_met)))

# =============================================================================
# STEP 9: Biomass & Carbon Estimation (summary)
# =============================================================================
cat("[Step 9] Biomass & Carbon Estimation (Silva2016 primary)...\n")
cat(sprintf("  Trees: %d | Area: %.4f ha\n\n", nrow(tree_met_csv), total_area_ha))

# =============================================================================
# STEP 9B: VALIDATION — Dalponte2016 vs Silva2016 (reference)
# =============================================================================
cat("[Step 9B] Cross-Algorithm Validation (Dalponte2016 vs Silva2016)...\n")
cat(strrep("-", 55), "\n")

if (nrow(tree_met_silv) < 3 || nrow(tree_met_dalp) < 3) {
  cat("  WARNING: Insufficient trees for validation (n < 3)\n")
  val_dalp <- list(
    metrics = data.frame(
      Algorithm="Dalponte2016", N_ref=NA, N_detected=NA,
      TP=NA, FP=NA, FN=NA,
      Precision=NA, Recall=NA, F1_Score=NA, Mean_IoU=NA,
      RMSE_H_m=NA, MAE_H_m=NA, RMSE_CA_m2=NA, MAE_CA_m2=NA),
    iou_vals = NULL)
} else {
  val_dalp <- compute_seg_metrics(tree_met_silv, tree_met_dalp,
                                   match_dist_m, "Dalponte2016")
}

val_master <- val_dalp$metrics

if (!all(is.na(val_master[1,]))) {
  write.csv(val_master,
            file.path(dir_valid, "Validation_CrossAlgorithm.csv"),
            row.names = FALSE)
  cat("\n  Saved: Validation_CrossAlgorithm.csv\n")
}

if (!is.null(val_dalp$iou_vals) && length(val_dalp$iou_vals) > 0) {
  iou_df <- data.frame(Algorithm = "Dalponte2016", IoU = val_dalp$iou_vals)
  write.csv(iou_df,
            file.path(dir_valid, "Validation_IoU_PerTree.csv"),
            row.names = FALSE)
  cat("  Saved: Validation_IoU_PerTree.csv\n\n")
} else {
  iou_df <- NULL
  cat("  Validation metrics not available.\n\n")
}

# =============================================================================
# STEP 10: Export AGB/AGC Raster Maps
# [v8.1-FIX1] Now uses make_biomass_raster(tree_met_sf, ...) — polygon-based
# =============================================================================

# === DEBUG: Check tree_met columns before rasterizing ===
cat("=== tree_met column check ===\n")
cat("Columns:", paste(names(tree_met), collapse=", "), "\n")
cat("Rows:", nrow(tree_met), "\n")
cat("Crown_area_m2 — exists:", "Crown_area_m2" %in% names(tree_met), "\n")
if ("Crown_area_m2" %in% names(tree_met)) {
  cat("Crown_area_m2 summary:\n")
  print(summary(tree_met$Crown_area_m2))
}
cat("AGB_teak_kg summary:\n")
print(summary(tree_met$AGB_teak_kg))
cat("AGB_Chave_kg summary:\n")
print(summary(tree_met$AGB_Chave_kg))
cat("==============================\n")

cat("[Step 10] Exporting Biomass/Carbon Rasters...\n")

if (nrow(tree_met_csv) == 0) {
  cat("  WARNING: No trees found. Skipping biomass raster generation.\n\n")
} else {
cat("  Rasterizing from crown polygon centroids (memory-safe)...\n")

# DBH — keep in cm, no normalization
r_DBH       <- make_biomass_raster(tree_met, chm_smooth, "DBH_cm",
                                    normalize_to_ha = FALSE)

# Biomass & Carbon — normalize to Mg/ha
r_AGB_teak  <- make_biomass_raster(tree_met, chm_smooth, "AGB_teak_kg",
                                    normalize_to_ha = TRUE)
r_AGC_teak  <- make_biomass_raster(tree_met, chm_smooth, "AGC_teak_kg",
                                    normalize_to_ha = TRUE)
r_AGB_Kom   <- make_biomass_raster(tree_met, chm_smooth, "AGB_Kom_kg",
                                    normalize_to_ha = TRUE)
r_AGC_Kom   <- make_biomass_raster(tree_met, chm_smooth, "AGC_Kom_kg",
                                    normalize_to_ha = TRUE)
r_AGB_Chave <- make_biomass_raster(tree_met, chm_smooth, "AGB_Chave_kg",
                                    normalize_to_ha = TRUE)
r_AGC_Chave <- make_biomass_raster(tree_met, chm_smooth, "AGC_Chave_kg",
                                    normalize_to_ha = TRUE)


  biomass_stack <- c(r_AGB_teak, r_AGC_teak, r_AGB_Kom, r_AGC_Kom,
                     r_AGB_Chave, r_AGC_Chave, r_DBH)
  names(biomass_stack) <- c("AGB_teak_kg","AGC_teak_kg","AGB_Kom_kg",
                             "AGC_Kom_kg","AGB_Chave_kg","AGC_Chave_kg","DBH_cm")

  terra::writeRaster(r_AGB_teak,    file.path(dir_rasters,"AGB_Teak_Soraya2025.tif"),  overwrite=TRUE)
  terra::writeRaster(r_AGC_teak,    file.path(dir_rasters,"AGC_Teak_Soraya2025.tif"),  overwrite=TRUE)
  terra::writeRaster(r_AGB_Kom,     file.path(dir_rasters,"AGB_Komiyama2005.tif"),     overwrite=TRUE)
  terra::writeRaster(r_AGC_Kom,     file.path(dir_rasters,"AGC_Komiyama2005.tif"),     overwrite=TRUE)
  terra::writeRaster(r_AGB_Chave,   file.path(dir_rasters,"AGB_Chave2014.tif"),        overwrite=TRUE)
  terra::writeRaster(r_AGC_Chave,   file.path(dir_rasters,"AGC_Chave2014.tif"),        overwrite=TRUE)
  terra::writeRaster(r_DBH,         file.path(dir_rasters,"DBH_estimated.tif"),        overwrite=TRUE)
  terra::writeRaster(biomass_stack, file.path(dir_rasters,"Biomass_Carbon_Stack.tif"), overwrite=TRUE)
  cat("  All rasters saved.\n\n")
}

# =============================================================================
# STEP 11: Full Forest LAS Export
# =============================================================================
cat("[Step 11] Exporting Full Forest LAS...\n")
options(lidR.progress = FALSE)
writeLAS(las_seg,  file.path(dir_vectors,"FullForest_Segmented_Silva2016.las"))
writeLAS(las_norm, file.path(dir_vectors,"FullForest_Normalized.las"))
options(lidR.progress = TRUE)
cat("  Saved: FullForest_Segmented_Silva2016.las\n")
cat("  Saved: FullForest_Normalized.las\n\n")

# =============================================================================
# STEP 12: Tree Crown Polygons & Tree Locations
# =============================================================================
cat("[Step 12] Exporting Crown Polygons & Tree Locations...\n")

safe_write_gpkg(tree_met_silv, file.path(dir_vectors,"TreeCrowns_Silva2016.gpkg"))
safe_write_gpkg(tree_met_dalp, file.path(dir_vectors,"TreeCrowns_Dalponte2016.gpkg"))

crowns_shp <- tree_met %>%
  dplyr::select(treeID,
                H_max_m=z_max, H_mean_m=z_mean, H_sd_m=z_sd,
                H_p25_m=z_p25, H_p50_m=z_p50, H_p75_m=z_p75, H_p95_m=z_p95,
                I_mean=i_mean, I_max=i_max, n_pts=n_points, cover, DBH_cm,
                AGB_teak_kg, AGC_teak_kg, AGB_Kom_kg, AGC_Kom_kg,
                AGB_Chave_kg, AGC_Chave_kg, Crown_m2=Crown_area_m2)
safe_write_shp(crowns_shp, file.path(dir_vectors,"TreeCrowns_Silva2016.shp"))

ttops_attr <- ttops_variable %>%
  dplyr::left_join(
    tree_met_csv %>%
      dplyr::select(treeID, z_max, z_mean, DBH_cm,
                    AGB_teak_kg, AGC_teak_kg, AGB_Kom_kg, AGC_Kom_kg,
                    AGB_Chave_kg, AGC_Chave_kg, Crown_area_m2, n_points),
    by = "treeID")

safe_write_gpkg(ttops_attr, file.path(dir_vectors,"TreeLocations_Points.gpkg"))
ttops_shp <- ttops_attr %>%
  dplyr::select(treeID, Z_apex_m=Z, H_max_m=z_max, H_mean_m=z_mean, DBH_cm,
                AGB_teak_kg, AGC_teak_kg, AGB_Kom_kg, AGC_Kom_kg,
                AGB_Chave_kg, AGC_Chave_kg,
                Crown_m2=Crown_area_m2, n_pts=n_points)
safe_write_shp(ttops_shp,
               file.path(dir_vectors,"TreeLocations_Points.shp"),
               force_2d = TRUE)
cat("  All vector files saved.\n\n")

# =============================================================================
# STEP 13: Individual Tree LAS + CSV Export (Silva2016 primary)
# =============================================================================
cat("[Step 13] Exporting Individual Tree Files...\n")

tree_pt_counts <- las_seg@data %>%
  dplyr::filter(!is.na(treeID)) %>%
  dplyr::group_by(treeID) %>%
  dplyr::summarise(n_pts_check = dplyr::n(), .groups = "drop")

n_skip   <- sum(tree_pt_counts$n_pts_check < min_pts_export)
tree_ids <- sort(unique(na.omit(las_seg$treeID)))
n_ids    <- length(tree_ids)
pb_step  <- max(1, floor(n_ids / 10))

cat(sprintf("  min_pts_export=%d | to export=%d | to skip=%d\n\n",
            min_pts_export, n_trees_seg, n_skip))

all_tree_attrs <- list()
n_exported <- 0L; n_skipped <- 0L; n_error <- 0L

options(lidR.progress = FALSE)
for (i in seq_along(tree_ids)) {
  tid <- tree_ids[i]
  if (i %% pb_step == 0)
    cat(sprintf("    [%3.0f%%] %d/%d (exp=%d skip=%d err=%d)\n",
                i/n_ids*100, i, n_ids, n_exported, n_skipped, n_error))

  las_tree   <- filter_poi(las_seg, treeID == tid)
  n_tree_pts <- nrow(las_tree@data)
  if (n_tree_pts < min_pts_export) { n_skipped <- n_skipped + 1L; next }

  z_vals   <- las_tree$Z; i_vals <- las_tree$Intensity
  h_max    <- max(z_vals,na.rm=TRUE); h_mean <- mean(z_vals,na.rm=TRUE)
  h_sd     <- sd(z_vals,na.rm=TRUE)
  h_p25    <- quantile(z_vals,0.25,na.rm=TRUE)
  h_p50    <- quantile(z_vals,0.50,na.rm=TRUE)
  h_p75    <- quantile(z_vals,0.75,na.rm=TRUE)
  h_p95    <- quantile(z_vals,0.95,na.rm=TRUE)
  i_mean_v <- mean(i_vals,na.rm=TRUE); i_max_v <- max(i_vals,na.rm=TRUE)
  x_ctr    <- mean(las_tree$X,na.rm=TRUE); y_ctr <- mean(las_tree$Y,na.rm=TRUE)

  dbh_cm    <- pmax((h_max / exp(a_feld)) ^ (1 / b_feld), 0.1)
  agb_teak  <- 0.153  * (dbh_cm ^ 2.382)
  agb_kom   <- 0.251  * wood_density * (dbh_cm ^ 2.46)
  bgb_kom   <- 0.199  * (wood_density ^ 0.899) * (dbh_cm ^ 2.22)
  agb_chave <- 0.0673 * ((wood_density * (dbh_cm ^ 2) * h_max) ^ 0.976)
  agc_teak  <- agb_teak * CF
  agc_kom   <- (agb_kom + bgb_kom) * CF
  agc_chave <- agb_chave * CF

  las_filename <- sprintf("tree_%04d.las", tid)
  tryCatch({
    writeLAS(las_tree, file.path(dir_trees_las, las_filename))
    n_exported <- n_exported + 1L
  }, error = function(e) {
    cat(sprintf("    WARNING: tree %d failed: %s\n", tid, conditionMessage(e)))
    n_error <<- n_error + 1L
  })

  attr_row <- data.frame(
    treeID=tid, las_file=las_filename,
    X_centroid=round(x_ctr,4), Y_centroid=round(y_ctr,4), n_points=n_tree_pts,
    H_max_m=round(h_max,4), H_mean_m=round(h_mean,4), H_sd_m=round(h_sd,4),
    H_p25_m=round(h_p25,4), H_p50_m=round(h_p50,4),
    H_p75_m=round(h_p75,4), H_p95_m=round(h_p95,4),
    DBH_cm=round(dbh_cm,4), I_mean=round(i_mean_v,4), I_max=round(i_max_v,4),
    AGB_teak_kg=round(agb_teak,4),  AGC_teak_kg=round(agc_teak,4),
    AGB_Kom_kg=round(agb_kom,4),    BGB_Kom_kg=round(bgb_kom,4),
    AGC_Kom_kg=round(agc_kom,4),
    AGB_Chave_kg=round(agb_chave,4), AGC_Chave_kg=round(agc_chave,4),
    stringsAsFactors=FALSE)
  write.csv(attr_row,
            file.path(dir_trees_csv, sprintf("tree_%04d_attributes.csv",tid)),
            row.names=FALSE)
  all_tree_attrs[[i]] <- attr_row
}
options(lidR.progress = TRUE)

cat(sprintf("\n  Exported=%d | Skipped=%d | Errors=%d\n\n",
            n_exported, n_skipped, n_error))
if (length(all_tree_attrs) > 0) {
  master_attrs <- do.call(rbind, Filter(Negate(is.null), all_tree_attrs))
  write.csv(master_attrs,
            file.path(dir_trees,"ALL_Trees_MasterAttributes.csv"), row.names=FALSE)
  cat("  Master CSV saved.\n\n")
} else {
  cat("  WARNING: No trees exported. Master CSV not created.\n\n")
}

# =============================================================================
# STEP 14: Cloud & Pixel Metrics
# =============================================================================
cat("[Step 14] Computing Cloud & Pixel Metrics...\n")

cloud_met <- cloud_metrics(las_norm, func = ~list(
  z_mean=mean(Z), z_max=max(Z), z_sd=sd(Z),
  z_cv=ifelse(mean(Z) > 0, sd(Z)/mean(Z), NA_real_),
  z_p10=quantile(Z,0.10), z_p25=quantile(Z,0.25), z_p50=quantile(Z,0.50),
  z_p75=quantile(Z,0.75), z_p90=quantile(Z,0.90),
  z_p95=quantile(Z,0.95), z_p99=quantile(Z,0.99),
  pct_above2=sum(Z>2)/length(Z)*100,
  pct_above5=sum(Z>5)/length(Z)*100,
  pct_above10=sum(Z>10)/length(Z)*100,
  i_mean=mean(Intensity), i_sd=sd(Intensity), n_points=length(Z)))
write.csv(as.data.frame(cloud_met),
          file.path(output_dir,"Metrics_CloudLevel.csv"), row.names=FALSE)

f_metrics <- function(z, i, rn) list(
  z_mean=mean(z), z_max=max(z), z_sd=sd(z),
  z_cv=ifelse(mean(z,na.rm=TRUE)>0, sd(z,na.rm=TRUE)/mean(z,na.rm=TRUE), NA_real_),
  z_p25=quantile(z,0.25,na.rm=TRUE), z_p50=quantile(z,0.50,na.rm=TRUE),
  z_p75=quantile(z,0.75,na.rm=TRUE), z_p95=quantile(z,0.95,na.rm=TRUE),
  i_mean=mean(i,na.rm=TRUE), i_sd=sd(i,na.rm=TRUE),
  pct_1st_rtn=ifelse(length(rn)>0, sum(rn==1L,na.rm=TRUE)/length(rn)*100, 0),
  cover_2m=sum(z>2,na.rm=TRUE)/length(z)*100,
  cover_5m=sum(z>5,na.rm=TRUE)/length(z)*100,
  n_points=length(z))
pixel_met <- pixel_metrics(las=las_norm,
                            func=~f_metrics(Z,Intensity,ReturnNumber), res=10)
terra::writeRaster(pixel_met,
                   file.path(dir_rasters,"Metrics_PixelLevel_10m.tif"),
                   overwrite=TRUE)
cat(sprintf("  Saved: CloudMetrics + PixelMetrics (%d layers)\n\n",
            terra::nlyr(pixel_met)))

# =============================================================================
# STEP 15: Plots 01-16
# [v8.2-FIX1] create_seg_raster() now uses crown polygon sf, not LAS points
# [v8.2-FIX2] gc() between heavy plots
# =============================================================================
cat("[Step 15] Generating plots 01-16...\n")

create_pastel_palette <- function(n_colors) {
  n_colors <- max(1L, as.integer(n_colors))
  hcl(h = seq(0, 360, length.out = n_colors + 1)[1:n_colors], c = 70, l = 75)
}

# --- Plot 01: CHM Comparison ---
png(file.path(dir_plots,"01_CHM_Comparison.png"), width=1400, height=600, res=120)
par(mfrow=c(1,3), mar=c(3,3,3,1))
plot(chm_p2r,     main="CHM - P2R",      col=height.colors(50))
plot(chm_pitfree, main="CHM - Pitfree",  col=height.colors(50))
plot(chm_smooth,  main="CHM - Smoothed", col=height.colors(50))
dev.off(); gc(); cat("  Plot 01: CHM Comparison\n")

create_seg_raster <- function(tree_met_sf, chm_template, label = "") {
  if (is.null(tree_met_sf) || nrow(tree_met_sf) == 0) return(NULL)
  crowns_vect <- terra::vect(tree_met_sf[, "treeID"])
  r_out       <- terra::rasterize(crowns_vect, chm_template,
                                  field = "treeID", fun = "max")
  return(r_out)
}

# --- Plot 02: Tree Detection Comparison ---
png(file.path(dir_plots,"02_TreeDetection_Comparison.png"), width=1400, height=500, res=120)
par(mfrow=c(1,3), mar=c(3,3,3,1))
plot(chm_smooth, main=paste0("Fixed 2m (n=",nrow(ttops_fixed_2m),")"), col=height.colors(50))
plot(sf::st_geometry(ttops_fixed_2m), add=TRUE, pch=3, col="red",   cex=0.5)
plot(chm_smooth, main=paste0("Fixed 3m (n=",nrow(ttops_fixed_3m),")"), col=height.colors(50))
plot(sf::st_geometry(ttops_fixed_3m), add=TRUE, pch=3, col="blue",  cex=0.5)
plot(chm_smooth, main=paste0("Variable WS (n=",nrow(ttops_variable),")"), col=height.colors(50))
plot(sf::st_geometry(ttops_variable), add=TRUE, pch=3, col="green", cex=0.5)
dev.off(); gc(); cat("  Plot 02: Tree Detection\n")


# --- Plot 02B: Segmentation CHM — 2 Algorithms (polygon-based, memory-safe) ---
# [v8.2-FIX1] Uses create_seg_raster(tree_met_sf) — no LAS points loaded
png(file.path(dir_plots,"02B_Segmentation_CHM_Algorithms.png"), width=1100, height=500, res=120)
par(mfrow=c(1,2), mar=c(3,3,3,2))

seg_rast_dalp <- create_seg_raster(tree_met_dalp, chm_smooth, "Dalponte2016")
seg_rast_silv <- create_seg_raster(tree_met_silv, chm_smooth, "Silva2016")

if (!is.null(seg_rast_dalp)) {
  n_ids_dalp <- length(unique(na.omit(terra::values(seg_rast_dalp))))
  plot(seg_rast_dalp,
       main  = paste0("A) Dalponte2016 (n=", n_dalp, ")"),
       col   = create_pastel_palette(n_ids_dalp),
       legend = FALSE)
} else {
  plot(chm_smooth, main="A) Dalponte2016 - No segmentation")
}

if (!is.null(seg_rast_silv)) {
  n_ids_silv <- length(unique(na.omit(terra::values(seg_rast_silv))))
  plot(seg_rast_silv,
       main  = paste0("B) Silva2016 (n=", n_silv, ")"),
       col   = create_pastel_palette(n_ids_silv),
       legend = FALSE)
} else {
  plot(chm_smooth, main="B) Silva2016 - No segmentation")
}

dev.off()
rm(seg_rast_dalp, seg_rast_silv); gc()
cat("  Plot 02B: Segmentation CHM (2 Algorithms)\n")

# --- Plot 03: AGB & AGC Maps ---
png(file.path(dir_plots,"03_AGB_AGC_Maps.png"), width=1600, height=1200, res=120)
par(mfrow=c(2,3), mar=c(3,3,3,2))
plot(r_AGB_teak,  main="AGB Teak-Soraya2025 [Mg/ha]",  col=viridis::viridis(50))
plot(r_AGC_teak,  main="AGC Teak-Soraya2025 [MgC/ha]", col=viridis::magma(50))
plot(r_AGB_Kom,   main="AGB Komiyama2005 [Mg/ha]",     col=viridis::viridis(50))
plot(r_AGC_Kom,   main="AGC Komiyama2005 [MgC/ha]",    col=viridis::magma(50))
plot(r_AGB_Chave, main="AGB Chave2014 [Mg/ha]",        col=viridis::viridis(50))
plot(r_AGC_Chave, main="AGC Chave2014 [MgC/ha]",       col=viridis::magma(50))
dev.off(); gc(); cat("  Plot 03: AGB/AGC Maps\n")

# --- Plot 04: DBH Map ---
png(file.path(dir_plots,"04_DBH_Map.png"), width=800, height=700, res=120)
plot(r_DBH, main="Estimated DBH (cm)", col=viridis::plasma(50))
dev.off(); gc(); cat("  Plot 04: DBH Map\n")

# --- Plot 05: Height Distribution ---
p05 <- ggplot(tree_met_csv, aes(x=z_max)) +
  geom_histogram(bins=40, fill="#2E8B57", color="white", alpha=0.85) +
  geom_vline(xintercept=mean(tree_met_csv$z_max,na.rm=TRUE),
             color="red", linetype="dashed", linewidth=1) +
  annotate("text", x=mean(tree_met_csv$z_max,na.rm=TRUE)*1.05, y=Inf, vjust=2,
           label=paste0("Mean=",round(mean(tree_met_csv$z_max,na.rm=TRUE),1),"m"),
           color="red", size=4) +
  labs(title="Tree Height Distribution",
       subtitle=paste0("n=",nrow(tree_met_csv)," trees | Silva2016"),
       x="Height z_max (m)", y="Count") +
  theme_minimal(base_size=13)
ggsave(file.path(dir_plots,"05_TreeHeight_Distribution.png"),
       p05, width=8, height=5, dpi=150)
cat("  Plot 05: Height Distribution\n")

# --- Plot 06: DBH Distribution ---
p06 <- ggplot(tree_met_csv, aes(x=DBH_cm)) +
  geom_histogram(bins=40, fill="#8B4513", color="white", alpha=0.85) +
  geom_vline(xintercept=mean(tree_met_csv$DBH_cm,na.rm=TRUE),
             color="red", linetype="dashed", linewidth=1) +
  annotate("text", x=mean(tree_met_csv$DBH_cm,na.rm=TRUE)*1.05, y=Inf, vjust=2,
           label=paste0("Mean=",round(mean(tree_met_csv$DBH_cm,na.rm=TRUE),1),"cm"),
           color="red", size=4) +
  labs(title="Estimated DBH Distribution",
       subtitle="Feldpausch et al. (2011) H-DBH allometry",
       x="DBH (cm)", y="Count") +
  theme_minimal(base_size=13)
ggsave(file.path(dir_plots,"06_DBH_Distribution.png"),
       p06, width=8, height=5, dpi=150)
cat("  Plot 06: DBH Distribution\n")

# --- Plot 07: AGB Model Comparison ---
tree_long <- tree_met_csv %>%
  dplyr::select(treeID, z_max, AGB_teak_kg, AGB_Kom_kg, AGB_Chave_kg) %>%
  tidyr::pivot_longer(cols=c(AGB_teak_kg, AGB_Kom_kg, AGB_Chave_kg),
                      names_to="Model", values_to="AGB_kg") %>%
  mutate(Model = dplyr::recode(Model,
    "AGB_teak_kg"  = "Teak (Soraya 2025)",
    "AGB_Kom_kg"   = "Komiyama (2005)",
    "AGB_Chave_kg" = "Chave (2014)"))
p07 <- ggplot(tree_long, aes(x=z_max, y=AGB_kg/1000, color=Model)) +
  geom_point(alpha=0.35, size=1.0) +
  geom_smooth(method="loess", se=FALSE, linewidth=1.2) +
  scale_color_manual(values=c("#2E8B57","#4682B4","#CD853F")) +
  labs(title="AGB Model Comparison",
       subtitle="Soraya 2025 | Komiyama 2005 | Chave 2014",
       x="Tree Height (m)", y="AGB (tons/tree)") +
  theme_minimal(base_size=13)
ggsave(file.path(dir_plots,"07_AGB_ModelComparison.png"),
       p07, width=9, height=6, dpi=150)
cat("  Plot 07: AGB Model Comparison\n")

# --- Plot 08: AGC Model Comparison ---
agc_long <- tree_met_csv %>%
  dplyr::select(treeID, z_max, AGC_teak_kg, AGC_Kom_kg, AGC_Chave_kg) %>%
  tidyr::pivot_longer(cols=c(AGC_teak_kg, AGC_Kom_kg, AGC_Chave_kg),
                      names_to="Model", values_to="AGC_kg") %>%
  mutate(Model = dplyr::recode(Model,
    "AGC_teak_kg"  = "Teak (Soraya 2025)",
    "AGC_Kom_kg"   = "Komiyama (2005)",
    "AGC_Chave_kg" = "Chave (2014)"))
p08 <- ggplot(agc_long, aes(x=z_max, y=AGC_kg, color=Model)) +
  geom_point(alpha=0.35, size=1.0) +
  geom_smooth(method="loess", se=FALSE, linewidth=1.2) +
  scale_color_manual(values=c("#2E8B57","#4682B4","#CD853F")) +
  labs(title="AGC Model Comparison",
       subtitle="Carbon fraction CF=0.47 (IPCC 2006)",
       x="Tree Height (m)", y="AGC (kg C/tree)") +
  theme_minimal(base_size=13)
ggsave(file.path(dir_plots,"08_AGC_ModelComparison.png"),
       p08, width=9, height=6, dpi=150)
cat("  Plot 08: AGC Model Comparison\n")

# --- Plot 09: AGB + BGB Top 20 Tallest Trees ---
top20 <- tree_met_csv %>%
  dplyr::arrange(desc(z_max)) %>%
  dplyr::slice(1:min(20, nrow(tree_met_csv))) %>%
  dplyr::select(treeID, AGB_Kom_kg, BGB_Kom_kg) %>%
  tidyr::pivot_longer(cols=c(AGB_Kom_kg, BGB_Kom_kg),
                      names_to="Component", values_to="Biomass_kg") %>%
  mutate(Component = dplyr::recode(Component,
    "AGB_Kom_kg" = "AGB",
    "BGB_Kom_kg" = "BGB"))
p09 <- ggplot(top20, aes(x=factor(treeID), y=Biomass_kg, fill=Component)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("AGB"="#2E8B57","BGB"="#8B4513")) +
  labs(title="AGB & BGB — Top 20 Tallest Trees",
       subtitle="Komiyama et al. (2005)",
       x="Tree ID", y="Biomass (kg)") +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(dir_plots,"09_AGB_BGB_Top20.png"),
       p09, width=10, height=6, dpi=150)
cat("  Plot 09: AGB+BGB Top 20\n")

# --- Plot 10: Tree Crown Spatial Map (Silva2016 primary) ---
p10 <- ggplot() +
  geom_sf(data=tree_met, aes(fill=z_max),
          color="white", alpha=0.7, linewidth=0.1) +
  geom_sf(data=sf::st_zm(ttops_attr, drop=TRUE),
          aes(color=DBH_cm), size=0.8) +
  scale_fill_viridis_c(name="Height\n(m)",  option="viridis") +
  scale_color_viridis_c(name="DBH\n(cm)",   option="plasma") +
  labs(title="Tree Crown Map — Silva2016",
       subtitle=paste0("n=",nrow(tree_met)," trees | ",
                       round(total_area_ha,3)," ha"),
       x="Easting", y="Northing") +
  theme_minimal(base_size=11) +
  theme(axis.text=element_text(size=7))
ggsave(file.path(dir_plots,"10_TreeCrownMap.png"),
       p10, width=10, height=8, dpi=150)
cat("  Plot 10: Tree Crown Map\n")

# --- Plot 11: Crown Area Distribution ---
p11 <- ggplot(tree_met_csv, aes(x=Crown_area_m2)) +
  geom_histogram(bins=40, fill="#4682B4", color="white", alpha=0.85) +
  geom_vline(xintercept=median(tree_met_csv$Crown_area_m2, na.rm=TRUE),
             color="orange", linetype="dashed", linewidth=1) +
  labs(title="Crown Area Distribution",
       subtitle=paste0("Median=",
                       round(median(tree_met_csv$Crown_area_m2,na.rm=TRUE),1),
                       " m2"),
       x="Crown Area (m2)", y="Count") +
  theme_minimal(base_size=13)
ggsave(file.path(dir_plots,"11_CrownArea_Distribution.png"),
       p11, width=8, height=5, dpi=150)
cat("  Plot 11: Crown Area Distribution\n")

# --- Plot 12: Pixel Mean Height ---
png(file.path(dir_plots,"12_PixelMetrics_ZMean.png"), width=800, height=700, res=120)
plot(pixel_met[["z_mean"]], main="Mean Canopy Height per 10m Pixel (ABA)",
     col=viridis::viridis(50))
dev.off(); gc(); cat("  Plot 12: Pixel Mean Height\n")

# --- Plot 13: Precision / Recall / F1 / IoU ---
val_long <- val_master %>%
  dplyr::select(Algorithm, Precision, Recall, F1_Score, Mean_IoU) %>%
  tidyr::pivot_longer(cols=c(Precision, Recall, F1_Score, Mean_IoU),
                      names_to="Metric", values_to="Value")

p13 <- ggplot(val_long, aes(x=Metric, y=Value, fill=Algorithm)) +
  geom_bar(stat="identity", position="dodge", width=0.55) +
  geom_text(aes(label=round(Value,3)),
            position=position_dodge(width=0.55), vjust=-0.4, size=3.5) +
  scale_fill_manual(values=c("Dalponte2016"="#4682B4")) +
  scale_y_continuous(limits=c(0, 1.15)) +
  labs(title="Segmentation Validation: Precision / Recall / F1 / IoU",
       subtitle=paste0("Reference = Silva2016 | match_dist = ",
                       match_dist_m," m"),
       x="Metric", y="Score (0-1)") +
  theme_minimal(base_size=13)
ggsave(file.path(dir_plots,"13_Validation_PRF1_IoU.png"),
       p13, width=8, height=6, dpi=150)
cat("  Plot 13: Validation PRF1+IoU\n")

# --- Plot 14: IoU Distribution histogram ---
if (!is.null(iou_df) && nrow(iou_df) > 0) {
  iou_plot_df   <- iou_df %>% dplyr::filter(!is.na(IoU))
  if (nrow(iou_plot_df) > 0) {
    mean_iou_dalp <- mean(iou_plot_df$IoU, na.rm=TRUE)
    p14 <- ggplot(iou_plot_df, aes(x=IoU)) +
      geom_histogram(bins=30, fill="#4682B4", color="white", alpha=0.8) +
      geom_vline(xintercept=mean_iou_dalp,
                 color="#1a3a5c", linetype="dashed", linewidth=1) +
      annotate("text",
               x=mean_iou_dalp, y=Inf,
               label=paste0("μ=",round(mean_iou_dalp,3)),
               vjust=2, hjust=-0.15, size=4, color="#1a3a5c") +
      labs(title="IoU Distribution — Dalponte2016 vs Silva2016",
           subtitle="Crown polygon overlap | Reference = Silva2016",
           x="IoU (Intersection over Union)", y="Count") +
      theme_minimal(base_size=13)
    ggsave(file.path(dir_plots,"14_Validation_IoU_Distribution.png"),
           p14, width=8, height=5, dpi=150)
    cat("  Plot 14: IoU Distribution\n")
  }
} else {
  cat("  Plot 14: Skipped (no IoU data available)\n")
}

# --- Plot 15: RMSE & MAE ---
err_df <- val_master %>%
  dplyr::select(Algorithm, RMSE_H_m, MAE_H_m, RMSE_CA_m2, MAE_CA_m2) %>%
  tidyr::pivot_longer(cols=c(RMSE_H_m, MAE_H_m, RMSE_CA_m2, MAE_CA_m2),
                      names_to="Metric", values_to="Value")

p15 <- ggplot(err_df, aes(x=Metric, y=Value, fill=Algorithm)) +
  geom_bar(stat="identity", position="dodge", width=0.55) +
  geom_text(aes(label=round(Value,2)),
            position=position_dodge(width=0.55), vjust=-0.4, size=3.5) +
  scale_fill_manual(values=c("Dalponte2016"="#4682B4")) +
  labs(title="RMSE & MAE — Tree Height and Crown Area",
       subtitle="Dalponte2016 matched trees vs Silva2016 reference",
       x="Metric", y="Value  (m  or  m2)") +
  theme_minimal(base_size=13)
ggsave(file.path(dir_plots,"15_Validation_RMSE_MAE.png"),
       p15, width=8, height=6, dpi=150)
cat("  Plot 15: RMSE & MAE\n")

# --- Plot 16: 2-Panel Crown Map Comparison ---
cat("  Plot 16: 2-Panel Crown Map Comparison...\n")

bbox_ext <- sf::st_bbox(tree_met_silv)

make_crown_map <- function(met, title_str,
                            fill_var   = "z_max",
                            fill_label = "Max Height (m)",
                            palette    = "YlGn") {
  ggplot(met) +
    geom_sf(aes(fill = .data[[fill_var]]),
            color="white", linewidth=0.12, alpha=0.88) +
    scale_fill_distiller(
      palette   = palette,
      direction = 1,
      name      = fill_label,
      na.value  = "grey80",
      guide     = guide_colorbar(
        barheight = unit(3.5,"cm"),
        barwidth  = unit(0.4,"cm"))) +
    coord_sf(
      xlim   = c(bbox_ext["xmin"], bbox_ext["xmax"]),
      ylim   = c(bbox_ext["ymin"], bbox_ext["ymax"]),
      expand = FALSE) +
    labs(title    = title_str,
         subtitle = paste0("n = ", nrow(met), " trees")) +
    theme_minimal(base_size=10) +
    theme(
      plot.title      = element_text(face="bold", hjust=0.5, size=11),
      plot.subtitle   = element_text(hjust=0.5, size=9, color="grey40"),
      legend.position = "right",
      axis.text       = element_text(size=6.5),
      axis.title      = element_blank(),
      panel.grid      = element_line(color="grey92", linewidth=0.25),
      panel.border    = element_rect(color="grey50", fill=NA, linewidth=0.5))
}

p16a <- make_crown_map(tree_met_dalp, "A) Dalponte2016",
                        fill_var="z_max", fill_label="Max Height (m)",
                        palette="YlGn")
p16b <- make_crown_map(tree_met_silv, "B) Silva2016  [Primary]",
                        fill_var="z_max", fill_label="Max Height (m)",
                        palette="YlOrRd")

p16_combined <- (p16a | p16b) +
  plot_annotation(
    title    = "Figure 16 — Individual Tree Crown Maps: Dalponte2016 vs Silva2016",
    subtitle = paste0(
      "Site: ", site_name, "  |  ",
      "Dalponte2016: ", nrow(tree_met_dalp), " trees  |  ",
      "Silva2016: ",    nrow(tree_met_silv), " trees"),
    caption  = paste0(
      "Crown polygons (convex hull) · lidR ITS algorithms · ",
      "Primary = Silva2016 · EPSG:", target_epsg),
    theme = theme(
      plot.title    = element_text(face="bold", size=13, hjust=0.5),
      plot.subtitle = element_text(size=9,  hjust=0.5, color="grey30"),
      plot.caption  = element_text(size=7.5, hjust=1,  color="grey50")))

ggsave(
  filename = file.path(dir_plots,"16_CrownMap_AlgoComparison_2Panel.png"),
  plot     = p16_combined,
  width=14, height=7, dpi=300, bg="white")
rm(p16a, p16b, p16_combined); gc()
cat("  Plot 16: 2-Panel Crown Map saved.\n\n")

# =============================================================================
# STEP 17: Final Summary Report (console + CSV)
# =============================================================================
cat("[Step 17] Writing Final Summary Report...\n")
cat(strrep("=", 60), "\n")

total_AGB_teak_Mg  <- sum(tree_met_csv$AGB_teak_kg,  na.rm=TRUE) / 1000
total_AGC_teak_Mg  <- sum(tree_met_csv$AGC_teak_kg,  na.rm=TRUE) / 1000
total_AGB_Kom_Mg   <- sum(tree_met_csv$AGB_Kom_kg,   na.rm=TRUE) / 1000
total_BGB_Kom_Mg   <- sum(tree_met_csv$BGB_Kom_kg,   na.rm=TRUE) / 1000
total_AGC_Kom_Mg   <- sum(tree_met_csv$AGC_Kom_kg,   na.rm=TRUE) / 1000
total_AGB_Chave_Mg <- sum(tree_met_csv$AGB_Chave_kg, na.rm=TRUE) / 1000
total_AGC_Chave_Mg <- sum(tree_met_csv$AGC_Chave_kg, na.rm=TRUE) / 1000

AGB_teak_Mgha  <- total_AGB_teak_Mg  / total_area_ha
AGC_teak_Mgha  <- total_AGC_teak_Mg  / total_area_ha
AGB_Kom_Mgha   <- total_AGB_Kom_Mg   / total_area_ha
AGC_Kom_Mgha   <- total_AGC_Kom_Mg   / total_area_ha
AGB_Chave_Mgha <- total_AGB_Chave_Mg / total_area_ha
AGC_Chave_Mgha <- total_AGC_Chave_Mg / total_area_ha

mean_dbh <- mean(tree_met_csv$DBH_cm, na.rm=TRUE)
if (mean_dbh < 10)
  cat("  *** WARNING: Mean DBH < 10 cm. Verify allometric model",
      "against local field inventory. ***\n")

cat(sprintf("\n  Site              : %s\n",  site_name))
cat(sprintf("  Area (ha)         : %.4f\n", total_area_ha))
cat(sprintf("  Trees — Silva2016 : %d  [PRIMARY]\n", nrow(tree_met_csv)))
cat(sprintf("  Trees — Dalponte  : %d\n",   nrow(tree_met_dalp)))
cat(sprintf("  Mean H_max (m)    : %.2f\n", mean(tree_met_csv$z_max,  na.rm=TRUE)))
cat(sprintf("  Mean DBH (cm)     : %.2f\n", mean_dbh))
cat(sprintf("  Mean Crown (m2)   : %.2f\n", mean(tree_met_csv$Crown_area_m2, na.rm=TRUE)))
cat("\n  --- Biomass & Carbon (Silva2016 primary) ---\n")
cat(sprintf("  [Soraya 2025]  AGB total: %8.2f Mg  |  %7.2f Mg/ha\n",
            total_AGB_teak_Mg,  AGB_teak_Mgha))
cat(sprintf("                 AGC total: %8.2f Mg  |  %7.2f MgC/ha\n",
            total_AGC_teak_Mg,  AGC_teak_Mgha))
cat(sprintf("  [Komiyama2005] AGB total: %8.2f Mg  |  %7.2f Mg/ha\n",
            total_AGB_Kom_Mg,   AGB_Kom_Mgha))
cat(sprintf("                 BGB total: %8.2f Mg\n", total_BGB_Kom_Mg))
cat(sprintf("                 AGC total: %8.2f Mg  |  %7.2f MgC/ha\n",
            total_AGC_Kom_Mg,   AGC_Kom_Mgha))
cat(sprintf("  [Chave 2014]   AGB total: %8.2f Mg  |  %7.2f Mg/ha\n",
            total_AGB_Chave_Mg, AGB_Chave_Mgha))
cat(sprintf("                 AGC total: %8.2f Mg  |  %7.2f MgC/ha\n",
            total_AGC_Chave_Mg, AGC_Chave_Mgha))
cat("\n  --- Cross-Algorithm Validation (Dalponte2016 vs Silva2016) ---\n")
print(val_master[, c("Algorithm","N_ref","N_detected",
                      "Precision","Recall","F1_Score","Mean_IoU",
                      "RMSE_H_m","MAE_H_m")])

# --- Save summary CSV ---
summary_df <- data.frame(
  Site                  = site_name,
  Area_ha               = round(total_area_ha, 4),
  N_trees_Silva2016     = nrow(tree_met_csv),
  N_trees_Dalponte2016  = nrow(tree_met_dalp),
  Mean_H_max_m          = round(mean(tree_met_csv$z_max,         na.rm=TRUE), 3),
  Mean_DBH_cm           = round(mean(tree_met_csv$DBH_cm,        na.rm=TRUE), 3),
  Mean_Crown_m2         = round(mean(tree_met_csv$Crown_area_m2, na.rm=TRUE), 3),
  AGB_Soraya_total_Mg   = round(total_AGB_teak_Mg,  3),
  AGB_Soraya_Mgha       = round(AGB_teak_Mgha,      3),
  AGC_Soraya_total_Mg   = round(total_AGC_teak_Mg,  3),
  AGC_Soraya_Mgha       = round(AGC_teak_Mgha,      3),
  AGB_Kom_total_Mg      = round(total_AGB_Kom_Mg,   3),
  AGB_Kom_Mgha          = round(AGB_Kom_Mgha,       3),
  BGB_Kom_total_Mg      = round(total_BGB_Kom_Mg,   3),
  AGC_Kom_total_Mg      = round(total_AGC_Kom_Mg,   3),
  AGC_Kom_Mgha          = round(AGC_Kom_Mgha,       3),
  AGB_Chave_total_Mg    = round(total_AGB_Chave_Mg, 3),
  AGB_Chave_Mgha        = round(AGB_Chave_Mgha,     3),
  AGC_Chave_total_Mg    = round(total_AGC_Chave_Mg, 3),
  AGC_Chave_Mgha        = round(AGC_Chave_Mgha,     3),
  Val_Dalp_Precision    = val_dalp$metrics$Precision,
  Val_Dalp_Recall       = val_dalp$metrics$Recall,
  Val_Dalp_F1           = val_dalp$metrics$F1_Score,
  Val_Dalp_IoU          = val_dalp$metrics$Mean_IoU,
  Val_Dalp_RMSE_H       = val_dalp$metrics$RMSE_H_m,
  Val_Dalp_MAE_H        = val_dalp$metrics$MAE_H_m,
  Val_Dalp_RMSE_CA      = val_dalp$metrics$RMSE_CA_m2,
  Val_Dalp_MAE_CA       = val_dalp$metrics$MAE_CA_m2,
  Wood_density_gcc      = wood_density,
  CF_IPCC               = CF,
  EPSG                  = target_epsg,
  Processed_at          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors      = FALSE
)

write.csv(summary_df,
          file.path(output_dir,"SUMMARY_FinalReport.csv"),
          row.names=FALSE)
cat("\n  Saved: SUMMARY_FinalReport.csv\n")

cat(strrep("=", 60), "\n")
cat("=== LiDAR Pipeline v8.2 COMPLETE ===\n")
cat(sprintf("    Site     : %s\n",  site_name))
cat(sprintf("    Finished : %s\n",  format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("    Output   : %s\n",  output_dir))
cat(strrep("=", 60), "\n")