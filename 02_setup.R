# ST920 | Childhood malaria in The Gambia
# 02_setup.R

# Analysis settings: change these before running the script.
refit_stan <- FALSE       # Reuse saved fit only for identical data/model/priors
parallel_chains <- 2L
M <- 1000L                # Number of joint predictive realisations
threshold <- 0.30         # Illustrative prevalence threshold
rho0 <- 30                # P(rho < rho0) = prior_tail
sigma0 <- 2               # P(sigma > sigma0) = prior_tail
prior_tail <- 0.05
beta_sd <- 5             
jitter <- 1e-8 
grid_spacing <- 5    # Prediction grid spacing (km)
map_buffer <- 5      # Margin around the full map (km)

dir.create("images", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
# Aggregate individual binary outcomes into village-level binomial counts.
data(gambia, package = "geoR")
stopifnot(all(c("x", "y", "pos") %in% names(gambia)),
          all(na.omit(gambia$pos) %in% c(0, 1)))

dat <- gambia |>
  filter(!is.na(pos), !is.na(x), !is.na(y)) |>
  group_by(x, y) |>
  summarise(Y = sum(pos), N = n(), prevalence = Y / N,
            .groups = "drop") |>
  mutate(x = x / 1000, y = y / 1000)

if (!exists("gambia.borders")) data(gambia.borders, package = "geoR")
border <- as.data.frame(gambia.borders)[, 1:2] / 1000
names(border) <- c("x", "y")
# Missing rows separate geographical features, including river lines.
missing <- !is.finite(border$x) | !is.finite(border$y)
border$segment <- cumsum(missing)
border <- border[!missing, , drop = FALSE]
stopifnot(nrow(border) > 0, all(is.finite(border$x)), all(is.finite(border$y)))

coords <- as.matrix(dat[, c("x", "y")])
n <- nrow(dat)
D <- as.matrix(dist(coords))
print(c(villages = n, children = sum(dat$N),
        prevalence = sum(dat$Y) / sum(dat$N)))

# Common map layers: reuse these throughout the analysis.
map_theme <- theme_void(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")
boundary <- geom_path(
  data = border, 
  aes(x, y, group = segment), 
  inherit.aes = FALSE,
  colour = "grey35", 
  linewidth = 0.35)
boundary_white <- geom_path(
  data = border,
  aes(x, y, group = segment),
  inherit.aes = FALSE,
  colour = "white",
  linewidth = 0.5
)
base_map <- ggplot() + boundary + coord_equal() + map_theme

# Compare positive counts, sample sizes and observed proportions.
p_counts <- base_map +
  geom_point(data = dat, aes(x, y, size = Y), colour = "darkorange3") +
  labs(title = "Number positive", size = "Positive")
p_sample <- base_map +
  geom_point(data = dat, aes(x, y, size = N), colour = "steelblue") +
  labs(title = "Number tested", size = "Tested")
p_prev <- base_map +
  geom_point(data = dat, aes(x, y, colour = prevalence), size = 2.5) +
  scale_colour_viridis_c(limits = c(0, 1), labels = label_percent()) +
  labs(title = "Observed prevalence", colour = "Prevalence")
ggsave("images/gambia_data.pdf", p_counts / p_sample / p_prev,
       width = 10, height = 5.5)
# PC priors in two dimensions: sigma ~ exponential(lambda_sigma),
# rho ~ inverse-gamma(1, lambda_rho). The intercept has a normal prior.
lambda_sigma <- -log(prior_tail) / sigma0
lambda_rho <- -log(prior_tail) * rho0

# Rectangular prediction domain covering all geographical features and villages.
# Coordinates are planar kilometres; no external CRS assumption is needed.
# Round outwards to grid multiples so the full map and margin are covered.
stopifnot(grid_spacing > 0, map_buffer >= 0)
xlim <- range(c(border$x, dat$x)) + c(-map_buffer, map_buffer)
ylim <- range(c(border$y, dat$y)) + c(-map_buffer, map_buffer)
xlim <- c(floor(xlim[1] / grid_spacing), ceiling(xlim[2] / grid_spacing)) * grid_spacing
ylim <- c(floor(ylim[1] / grid_spacing), ceiling(ylim[2] / grid_spacing)) * grid_spacing

study_domain <- st_as_sfc(st_bbox(c(
  xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]
)))

# Keep every grid point: no national-boundary or concave-hull mask.
grid <- expand.grid(x = seq(xlim[1], xlim[2], by = grid_spacing),
                    y = seq(ylim[1], ylim[2], by = grid_spacing))
pred_coords <- as.matrix(grid[, c("x", "y")])
n_pred <- nrow(grid)
stopifnot(n_pred > 0, all(is.finite(pred_coords)))
cat("Number of prediction locations:", n_pred, "\n")

# Inspect the prediction domain before running predictions.
p_domain <- ggplot() +
  geom_point(data = grid, aes(x, y), colour = "steelblue", size = 0.8) +
  boundary +
  geom_point(data = dat, aes(x, y), colour = "red", size = 1.5) +
  coord_equal() + theme_bw() +
  labs(title = "Rectangular study domain",
       subtitle = "Red: observed villages; blue: prediction locations",
       x = "Easting (km)", y = "Northing (km)")
print(p_domain)
ggsave("images/gambia_prediction_domain.pdf", p_domain, width = 10, height = 4.5)


# Common raster-map functions for both dense GP and INLA results.
raster_map <- function(data, variable, title) {
  ggplot(data, aes(x, y)) + geom_raster(aes(fill = .data[[variable]])) +
    boundary_white + coord_equal() + labs(title = title, subtitle = "Rectangular study domain (not restricted to national land boundary)") + map_theme
}
make_map <- function(data, variable, title, legend, limits = NULL) {
  raster_map(data, variable, title) +
    scale_fill_viridis_c(limits = limits, labels = label_percent(), name = legend)
}
signed_map <- function(data, variable, title) {
  lim <- max(abs(data[[variable]]), 1e-10)
  raster_map(data, variable, title) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, limits = c(-lim, lim),
                         labels = label_percent(), name = "Difference")
}

# Shared settings recorded with prediction results for comparison.
analysis_spec <- list(data = dat, grid = grid, rho0 = rho0, sigma0 = sigma0,
                      prior_tail = prior_tail, beta_sd = beta_sd,
                      threshold = threshold, jitter = jitter)


