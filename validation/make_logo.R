# validation/make_logo.R
# Generates the rnaSentry hex-sticker logo: perfect regular hexagon, radial
# electric-blue glow fading to deep navy, white hex frame line, minimalist
# single-strand RNA curve through 4 glowing data nodes, professional layout.
# Pure ggplot2 + ragg + grid (fully offline; no hexSticker needed).
# Output: man/figures/logo.png (1732 x 2000 px, transparent corners).

library(ggplot2)
library(ragg)
library(grid)

pdf(NULL)   # metrics device; prevents a stray Rplots.pdf in the repo root

# --- canvas / geometry -----------------------------------------------------
px_per_unit <- 1000
res         <- 300
canvas_w    <- round(sqrt(3) * 1000)   # 1732
canvas_h    <- 2000

hex_pts <- function(r) {
  th <- pi / 180 * (90 + 60 * (0:5))   # exact 6 vertices, no arcs
  data.frame(x = r * cos(th), y = r * sin(th))
}

# --- radial gradient (electric blue centre -> deep navy edge, hex-masked) ----
make_hex_gradient <- function(n = 800, from = "#8CCDF5", to = "#0B2A4A") {
  m  <- round(n * 0.8660)
  ys <- seq(1, -1, length.out = n)
  xs <- seq(-0.8660, 0.8660, length.out = m)
  g  <- expand.grid(x = xs, y = ys)
  d  <- sqrt(g$x^2 + g$y^2)            # radial distance, 1 at top/bottom vertex
  t  <- pmin(1, d)
  t  <- t * t * (3 - 2 * t)            # smoothstep for smooth falloff
  c1 <- col2rgb(from) / 255
  c2 <- col2rgb(to)   / 255
  col_r <- (1 - t) * c1[1] + t * c2[1]
  col_g <- (1 - t) * c1[2] + t * c2[2]
  col_b <- (1 - t) * c1[3] + t * c2[3]
  absy <- abs(g$y)
  hw   <- 0.8660 * ifelse(absy <= 0.5, 1, (1 - absy) / 0.5)
  inside <- abs(g$x) <= hw
  col   <- rgb(col_r, col_g, col_b, alpha = ifelse(inside, 1, 0))
  as.raster(matrix(col, nrow = n, ncol = m, byrow = TRUE))
}
grad_raster <- make_hex_gradient()

# --- text metrics -----------------------------------------------------------
fit_size_pt <- function(txt, target_units, bold = FALSE) {
  g <- grid::textGrob(txt, gp = grid::gpar(fontsize = 100, fontfamily = "sans",
                                           fontface = if (bold) "bold" else "plain"))
  w100 <- grid::convertWidth(grid::grobWidth(g), "pt", valueOnly = TRUE)
  target_units * px_per_unit * 100 / (w100 * res / 72)
}
size_mm <- function(S_pt) S_pt / 2.845276

w_pt <- fit_size_pt("rnaSentry",          0.55, bold = TRUE)
s_pt <- fit_size_pt("guarded & auditable", 0.55)
cat(sprintf("wordmark %.1f pt, subtitle %.1f pt\n", w_pt, s_pt))

# --- palette ---------------------------------------------------------------
white <- "#FFFFFF"
ice   <- "#B8DCE6"
cyan  <- "#7CCBE0"
navy  <- "#0B2A4A"

# --- shield (scaled down, held high) ----------------------------------------
shield <- data.frame(
  x = c(-0.26, -0.26, -0.15,  0.00,  0.15,  0.26,  0.26),
  y = c( 0.72,  0.50,  0.34,  0.18,  0.34,  0.50,  0.72)
)

# --- single-strand RNA curve through 4 glowing nodes -------------------------
xs_n <- c(-0.15, -0.05, 0.05, 0.15)
y_n  <- 0.45 + 0.095 * cos(pi * xs_n / 0.30)     # gentle hump through nodes
nodes <- data.frame(x = xs_n, y = y_n)

xs_c <- seq(-0.17, 0.17, length.out = 121)       # smooth curve (same formula)
curve_df <- data.frame(x = xs_c,
                       y = 0.45 + 0.095 * cos(pi * xs_c / 0.30))

node_r <- function(r_units) r_units * px_per_unit / (res / 25.4)   # mm diameter

# --- build -------------------------------------------------------------------
p <- ggplot() +
  coord_fixed(xlim = c(-0.866, 0.866), ylim = c(-1, 1), expand = FALSE) +
  theme_void() +
  theme(plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.margin      = margin(0, 0, 0, 0)) +
  # radial gradient face, masked to the perfect hexagon
  annotation_custom(rasterGrob(grad_raster, interpolate = TRUE),
                    xmin = -0.866, xmax = 0.866, ymin = -1, ymax = 1) +
  # crisp white frame line
  geom_polygon(data = hex_pts(1.0), aes(x, y), fill = NA,
               colour = white, linewidth = 3.2) +
  # shield (faint fill + clean white outline)
  geom_polygon(data = shield, aes(x, y), fill = "#FFFFFF10",
               colour = white, linewidth = 2.2) +
  # soft glow under the RNA curve
  geom_line(data = curve_df, aes(x, y), colour = cyan, alpha = 0.20,
            linewidth = 6.5, lineend = "round") +
  geom_line(data = curve_df, aes(x, y), colour = "#EAF6FB",
            linewidth = 2.3, lineend = "round") +
  # 4 glowing data nodes (outer glow + halo + white core)
  geom_point(data = nodes, aes(x, y), colour = cyan,  alpha = 0.30,
             size = node_r(0.034)) +
  geom_point(data = nodes, aes(x, y), colour = cyan,  alpha = 0.55,
             size = node_r(0.021)) +
  geom_point(data = nodes, aes(x, y), colour = white, fill = white,
             size = node_r(0.011), shape = 21, stroke = 0) +
  # typography: title at the optical centre, subtitle tucked just below
  annotate("text", x = 0, y = -0.05, label = "rnaSentry", family = "sans",
           size = size_mm(w_pt), fontface = "bold", colour = white) +
  annotate("text", x = 0, y = -0.30, label = "guarded & auditable",
           family = "sans", size = size_mm(s_pt), fontface = "plain",
           colour = ice)

out <- file.path("man", "figures", "logo.png")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
agg_png(out, width = canvas_w, height = canvas_h, units = "px", res = res,
        background = "transparent")
print(p)
dev.off()

# --- verification (programmatic) ---------------------------------------------
mw <- function(txt, S, bold = FALSE) grid::convertWidth(grid::grobWidth(
  grid::textGrob(txt, gp = grid::gpar(fontsize = S, fontfamily = "sans",
                                      fontface = if (bold) "bold" else "plain"))),
  "pt", valueOnly = TRUE)
pt2u <- function(pt) pt * (res / 72) / px_per_unit

sub_units <- pt2u(mw("guarded & auditable", s_pt))
half      <- sub_units / 2
edge_y    <- -1 + 0.5773503 * half
line_s    <- pt2u(72.07 * s_pt / 100)          # subtitle line height
sub_bot   <- -0.30 - line_s / 2
gap_edge  <- sub_bot - edge_y
line_w    <- pt2u(72.07 * w_pt / 100)
title_bot <- -0.05 + line_w / 2
gap_t_s   <- title_bot - (-0.30 - line_s / 2)   # title bottom - subtitle top
cat(sprintf("subtitle half-width %.3f u, clear gap to hex edge %.3f u (>0.10 ok)\n",
            half, gap_edge))
cat(sprintf("title-to-subtitle gap %.3f u (>0.08 ok)\n", gap_t_s))
stopifnot(gap_edge > 0.10, gap_t_s > 0.08, half < 0.866)
# nodes must stay inside the shield
stopifnot(all(abs(nodes$x) <= 0.26),
          all(nodes$y <= 0.72), all(nodes$y >= 0.18))
message("Wrote ", out, " (", canvas_w, " x ", canvas_h, " px)")
