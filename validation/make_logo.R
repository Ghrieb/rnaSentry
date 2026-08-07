# validation/make_logo.R
#
# rnaSentry hex logo — redesign.
#
# Design: one idea, not a collage. A shield (the package's "guardrail"
# philosophy) containing a Kaplan-Meier step curve (the specific thing
# km_curve() produces — this is what makes the icon ownable rather than
# generic bio clip-art). Flat fills only, one accent color, no gradients
# or glow effects that degrade at small (favicon/README-badge) sizes.
#
# All coordinates below are copied directly from rnaSentry-logo-preview.svg
# (same 900 x 1000, y-down convention) — coord_fixed()'s ylim is given
# reversed (1000, 0) so R's plot matches the SVG without manual flipping.
#
# Output: man/figures/logo.png (transparent corners)

library(ggplot2)
library(ragg)

# ---- canvas -------------------------------------------------------------
res      <- 300
scale_px <- 1.5
canvas_w <- 900 * scale_px
canvas_h <- 1000 * scale_px

# ---- palette (shared with the rnaSentry web-app mockup / README) -------
navy_bottom <- "#0D1319"
paper       <- "#F5F6F4"
line_white  <- "#EDEFEC"
teal        <- "#4FC3C7"
tagline_col <- "#8FA5A8"

# ---- hexagon: perfect pointy-top hexagon, same vertices as the SVG -----
hexagon <- data.frame(
  x = c(450, 848, 848, 450, 52, 52),
  y = c(40, 270, 730, 960, 730, 270)
)

# ---- shield (the guardrail) — quadratic-Bezier sides, matches the SVG --
bezier <- function(x0, y0, cx, cy, x1, y1, n = 24) {
  t <- seq(0, 1, length.out = n)
  data.frame(
    x = (1 - t)^2 * x0 + 2 * (1 - t) * t * cx + t^2 * x1,
    y = (1 - t)^2 * y0 + 2 * (1 - t) * t * cy + t^2 * y1
  )
}
shield_df <- rbind(
  data.frame(x = c(390, 510), y = c(290, 290)),          # top edge
  data.frame(x = 510, y = 380),                           # right shoulder
  bezier(510, 380, 510, 452, 450, 502),                   # right curve to point
  bezier(450, 502, 390, 452, 390, 380)                    # left curve back up
)

# ---- Kaplan-Meier step curve inside the shield --------------------------
km_df <- data.frame(
  x = c(403, 434, 434, 460, 460, 488, 488, 499),
  y = c(328, 328, 360, 360, 398, 398, 433, 433)
)

# ---- build ----------------------------------------------------------------
p <- ggplot() +
  coord_fixed(xlim = c(0, 900), ylim = c(1000, 0), expand = FALSE) +
  theme_void() +
  theme(plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.margin      = margin(0, 0, 0, 0)) +
  geom_polygon(data = hexagon, aes(x, y), fill = navy_bottom,
               colour = line_white, linewidth = 4.2 * scale_px,
               linejoin = "round") +
  geom_polygon(data = shield_df, aes(x, y), fill = "#FFFFFF0D",
               colour = line_white, linewidth = 3 * scale_px,
               linejoin = "round") +
  geom_path(data = km_df, aes(x, y), colour = teal,
            linewidth = 3.3 * scale_px, lineend = "round",
            linejoin = "round") +
  annotate("text", x = 450, y = 705, label = "rnaSentry",
           family = "sans", fontface = "bold", size = 15.5, colour = paper) +
  annotate("text", x = 450, y = 752, label = "GUARDED  \u00b7  AUDITABLE",
           family = "mono", size = 4.4, colour = tagline_col)

# ---- optional: subtle native gradient background (ggplot2 >= 3.5) -------
# Replace the hex geom_polygon's `fill = navy_bottom` above with:
#   fill = grid::linearGradient(colours = c("#182430", navy_bottom),
#                                x1 = 0.5, y1 = 0, x2 = 0.5, y2 = 1)
# for a soft top-to-bottom depth instead of the flat fill (no manual
# raster/pixel math needed — this is native grid gradient support).

out <- file.path("man", "figures", "logo.png")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
agg_png(out, width = canvas_w, height = canvas_h, units = "px", res = res,
        background = "transparent")
print(p)
dev.off()

message("Wrote ", out, " (", canvas_w, " x ", canvas_h, " px)")
