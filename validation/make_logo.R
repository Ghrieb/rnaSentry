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
# Coordinates follow the rnaSentry-logo-preview.svg convention (900 x 1000,
# y-down: y = 0 is the TOP, y = 1000 is the BOTTOM). The y axis is reversed
# explicitly with scale_y_reverse() so R's plot matches the SVG.
#
# Round 25.7: fixed vertical orientation (the previous coord_fixed(ylim =
# c(1000, 0)) trick silently rendered the axis NOT reversed, mirroring the
# whole composition top/bottom) and enlarged the emblem ~2x / text ~1.3x so
# the shield + KM curve fill the hexagon instead of floating tiny inside it.
#
# Output: man/figures/logo.png (1350 x 1500 px, transparent corners)

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
# ~2x the original 120x212 shield: width 260 (x 320..580), height 420
# (y 190..610), shoulder at y 368, Beziér controls pulled down to y 511 so
# the point keeps its guardrail taper.
shield_df <- rbind(
  data.frame(x = c(320, 580), y = c(190, 190)),          # top edge
  data.frame(x = 580, y = 368),                           # right shoulder
  bezier(580, 368, 580, 511, 450, 610),                   # right curve to point
  bezier(450, 610, 320, 511, 320, 368)                    # left curve back up
)

# ---- Kaplan-Meier step curve inside the shield (~2x, same 4-step shape) --
km_df <- data.frame(
  x = c(382, 427, 427, 464, 464, 505, 505, 521),
  y = c(285, 285, 355, 355, 438, 438, 515, 515)
)

# ---- build ----------------------------------------------------------------
p <- ggplot() +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits = c(0, 900), expand = c(0, 0)) +
  scale_y_reverse(limits = c(0, 1000), expand = c(0, 0)) +  # y-down, as the SVG
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
  annotate("text", x = 450, y = 695, label = "rnaSentry",
           family = "sans", fontface = "bold", size = 18.5, colour = paper) +
  annotate("text", x = 450, y = 745, label = "GUARDED  \u00b7  AUDITABLE",
           family = "mono", size = 5.7, colour = tagline_col)

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

# ---- verification (programmatic — round 25.7 fix) -----------------------
# scale_y_reverse => data y renders at row = y * 1.5 (y=0 top, y=1000 bottom),
# x renders at col = x * 1.5.
stopifnot(requireNamespace("png", quietly = TRUE))
im  <- png::readPNG(out)
nr  <- nrow(im); nc <- ncol(im)
R <- im[,,1]; G <- im[,,2]; B <- im[,,3]; A <- im[,,4]
at <- function(x, y) { r <- round(y * 1.5) + 1; c <- round(x * 1.5) + 1; im[r, c, ] }
near <- function(k, ref, tol = 0.06) all(abs(k[1:3] - col2rgb(ref)/255) < tol)

# 1. canvas + transparent corners
stopifnot(nr == 1500, nc == 1350)
stopifnot(A[2, 2] < 0.5, A[2, nc-1] < 0.5, A[nr-1, 2] < 0.5, A[nr-1, nc-1] < 0.5)

# 2. hexagon orientation: top vertex (450,40) outline near the TOP (row ~61),
#    bottom vertex (450,960) outline near the BOTTOM (row ~1441).
stopifnot(near(at(450, 40), line_white), near(at(450, 960), line_white))

# 3. shield: top edge (450,190) white outline; point (450,610) outline below it.
stopifnot(near(at(450, 190), line_white), near(at(450, 610), line_white))

# 4. KM curve teal, descending left->right (y increases with x in y-down).
stopifnot(near(at(382, 285), teal, tol = 0.10), near(at(521, 515), teal, tol = 0.10))

# 5. wordmark paper text sits BELOW the shield point and inside the hexagon.
paper_ok <- abs(R - 245/255) < 0.03 & abs(G - 246/255) < 0.03 &
            abs(B - 244/255) < 0.03 & A > 0.5
if (sum(paper_ok) > 0) {
  rr <- which(paper_ok, arr.ind = TRUE)
  stopifnot(min(rr[, 1]) > round(610 * 1.5) + 1,   # wordmark starts below shield point
            max(rr[, 2]) < 848 * 1.5, min(rr[, 2]) > 52 * 1.5)  # inside hexagon
  message("wordmark bbox rows ", min(rr[,1]), "-", max(rr[,1]),
          " cols ", min(rr[,2]), "-", max(rr[,2]))
} else {
  stop("no paper (wordmark) pixels found")
}

# 6. KM curve stays inside the hexagon interior (stroke extends ~3 px).
tl_ok <- abs(R - 79/255) < 0.06 & abs(G - 195/255) < 0.06 &
         abs(B - 199/255) < 0.06 & A > 0.5
rr <- which(tl_ok, arr.ind = TRUE)
stopifnot(nrow(rr) > 100,
          max(rr[, 2]) < (848 * 1.5) + 8, min(rr[, 2]) > (52 * 1.5) - 8,
          max(rr[, 1]) < (960 * 1.5) + 8, min(rr[, 1]) > (40 * 1.5) - 8)
message("All logo assertions passed (round 25.7 orientation + scale).")
