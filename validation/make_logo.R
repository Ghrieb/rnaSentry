# validation/make_logo.R
# Generates the rnaSentry hex-sticker logo with base R graphics only (no
# external packages), so the asset is fully reproducible offline. Output:
# man/figures/logo.png (used by pkgdown and shown on the GitHub repo page).

rounded_hex <- function(cx, cy, r, round_r = NULL, n_arc = 24) {
  # Vertices of a pointy-top hexagon (corner at 90 deg), ordered counterclockwise.
  ang <- pi / 180 * (90 + 60 * (0:5))
  corners <- cbind(cx + r * cos(ang), cy + r * sin(ang))
  if (is.null(round_r)) round_r <- 0.16 * r
  rr <- min(round_r, 0.5 * r)
  arcs <- vector("list", 6)
  for (i in 0:5) {
    p_in  <- corners[i %% 6 + 1L, ]     # previous corner
    p_cor <- corners[(i + 1) %% 6 + 1L, ]  # corner being rounded
    p_out <- corners[(i + 2) %% 6 + 1L, ]  # next corner
    v1 <- p_in - p_cor
    v2 <- p_out - p_cor
    v1 <- v1 / sqrt(sum(v1^2))
    v2 <- v2 / sqrt(sum(v2^2))
    a0 <- atan2(v1[2], v1[1])
    th <- seq(a0, a0 + pi / 3, length.out = n_arc)
    arcs[[i + 1L]] <- cbind(p_cor[1] + rr * cos(th), p_cor[2] + rr * sin(th))
  }
  do.call(rbind, arcs)
}

out <- file.path("man", "figures", "logo.png")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

png(out, width = 2000, height = 2000, res = 500)
par(mar = c(0, 0, 0, 0), bg = "transparent")
plot.new()
plot.window(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), asp = 1)

# --- radial gradient background (deep navy edge -> teal core) ----------------
edge <- "#081C33"
core <- c("#0E3A5C", "#1282A2")
pal <- colorRampPalette(c(edge, core[1], core[2]))(120)
for (i in seq_along(pal)) {
  r <- 1.00 * (1 - (i - 1) / length(pal))
  polygon(rounded_hex(0, 0, r, round_r = 0.16 * 1.00), col = pal[i], border = NA)
}

# --- subtle inner rim --------------------------------------------------------
polygon(rounded_hex(0, 0, 0.985), col = NA, border = "#FFFFFF40", lwd = 5)

# --- shield (sentry / guardrail motif) ---------------------------------------
shx <- c(-0.31, -0.31, -0.13,  0.00,  0.13,  0.31,  0.31)
shy <- c( 0.46,  0.10, -0.24, -0.33, -0.24,  0.10,  0.46)
shx <- shx * 0.92
shy <- shy * 0.92 + 0.10
polygon(shx, shy, col = "#FFFFFF12", border = "#F2F7FB", lwd = 7,
        ljoin = "round")

# guardrail notches across the shield's top edge
notch_x <- seq(-0.16, 0.16, length.out = 7)
segments(notch_x, rep(0.555, 7), notch_x, rep(0.555, 7) + 0.07,
         col = "#F2F7FB", lwd = 4, lend = "butt")

# --- RNA strand (A/C/G/U nodes) inside the shield ----------------------------
base_cols <- c(A = "#4EC9B0", C = "#F45B69", G = "#FFC53D", U = "#64D2FF")
nx <- c(-0.19, -0.07, 0.05, 0.17)
ny <- c(0.20, 0.335, 0.18, 0.335)
lines(nx, ny, col = "#FFFFFFB3", lwd = 5)
points(nx, ny, pch = 21, cex = 4.2, bg = base_cols, col = "#FFFFFF", lwd = 2.5)

# --- wordmark ----------------------------------------------------------------
text(0, -0.52, "rnaSentry", col = "#FFFFFF", cex = 2.6, font = 2, family = "sans")
text(0, -0.715, "guarded & auditable", col = "#B8DCE6", cex = 1.05,
     font = 1, family = "sans")

dev.off()
message("Wrote ", out)
