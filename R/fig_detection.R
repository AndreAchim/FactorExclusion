# fig_detection.R — Observed and fitted detection rates (quasibinomial GLM)
#
# Prerequisites: df_det and m3q must exist in the environment.
#   source("R/build_detection_data.R")   # builds df_det, m0q, m3q from efg
# Usage:
#   source("R/fig_detection.R")   # produces Figure2.pdf and Figure2.png
#
# Note: curves are from the quasibinomial GLM (m3q); the formal test
# exploiting the paired structure of simulations uses clogit (cl0-cl2).
# Duplicate att values simply contribute additional paired cases to the
# regression; the fitted curves already summarise them.

library(ggplot2)

df_det$fitted <- fitted(m3q)
df_det$prop   <- df_det$n_det / df_det$n_tot
df_det$N_lab  <- factor(paste0("N = ", df_det$N_val),
                        levels = paste0("N = ", c(125, 500, 2000)))

# 5 saturation levels: grey shades + line types for B&W legibility
greys  <- c("(0.4,0.5)" = "grey75", "(0.5,0.6)" = "grey58",
            "(0.6,0.7)" = "grey42", "(0.7,0.8)" = "grey25",
            "(0.8,0.9)" = "grey8")
ltypes <- c("(0.4,0.5)" = "dotted",  "(0.5,0.6)" = "longdash",
            "(0.6,0.7)" = "dashed",  "(0.7,0.8)" = "dotdash",
            "(0.8,0.9)" = "solid")
shapes <- c("(0.4,0.5)" = 1, "(0.5,0.6)" = 2, "(0.6,0.7)" = 0,
            "(0.7,0.8)" = 6, "(0.8,0.9)" = 16)

fig_det <- ggplot(df_det, aes(x = att, color = sat, linetype = sat, shape = sat)) +
  geom_point(aes(y = prop), size = 0.9, alpha = 0.6,
             position = position_jitter(width = 0.002, seed = 4721)) +
  geom_line(aes(y = fitted), linewidth = 0.5) +
  facet_wrap(~ N_lab, ncol = 1) +
  scale_color_manual(values = greys,    name = "Anchor\nloadings") +
  scale_linetype_manual(values = ltypes, name = "Anchor\nloadings") +
  scale_shape_manual(values = shapes,   name = "Anchor\nloadings") +
  guides(color    = guide_legend(reverse = TRUE),
         linetype = guide_legend(reverse = TRUE),
         shape    = guide_legend(reverse = TRUE)) +
  labs(x = "|Expected residual correlation|",
       y = expression(paste("Detection rate (", alpha, " = .05)"))) +
  theme_bw(base_size = 9) +
  theme(legend.position  = "right",
        legend.key.width  = unit(1.0, "cm"),
        legend.text       = element_text(size = 7),
        legend.title      = element_text(size = 8),
        panel.grid.minor  = element_blank(),
        strip.background  = element_rect(fill = "grey90", color = NA),
        strip.text        = element_text(colour = "black"))

ggsave("Figure2.pdf", fig_det, width = 12, height = 18, units = "cm")
ggsave("Figure2.png", fig_det, width = 12, height = 18, units = "cm", dpi = 300)
message("Figure2.pdf et Figure2.png enregistrés.")
