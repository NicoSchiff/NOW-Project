# Depth-layer contrasts within each station -------------------------------
# ──────────────────────────────────────────────────────────────
# 1)  VERTICAL STRUCTURE: Dunn tests between depth strata
# ──────────────────────────────────────────────────────────────
Strate_dunn_results <- NOW %>% # raw data
  select(Station_Cluster, eventDate, classification, any_of(var_levels)) %>% # keep what we need
  pivot_longer(all_of(var_levels), # → long format
               names_to  = "variable",
               values_to = "value"
  ) %>%
  # filter(!is.na(value)) %>%                                # Dunn requires complete cases
  group_by(variable, Station_Cluster) %>% # one test per variable * cluster
  rstatix::dunn_test(value ~ classification, # Dunn on ranks
                     p.adjust.method = "bonferroni"
  ) %>% # Bonferroni within each group
  ungroup()

# Heat-map of pair-wise strata differences
Strate_dunn.plot <- ggplot(
  Strate_dunn_results,
  aes(group1, group2, fill = p.adj.signif)
) +
  geom_tile(colour = "grey90", linewidth = 0.3) +
  facet_grid(Station_Cluster ~ variable, labeller = label_wrap_gen()) +
  scale_fill_manual(
    values = c(rev(pals::ocean.deep(5))[1:4], "grey75"),
    drop = FALSE,
    name = "Bonferroni\nsignificance"
  ) +
  theme(
    axis.text.x = element_text(
      size = 10, colour = "black",
      angle = 45, hjust = 1
    ),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

Strate_dunn.plot
#
# ggsave(
#   plot = Strate_dunn.plot,
#   filename = "Depth_contrasts.png", path = "DataBase_Nut/20250603/Figure/",
#   width = 12.5, height = 12.5, dpi = 300
# )

# ──────────────────────────────────────────────────────────────
# 2)  LATERAL STRUCTURE: Dunn tests between stations
# ──────────────────────────────────────────────────────────────
Station_dunn_results <- NOW %>%
  select(Station_Cluster, eventDate, classification, any_of(var_levels)) %>%
  pivot_longer(all_of(var_levels),
               names_to  = "variable",
               values_to = "value"
  ) %>%
  # filter(!is.na(value)) %>%                                # Dunn requires complete cases
  group_by(variable, classification) %>% # one test per variable * depth
  rstatix::dunn_test(value ~ Station_Cluster,
                     p.adjust.method = "bonferroni"
  ) %>%
  ungroup()

# Prepare the tile data
# Station_tile_dat <- Station_dunn_results %>%
#   mutate(
#     group2   = factor(group2, levels = rev(St_order)), # reverse y axis
#     group1   = factor(group1, levels = St_order), # x axis
#     variable = factor(variable, levels = var_levels), # keep canonical order
#     Var.plot = factor(var_labels[as.character(variable)], levels = var_labels)
#   )
#
# # Heat-map of pair-wise station differences
Station_dunn.plot <- ggplot(
  Station_dunn_results,
  aes(group1, group2, fill = p.adj.signif)
) +
  geom_tile(colour = "grey90", linewidth = 0.3) +
  facet_grid(classification ~ variable, labeller = label_wrap_gen()) +
  scale_fill_manual(
    values = c(rev(pals::ocean.deep(5))[1:4], "grey75"),
    drop = FALSE,
    name = "Bonferroni\nsignificance"
  ) +
  theme(
    axis.text.x = element_text(
      size = 10, colour = "black",
      angle = 45, hjust = 1
    ),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

Station_dunn.plot
#
# ggsave(
#   plot = Station_dunn.plot,
#   filename = "Station_contrasts.png", path = "DataBase_Nut/20250603/Figure/",
#   width = 17.5, height = 10, dpi = 300
# )
