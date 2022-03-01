source("./main.R")
# Read files
path="./data"
fqF <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fqR <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

params <- optimize_dada2_qc(fqF, fqR, 1:10, 273, 175)
# Visualize results
plot_data <- pivot_longer(params, 
                          cols = c("fo_diff", "pass_ratio"),
                          names_to = "names",
                          values_to = "values"
)
plot_data$names <- 
  factor(plot_data$names, levels = c("pass_ratio", "fo_diff"),
         labels = c("Remained read ratio", "First-order differences"))


p <- 
  ggplot(data = plot_data,
            aes(x = ranking_index, y = values)) +
  geom_vline(aes(xintercept = params[which.max(params$fo_diff), "ranking_index"]),
             linetype = "dashed", color = "red", size = 1) +
  geom_text(data = data.frame(names = factor("Remained read ratio", levels = c("Remained read ratio", "First-order differences"))),
            aes(x = params[1, "ranking_index"] + 0.5, y = params[1, "pass_ratio"]),
            label = paste0("Trim position\n",
                           "Forward: ", params$lenF[1],"; Reverse: ", params$lenR[1], "\n",
                           "Max expected errors\n",
                           "Forward: ", params$eeMaxF[1],"; Reverse: ", params$eeMaxR[1], "\n"),
            hjust = 0, vjust = 1, size = 3) +
  labs(x = "Ranking index") + 
  geom_line(size = 1) +
  facet_wrap(~names, nrow = 2, scales = "free_y", strip.position = "left", shrink = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold", color = "black"),
        strip.placement = "outside")

ggsave(p, filename = "./example.jpg", dpi = 300, units = "cm", width = 16, height = 12)
