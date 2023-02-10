# this script creates from scratch the toy plot in Fig 4
# illustrates the difference between negative and positive correlations
library(tidyverse)
library(scales)
library(colorspace)
library(Cairo)

toy_plot <- ggplot(data.frame(x=c(-2,2)), aes(x)) +
  stat_function(fun=function(x)2*x+4, geom = "line", 
                colour = "#C2E1C2", xlim = c(-1.85, 0.5), size = 1.8) + 
  stat_function(fun=function(x)-2.3*x + 0.5, geom = "line", 
                colour = "#F7B05B", xlim = c(0.2, -1.7), size = 1.8) +
  geom_vline(xintercept = -0.6, linetype = "dashed", colour = "firebrick", size = 1) + 
  geom_segment(aes(x = -0.57, y = 5, xend = 0.15, yend = 5),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = -0.63, y = 5, xend = -1.4, yend = 5),
               arrow = arrow(length = unit(0.3, "cm"))) +
  annotate(geom = "text", x = -1, y = 4.75, label = "non-essential", size = 4.5) +
  annotate(geom = "text", x = -0.25, y = 4.75, label = "essential", size = 4.5) +
  annotate(geom = "text", x = 0.65, y = 3.5, colour = darken("#C2E1C2", amount = 0.4),
           label = expression("\u2191 expression \n \u2191 sensitivity"), size = 5) +
  annotate(geom = "text", x = -1.67, y = 2.5, colour = darken("#F7B05B", amount = 0.4),
           label = expression("\u2191 expression \n \u2193 sensitivity"), size = 5) +
  scale_x_continuous(limits = c(-2, 1)) +
  labs(x = "Gene inhibition sensitivity", y = "Gene expression/activity") + 
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.text.y = element_blank(), 
        axis.text.x = element_blank())

ggsave(filename = "./plots/pos_vs_neg_toy_plot.pdf", plot = toy_plot, 
       width = 20, height = 10, dpi = 1000, units = "cm", device = cairo_pdf)
