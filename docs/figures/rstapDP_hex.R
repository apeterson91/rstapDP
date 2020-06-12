library(hexSticker)

# model <- readRDS(file = "~/Documents/NDP_CA/server_model/normal_fit.rds")

p <- bendr::plot_cluster_densities(model,pi_threshold=0.01,switch="color")
p <- p + ggplot2::theme_void() + ggplot2::theme(legend.position = 'none') + ggplot2::ggtitle("")
library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Raleway")
## Automatically use showtext to render text for future devices
showtext_auto()

hexSticker::sticker(p, package="bendr", p_size = 7, p_x = 1, p_y = 1.5,
                    s_x=1, s_y=1, s_width=1.3, s_height=1,
                    h_fill = "#f9fcfc", h_color = "#060707",
                    p_color = "black",# p_family = "Raleway",
                    filename="Desktop/bendr_hex.png")
