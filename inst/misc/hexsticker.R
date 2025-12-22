col <- list(trueblue=rgb(0/255,25/255,101/255),
            lightblue=rgb(59/255,151/255,222/255),
            oceangreen=rgb(42/255,145/255,139/255),
            rosepink=rgb(238/255,167/255,191/255))

## install.packages(c("rsvg", "svglite", "systemfonts", "hexSticker", "meme"))
## meme::font_pokemon()
library("hexSticker")
library("extrafont")
library("showtext")
font_add_google(name = "Pacifico", family = "pacifico")
font_add_google(name = "Orbitron", family = "orbitron")
imgurl <- "../../man/figures/timerrr.png"
col1 <- col[[2]] # "#fcb100"
fntcol <- lava::Col(col[[2]], 1) # col1 # "gray"
bgcol <- "#fffaed"
s <- sticker(imgurl,
  package = "mets",
  s_x = .97, s_y = .9, s_width = .6,
  p_color = fntcol,
  p_size = 75,
  p_y = 1.5,
  p_family = "orbitron",
  white_around_sticker = FALSE,
  h_fill = gray(0.95),
  h_color = "darkblue",
  h_size = 2,
  dpi = 300,
  filename = "logo.jpg",
  ## asp = 2,
  )
k <- 50
ggplot2::ggsave(
  "logo.png",
  s,
  #+
  ## ggplot2::theme(plot.margin = ggplot2::margin(k, k, k, k), "pt")
, scale=0.6
)
## ggplot2::ggsave("logo.png", s)
p <- magick::image_read("logo.png")
p2 <- magick::image_scale(p, "500")
magick::image_write(image = p2, path = "../../man/figures/logo.png")

## system("magick logo.png -transparent white logo2.png")

library(magick)
p <- magick::image_read("logo.png")
xx <- c(1,2,image_info(p)$width-2,image_info(p)$width-1)
yy <- c(1,2,image_info(p)$height-2,image_info(p)$height-1)
for (x in xx)
  for (y in yy)
    p <- magick::image_fill(p, "transparent", point=paste0("+",x,"+",y))
image_write(image = p, path = "../../man/figures/logohex1.png")
p2 <- magick::image_scale(p, "500")
image_write(image = p2, path = "../../man/figures/logohex.png")
