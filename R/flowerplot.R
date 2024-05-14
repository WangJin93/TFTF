#' @title Plotting folower plot
#' @description
#' Plotting the folower plot for visualizing the intersection of more than 5 datasets.
#' @import RColorBrewer
#' @param flower_dat List of data.
#' @param ellipse_col_pal Color palette in RColorBrewer, default "Set1".
#' @param label_text_cex Font size, default 1.
#' @export
#'
flowerplot <- function(flower_dat, angle = 90,
                       a = 1, b = 2, r = 1,
                       ellipse_col_pal = "Set1",
                       circle_col = "white",
                       label_text_cex = 1)
{
  set_name <- names(flower_dat)
  item_id <- unique(flower_dat[[1]])
  item_id <- item_id[item_id != '']
  core_item_id <- item_id
  item_num <- length(item_id)

  for (i in 2:length(flower_dat)) {
    item_id <- unique(flower_dat[[i]])
    item_id <- item_id[item_id != '']
    core_item_id <- intersect(core_item_id, item_id)
    item_num <- c(item_num, length(item_id))
  }
  core_num <- length(core_item_id)

  graphics::par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(2,2,2,2))
  graphics::plot(c(0,9), c(0,9), type='n')
  n   <- length(set_name)
  # set the angle of degress
  deg <- 360 / n
  # set the ellipse filling color
  colors <- RColorBrewer::brewer.pal(8, ellipse_col_pal)
  ellipse_col <- grDevices::colorRampPalette(colors)(n)

  res <- lapply(1:n, function(t){
    plotrix::draw.ellipse(x = 5 + cos((angle + deg * (t - 1)) * pi / 180),
                          y = 5 + sin((angle + deg * (t - 1)) * pi / 180),
                          col = ellipse_col[t],
                          border = ellipse_col[t],
                          a = a, b = b,
                          angle = deg * (t - 1))
    graphics::text(x = 5 + 2.5 * cos((angle + deg * (t - 1)) * pi / 180),
                   y = 5 + 2.5 * sin((angle + deg * (t - 1)) * pi / 180),
                   item_num[t],
                   cex = label_text_cex)

    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      graphics::text(x = 5 + 3.3 * cos((angle + deg * (t - 1)) * pi / 180),
                     y = 5 + 3.3 * sin((angle + deg * (t - 1)) * pi / 180),
                     set_name[t],
                     srt = deg * (t - 1) - angle,
                     adj = 1,
                     cex = label_text_cex
      )
    } else {
      graphics::text(x = 5 + 3.3 * cos((angle + deg * (t - 1)) * pi / 180),
                     y = 5 + 3.3 * sin((angle + deg * (t - 1)) * pi / 180),
                     set_name[t],
                     srt = deg * (t - 1) + angle,
                     adj = 0,
                     cex = label_text_cex
      )
    }
  })
  plotrix::draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
  graphics::text(x = 5, y = 5, paste('Overlap:\n', core_num), cex = label_text_cex)
}
