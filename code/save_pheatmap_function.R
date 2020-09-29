# Usage: 
# 
# Step1: create heatmap 
# my_heatmap <- pheatmap(df_for_heatmap)  
#
# Step 2: save it
# save_pheatmap_png(my_heatmap, "heatmap.png")

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
