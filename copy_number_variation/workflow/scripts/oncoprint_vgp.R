library(ComplexHeatmap)

column_title = snakemake@params[[1]]

table = read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
mat = as.matrix(table)
mat[is.na(mat)]=""
mat = t(as.matrix(mat))

mat=as.matrix(mat)
mat[is.na(mat)] = ""


col = c("amp1" = "darkgreen", "amp2" = "green4", "amp3" = "green1","LOSS" = "gray34", "HOMDEL" = "black")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#FFFFFF", col = NA))
  },
 # small green
  amp1 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.4, 
              gp = gpar(fill = col["amp1"], col = NA))
  },
  # small green
  amp2 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.66, 
              gp = gpar(fill = col["amp2"], col = NA))
  },
  # small green
  amp3 = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.9, 
              gp = gpar(fill = col["amp3"], col = NA))
  },
  # small green
  LOSS = function(x, y, w, h) {
    grid.polygon(
        unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
        unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
              gp = gpar(fill = col["LOSS"], col = NA))
  },
  # small green
  HOMDEL = function(x, y, w, h) {
    grid.polygon(
        unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
        unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
              gp = gpar(fill = col["HOMDEL"], col = NA))
  } 
)


heatmap_legend_param = list(title = "Legend", at = c("amp1", "amp2", "amp3", "LOSS", "HOMDEL"), 
                            labels = c("amp 3-5", "amp 5-10", "amp >10", "LOSS", "HOMDEL"))
if (ncol(mat) > 1 ){
  mat <- mat[order(apply(mat, 1, function(row) sum(row != "")), decreasing = T), ]
}

if (nrow(mat) > 2000) {
  mat <- mat[1:2000,]
}
rows_matrix <- nrow(mat)
height_plot <- (rows_matrix/5)
if (height_plot < 4) {
  height_plot <- 4
}
pdf(file = snakemake@output[[1]], height=height_plot)
if (rows_matrix > 0) {
  oncoprint <- oncoPrint(mat,
                         alter_fun = alter_fun, col = col,
                         remove_empty_columns = FALSE, remove_empty_rows = TRUE, 
                         column_title = column_title,show_column_names=T, 
                         row_names_gp = gpar(fontsize = 5.5),heatmap_legend_param=heatmap_legend_param)
  draw(oncoprint, newpage=F)
}
dev.off()