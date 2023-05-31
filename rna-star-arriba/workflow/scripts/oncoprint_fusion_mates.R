library(ComplexHeatmap)

oncoprint_title = snakemake@params[[1]]
fusions= "\n(fusion mates)"
column_title = paste(oncoprint_title, fusions)

table = read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
mat = as.matrix(table)
mat[is.na(mat)]=""
mat = t(as.matrix(mat))

mat=as.matrix(mat)
mat[is.na(mat)] = ""


col = c("FUS_mate" = "blue")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big black
  FUS_mate = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["FUS_mate"], col = NA))
  }
)


heatmap_legend_param = list(title = "Legend", at = c("FUS_mate"), 
                            labels = c("Fusion Mate"))


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
                         row_names_gp = gpar(fontsize = 5.5), heatmap_legend_param=heatmap_legend_param)
  draw(oncoprint, newpage=F)
}
dev.off()
