library(ComplexHeatmap)

column_title = "OncoPrint CCS"

in_tab="/Users/dawid/sciebo/drArbeit/diverse/data/arriba/matrix_filtriert.tsv.txt"
pdf_out="/Users/dawid/sciebo/drArbeit/diverse/data/arriba/test.pdf"




table = read.table(in_tab, sep="\t", header=TRUE, row.names=1)
mat = as.matrix(table)
mat[is.na(mat)]=""
mat = t(as.matrix(mat))

mat=as.matrix(mat)
mat[is.na(mat)] = ""


col = c("translocation" = "red", 
        "deletion_read_through" = "blue",
        "inversion" = "#9ACD32",
        "duplication" = "#CD9B1D",
        "deletion" = "#1C86EE",
        "translocation_truncation"="gray40",
        "inversion_truncation"="gray40",
        "known" = "#893DF6", 
        "duplication_truncation"="gray40",
        "deletion_read_through_truncation"="gray40",
        "deletion_truncation"="gray40",
        "duplication_non_canonical_splicing" = "#FFC125",
        "duplication_ITD" = "#DAA520")

# translocation                         316
# deletion/read-through                 188
# inversion                             110
# duplication                            99
# deletion                               84
# translocation_truncation               78
# inversion_truncation                   54
# known                                  38
# duplication_truncation                 37
# deletion/read-through_truncation       28
# deletion_truncation                    22
# duplication/non-canonical_splicing     11
# duplication/ITD                        10

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#FFFFFF", col = NA))
  },
  # known
  known = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["known"], col = NA))
  },
  # translocation
  translocation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["translocation"], col = NA))
  },
  # deletion
  deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["deletion"], col = NA))
  },
  # duplication_non_canonical_splicing
  duplication_non_canonical_splicing = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["duplication_non_canonical_splicing"], col = NA))
  },
  # deletion_read_through
  deletion_read_through = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["deletion_read_through"], col = NA))
  },
  # duplication
  duplication = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["duplication"], col = NA))
  },
  # inversion
  inversion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["inversion"], col = NA))
  },
  # duplication_ITD
  duplication_ITD = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["duplication_ITD"], col = NA))
  },
  # translocation_truncation
  translocation_truncation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["translocation_truncation"], col = NA))
  },
  # inversion_truncation
  translocation_truncation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["inversion_truncation"], col = NA))
  },
  # duplication_truncation
  translocation_truncation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["duplication_truncation"], col = NA))
  },
  # deletion/read-through_truncation
  deletion_read_through_truncation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["deletion_read_through_truncation"], col = NA))
  },
  # inversion_truncation
  inversion_truncation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["inversion_truncation"], col = NA))
  },
  # deletion_truncation
  deletion_truncation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["deletion_truncation"], col = NA))
  },
  # duplication_truncation
  duplication_truncation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["duplication_truncation"], col = NA))
  }
)


heatmap_legend_param = list(title = "Legend", at = c("known",
                                                     "translocation",
                                                     "deletion_read_through",
                                                     "deletion",
                                                     "duplication",
                                                     "duplication_non_canonical_splicing",
                                                     "duplication_ITD",
                                                     "inversion",
                                                     "translocation_truncation"), 
                            labels = c("known fusion", 
                                       "translocation", 
                                       "deletion/read-through",
                                       "deletion","duplication",
                                       "duplication_non_canonical_splicing",
                                       "duplication_ITD",
                                       "inversion",
                                       "truncation"))


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
pdf(file = pdf_out, height=height_plot)
if (rows_matrix > 0) {
  oncoprint <- oncoPrint(mat,
                         alter_fun = alter_fun, col = col,
                         remove_empty_columns = FALSE, remove_empty_rows = TRUE, 
                         column_title = column_title,show_column_names=T, 
                         row_names_gp = gpar(fontsize = 3),heatmap_legend_param=heatmap_legend_param, column_names_side="top",
                         column_names_gp = gpar(fontsize = 5),show_heatmap_legend=TRUE)
  draw(oncoprint, newpage=F)
}
dev.off()


