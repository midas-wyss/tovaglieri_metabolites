gene_summarizing <- function(expression_data, gene_id) {
  x1 <- unique(gene_id)
  
  summ_data <- vector(mode = "list", length = length(x1))
  
  for (i in seq(1, length(x1))) {
    x2 <- which(gene_id == x1[i])
    if (length(x2) == 1) {
      summ_data[[i]] <- expression_data[i, ]
    } else {
      x3 <- expression_data[x2, ]
      summ_data[[i]] <- colSums(x3)
    }
  }
  
  summ_data <- do.call(what = rbind, args = summ_data)
  
  return(list(gene_id = x1, summarized_expression = summ_data))
}


# x1 <- unique(human_genes$symbol)
# summ <- vector(mode = "list", length = length(x1))
# for (i in seq(1, length(x1))) {
#   print(i)
#   x2 <- which(human_genes$symbol == x1[i])
#   if (length(x2) == 1) {
#     summ[[i]] <- human_counts[i, ]
#   } else {
#     x3 <- human_counts[x2, ]
#     summ[[i]] <- colSums(x3)
#   }
# }