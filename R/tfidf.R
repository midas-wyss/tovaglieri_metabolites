#' Combined tf-idf matrix
#'
#' Returns a tf-idf sparse matrix given a pathway set matrix.
#' 
#' @param data_matrix Sparse document-term matrix (pathway by gene matrix)
#' @return A sparse matrix with the computed tf-idf
tfidf <- function(data_matrix)
{
  # tidy data ----
  dtm1 <- tidytext::tidy(data_matrix)

  # count words in both dtms ----
  words1 <- dtm1 %>% 
    dplyr::count(row, column) %>% 
    dplyr::ungroup()
  
  # count all terms ----
  terms1 <- words1 %>% 
    dplyr::group_by(row) %>% 
    dplyr::summarize(total = sum(n))
  
  # calculate tf-idf ----
  tfidf1 <- words1 %>% 
    tidytext::bind_tf_idf(column, row, n)
  
  m1 <- matrix(0, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  for (i in seq(1, nrow(data_matrix))) {
    a1 <- which(tfidf1$row == rownames(data_matrix)[i])
    a2 <- match(tfidf1$column[a1], colnames(data_matrix))
    m1[i, a2] <- tfidf1$tf_idf[a1]
  }

  colnames(m1) <- colnames(data_matrix)
  rownames(m1) <- rownames(data_matrix)

  # make sparse matrix ----
  M <- Matrix::Matrix(m1, sparse = TRUE)
  
  return(M)
  
}