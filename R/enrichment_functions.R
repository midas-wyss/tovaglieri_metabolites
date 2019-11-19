#' Functions for pathway enrichment
#' These functions will do a pathway enrichment based on tf-idf, similar to the rPEGEOS package I have developed (https://github.com/diogocamacho/rpegeos).  Since that package is still under some development, the functions here offer more stability.
cosine_similarity <- function(x, y, tfidf_crossprod_mat)
{
  # cs <- crossprod(x,y)/sqrt(crossprod(x) * crossprod(y))
  
  # x1 <- apply(x,1,crossprod)
  x2 <- tcrossprod(y) # <-- y is a 1xG vector!!
  x3 <- sqrt(tfidf_crossprod_mat * as.vector(x2))
  
  y1 <- x %*% t(y)
  
  cs <- y1 / x3
  
  return(cs)
}

crossprod_matrix <- function(tfidf_matrix)
{
  cpm <- apply(tfidf_matrix,1,crossprod) # <-- crossproduct matrix for tfidf
  return(cpm)
}

random_sets <- function(number_sets,universe_size,gene_set_size)
{
  # rnd_set <- matrix(0,nrow=ncol(target_tfidf),ncol=number_sets)
  # universe_size = number of columns in tfidf matrix
  # gene_set_size = number of genes that we obtained as differntially expressed
  
  x1 <- matrix(0,nrow=number_sets,ncol=universe_size)
  x2 <- t(replicate(number_sets,sample(universe_size,size=gene_set_size,replace=FALSE)))
  # yy <- apply(x2,1,function(s) {x1[s] <- 1})
  
  for(i in seq(1,number_sets))
  {
    # a1 <- matrix(0,ncol=1,nrow=universe_size)
    # a1[sample(universe_size,size = gene_set_size,replace = FALSE)] <- 1
    # x <- cbind(x,a1)
    x1[i,x2[i,]] <- 1
  }
  
  # rnd_set <- replicate(n=number_sets,sample(x=ncol(target_tfidf),size=gene_set_size,replace=FALSE))
  # return(rnd_set)
  return(x1)
}

random_tfidf <- function(random_set,target_tfidf,tfidf_crossprod_mat)
{
  
  rnd_res <- apply(random_set,1,function(x) cosine_similarity(target_tfidf,t(as.matrix(x)),tfidf_crossprod_mat))
  
  return(rnd_res)
  
}

res_stats <- function(true_results,random_results)
{
  res_stats <- sapply(seq_along(true_results), function(x, y, i) length(which(y[i] > x[i, ]))/length(x[i,]), y = as.matrix(true_results) ,x = random_results)
  return(res_stats)
}