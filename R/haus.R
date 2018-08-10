#' Hausdorff distance
#'
#' Calculates the Hausdorff distance between two sets.
#'
#' @param A A set.
#' @param B A set.
#' 
#' @return Returnst the Hausdorff distance between the two sets A and B.
#' 
#' @details Given A and B, two non-empty subsets of a metric space \eqn{(M, d)}, their Hausdorff distance \eqn{d_H (A, B)} is:
#' 
#' \deqn{d_{ H}(A,B) = \max{\sup_{x \in A}\inf_{y \in B} d(x,y),\, \sup_{y \in B} \inf_{x \in A} d(x,y)}.}{d_H(A,B) = \max{\sup_{x \in A} inf_{y \in B} d(x,y), sup_{y \in B} inf_{x \in A} d(x,y)}.}
#' 
#' Given that the curves in study are a subsets of \eqn{R^2} , \eqn{\{(t_1 , x_1 ),..., (t_n , x_n )\}}. Hence a natural choice for the metric \eqn{d} is the Euclidean distance on \eqn{R^2} .
#' 
#' @export
hausdorff <- function(A,B)
{
  H <- max(distance(A,B),distance(B,A))
  return(H)
}

eucl.dist <- function(x) {sqrt(sum(x^2))}

distance <- function(A,B)
{
   n <- nrow(A)
   m <- nrow(B)
   a <- numeric(n)
   b <- numeric(m)

   for (i in 1:n)
      {
        for (j in 1:m)
        {
          b[j] <- eucl.dist(A[i,]-B[j,])
        }
          a[i] <- min(b)
      }
   return(max(a))
}
