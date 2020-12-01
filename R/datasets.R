#' @name process
#' @aliases process_data
#' @title Cheese production process data.
#'
#' @description A dataset containing process data for 4 connected process (A,B,C,D) stages,
#' three unconnected inputs (A,B,C) and a set of endpoint qualities (E). 
#' The stages A, B and C are influencing D and A, B, C and D are influencing E. 
#' This can be described as a directed asyclic graph (sketched below). \cr
#' 
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} B --- }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \ \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \ }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} --+-- \Sexpr{"\u200B"} | }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} /\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} |\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \| }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} A --- D -- E }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} /\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} / }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} C --- }\cr
#' 
#' Subpaths include: ADE, AE, BDE, BE, CDE, and CE, and SO-PLS-PM models model
#' A+B+C->D and A+B+C+D->E.
#'
#' @docType data
#' @usage data(process)
#' 
#' @format A list of matrices having 795 rows and varying numbers of variables:
#' \describe{
#'   \item{A}{Process stage A}
#'   \item{B}{Process stage B}
#'   ...
#' }
NULL

#' @name simulated
#' @title Data simulated to have certain characteristics.
#'
#' @description A dataset containing simulated data for 4 connected events where A is the
#' starting point and D is the end point. This can be described as a directed
#' asyclic graph (sketched below, moving left->right). \cr
#' 
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} C    }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} /^\   }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} / \Sexpr{"\u200B"}| \  }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} A--+->D }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \ \Sexpr{"\u200B"}| /  }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \|/   }\cr
#' \code{\Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} B    }\cr
#' 
#' Subpaths include: ABD, AD, ABCD, ACD
#'
#' @docType data
#' @usage data(simulated)
#' 
#' @format A list of matrices having 200 rows and 10 variables:
#' \describe{
#'   \item{A}{Simulated matrix A}
#'   \item{B}{Simulated matrix B}
#'   ...
#' }
#' 
#' @references Næs, Romano, ..., Liland, In preparation.
NULL

#' @name candies
#' @title Sensory assessment of candies.
#'
#' @description A dataset containing 9 sensory attributes for 5 candies assessed
#' by 11 trained assessors.
#'
#' @docType data
#' @usage data(candies)
#' 
#' @format A data.frame having 165 rows and 3 variables:
#' \describe{
#'   \item{assessment}{Matrix of sensory attributes}
#'   \item{assessor}{Factor of assessors}
#'   \item{candy}{Factor of candies}
#' }
#' 
#' @references Luciano G, Næs T. Interpreting sensory data by combining principal 
#' component analysis and analysis of variance. Food Qual Prefer. 2009;20(3):167‐175.
NULL
