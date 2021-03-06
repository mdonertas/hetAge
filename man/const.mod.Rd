% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{const.mod}
\alias{const.mod}
\title{Construct model matrix}
\usage{
const.mod(covs)
}
\arguments{
\item{covs}{a list of covarietes}
}
\value{
a dataframe that can be used as a design matrix in a model
}
\description{
Using a list of covarietes, constructs a model matrix to be used in a model
}
\examples{
mycovs <- list(batch = sample(c(0, 1), 6, replace = TRUE),
array = as.factor(sample(c('a', 'b', 'c'), 6, replace = TRUE)),
day = sample(c(0, 1, 3), 6, replace = TRUE))
des_mat <- const.mod(mycovs)
print(mycovs)
print(des_mat)

}
