% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{plTryRead}
\alias{plTryRead}
\title{plTryRead}
\usage{
plTryRead(patterns)
}
\arguments{
\item{patterns}{a vector of regular expressions.}
}
\value{
The content of a tab-delimited file, if any matches the set of expressions, otherwise NULL.
}
\description{
Tries to find and read a file in the working directory, using an ordered list of regular expression and returning the first matching file.
}
\examples{
plTryRead(c("test"))

}

