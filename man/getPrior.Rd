\name{getPrior}

\alias{getPrior}
\alias{getPrior,aldex.clr-method}

\title{getPrior}

\description{
  Returns the prior for the \code{aldex.clr} object.
}

\usage{
  getPriors(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirichlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns a numeric value, matrix, or array.
}

\value{
  A numeric value, matrix, or array representing the prior (pseudocount)
  used.
}

\seealso{
  \code{aldex.clr}
}

\examples{
data(selex)
    #subset for efficiency
    selex <- selex[1201:1600,]
conds <- c(rep("NS", 7), rep("S", 7))
x <- aldex.clr(selex, conds, mc.samples = 2, denom = "all", verbose = FALSE, prior=0.7)
scale.samps <- getPrior(x)
}
