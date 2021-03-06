% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_methods.R
\name{getPredMort}
\alias{getPredMort}
\title{get predation mortality rate}
\usage{
getPredMort(object, n, n_pp, pred_rate, time_range, drop = TRUE)
}
\arguments{
\item{object}{A \code{MizerParams} or \code{MizerSim} object.}

\item{n}{A matrix of species abundance (species x size). Only used if
\code{object} argument is of type \code{MizerParams}.}

\item{n_pp}{A vector of the plankton abundance by size. Only used if
\code{object} argument is of type \code{MizerParams}.}

\item{pred_rate}{An array of predation rates of dimension no. sp x no.
community size bins x no. of size bins in whole spectra (i.e. community +
plankton, the w_full slot). The array is optional. If it is not provided
it is calculated by the \code{getPredRate()} method.}

\item{time_range}{Subset the returned fishing mortalities by time. The time
range is either a vector of values, a vector of min and max time, or a
single value. Default is the whole time range. Only used if the
\code{object} argument is of type \code{MizerSim}.}

\item{drop}{Only used when object is of type \code{MizerSim}. Should
dimensions of length 1 in the output be dropped, simplifying the output.
Defaults to TRUE}
}
\value{
If a \code{MizerParams} object is passed in, the method returns a two
  dimensional array (prey species x prey size) based on the abundances also
  passed in. If a \code{MizerSim} object is passed in, the method returns a
  three dimensional array (time step x prey species x prey size) with the
  predation mortality calculated at every time step in the simulation.
}
\description{
Calculates the total predation mortality rate \eqn{\mu_{p,i}(w_p)} on each
prey species by prey size. This method is used by the \code{\link{project}}
method for performing simulations.
}
\examples{
\dontrun{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
# With constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get predation mortality at one time step
n <- sim@n[21,,]
n_pp <- sim@n_pp[21,]
getPredMort(params,n,n_pp)
# Get predation mortality at all saved time steps
getPredMort(sim)
# Get predation mortality over the time 15 - 20
getPredMort(sim, time_range = c(15,20))
}
}
\seealso{
\code{\link{getPredRate}} and \code{\link{project}}.
}
