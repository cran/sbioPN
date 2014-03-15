\name{b) Simulation Functions}
\alias{sGillespieOptimDirect}
\alias{sGillespieDirectCR}
\alias{sRungeKuttaDormandPrince45}
\alias{sHaseltineRawlings}
\title{Simulation of a biochemical system with spatial effects}
\description{These functions simulate a biochemical reacton system
  with spatial effects parameterized as a Petri Net. 
  \code{sGillespieOptimDirect} and \code{sGillespieDirectCR} performs pure
  stochastic simulations, \code{sRungeKuttaDormandPrince45} a pure
  deterministic integration, \code{sHaseltineRawlings} a hybrid of the
  above. Multiple runs can be performed at once.

  See \code{\link{init}} for a way of defining the model that is close
  to the way reactions are written.
}
\usage{
## Exact stochastic simulation: 
sGillespieOptimDirect(model, timep, delta=1, runs=1)
sGillespieDirectCR(model, timep, delta=1, runs=1)

## Pure deterministic:
sRungeKuttaDormandPrince45(model, timep, delta=1, ect = 1e-09)

## Hybrid stochastic/deterministic:
sHaseltineRawlings(model, timep, delta=1, runs=1, ect = 1e-09)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{list containing named elements:}
  \item{timep}{It can be either a numeric, indicating for how long (in
    the same time units as the propensity constants) the process will
    run, or a functions (R or C), in which case can be used to change
    the protocol at time intervals. See details.}
  \item{delta}{Interval time at which the state will be saved.}
  \item{runs}{How many runs will be performed.}
  \item{ect}{Precision for the fast reactions.}
}
\details{
  model is a list containing the following elements:
  \itemize{
    \item{model$pre}{: pre matrix, with as many rows as transitions
      (reactions), and columns as places (reactants). It has
      the stoichiometrics of the left sides of the reactions.}
    \item{model$post}{: post matrix, with as many rows as transitions,
      and columns as places (products). It has the stoichiometrics of
      the right sides of the reactions.}
    \item{model$h}{: list of propensity constants or functions returning
      the propensity (with as many elements as transitions).}
    \item{model$slow}{: vector of zeros for slow transitions and ones
      for fast transitions. Only needed for
      \code{HaseltineRawlings}. Ignored otherwise.}
    \item{model$M}{: initial marking (state) of the system.}
    \item{model$place}{: vector with names of the places.}
    \item{model$transition}{: vector with names of the transitions.}
  }
}
\value{
  The functions return a list with the following elements:  
  \item{place}{vector with the names of the places if supplied. If not, the function creates names as follows: P1, P2, ...}
  \item{transition}{vector with the names of the transitions if supplied. If not, the function creates names as follows: T1, T2, ...}
  \item{dt}{vector containing the discretized times at which the state is saved (according to delta)}
  \item{run}{list with as many elements as runs. We will describe the first element, run[[1]], as the rest have exactly the same structure. It is also a list, with the following elements:}
  \item{run[[1]]$M}{list with as many elements as places, each of them containing the state of the system sampled according to delta.}
  \item{run[[1]]$transitions}{vector with as many elements as transitions, with the total of time each slow reaction fired.}
  \item{run[[1]]$tot.transitions}{numeric with the summ of run[[1]]$transitions.}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{init}}, \code{\link{atr}}
}
\examples{
## sbioPN has been tested only on 64 bits machines.
## It may fail in 32 bits architecture.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
