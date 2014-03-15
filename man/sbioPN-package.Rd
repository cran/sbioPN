\name{a) sbioPN package}
\alias{sbioPN-package}
\docType{package}
\title{Simulation of deterministic and stochastic biochemical reaction
  networks with spatial effects using Petri Nets}
\description{
  sbioPN is a package of C functions that can be used to simulate
  time-dependent evolution of biochemical reaction networks with spatial
  effects. The model
  is defined as a place/transition Petri Net, which is close to how
  biochemical reactions are defined. The model can be either
  deterministically solved using an explicit Runge Kutta Dormand Prince
  45 method, simulated using two highly optimized variants of the
  stochastic simulation algorithm, or as a deterministic/stochastic
  hybrid, according to the Haseltine and Rawlings' algorithm. The library has been optimized for
  speed and flexibility.

  bioPN has been tested only on 64 bits machines, relying on integers of
  64 bits. The behavior on 32 bits architectures is untested and not supported.
}
\details{
\tabular{ll}{
Package: \tab sbioPN\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2014-03-15\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Roberto Bertolusso and Marek Kimmel

Maintainer: Roberto Bertolusso <rbertolusso@rice.edu>
}
