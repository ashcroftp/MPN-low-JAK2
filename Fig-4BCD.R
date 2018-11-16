#' Fig-4BCD.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulateModel.R")
source("R/plotFunctions.R")

#' Linear population structure -- how do changes in the mutant's birth rate,
#' death probability, or self-renewal probability affect the clonal burden.

#' Parameter definitions for this section
n <- 6
adj <- createAdjacency(n)
#' Initial conditions
wt.ic <- c(0.9, rep.int(0, n - 1))
mut.ic <- c(0.1, rep.int(0, n - 1))
#' Wild-type cell parameters
wt.div <- c(seq(from = 1.0, to = 3.0, length.out = n - 1), 1.0)
wt.sr <- c(0.5, rep.int(0.4, n - 2), 0.0)
wt.adj <- adj
wt.death <- c(0, rep.int(0.2, n - 1))


#' Plot for faster dividing mutant cells in a non-HSC compartment,
#' resulting in a lower clonal fraction in the modified compartment.
#' All other fractions are preserved.
#' Mutant cell parameters
mut.div <- wt.div
mut.div[4] <- 2 * wt.div[4]
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.death <- wt.death
#' Steady state and clonal fraction
steadyState <- numericSteadyState(
    wt.ic = wt.ic, mut.ic = mut.ic,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death
)
clonal.fraction <- getClonalFraction(steadyState)
#' Plot the data
plotFractionParams(
    fraction.df = clonal.fraction,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death,
    show.var = c("div"),
    colour.vec = rep.int("#FA2D1C", n), alpha.vec = seq.int(0.2, 1, length.out = n),
    show.tree = T, y.pos = c(1,1,1,1,1,1), tree.height = 0.2, y.lim = 0.5
)
#savePlot(filename="4B-left.pdf", width = 60, height.mult = 1)

#' Plot for faster dividing cells from HSC through to penultimate
#' Mutant cell parameters
mut.div <- 2 * wt.div
mut.div[6] <- wt.div[6]
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.death <- wt.death
#' Steady state and clonal fraction
steadyState <- numericSteadyState(
    wt.ic = wt.ic, mut.ic = mut.ic,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death
)
clonal.fraction <- getClonalFraction(steadyState)
#' Plot
plotFractionParams(
    fraction.df = clonal.fraction,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death,
    show.var = c("div"),
    colour.vec = rep.int("#FA2D1C", n), alpha.vec = seq.int(0.2, 1, length.out = n),
    show.tree = T, y.pos = c(1,1,1,1,1,1), tree.height = 0.2, y.lim = 0.5
)
#savePlot(filename="4B-right.pdf", width = 60, height.mult = 1)


#' Plot less-death-prone cells early
#' All other fractions are preserved.
#' Mutant cell parameters
mut.div <- wt.div
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.death <- wt.death
mut.death[2] <- 0.05

#' Steady state and clonal fraction
steadyState <- numericSteadyState(
    wt.ic = wt.ic, mut.ic = mut.ic,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death
)
clonal.fraction <- getClonalFraction(steadyState)
#' Plot the data
plotFractionParams(
    fraction.df = clonal.fraction,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death,
    show.var = c("death"),
    colour.vec = rep.int("#FA2D1C", n), alpha.vec = seq.int(0.2, 1, length.out = n),
    show.tree = T, y.pos = c(1,1,1,1,1,1), tree.height = 0.2, y.lim = 0.5
)
#savePlot(filename="4C-left.pdf", width = 60, height.mult = 1)

#' Plot less-death-prone cells late
#' All other fractions are preserved.
#' Mutant cell parameters
mut.div <- wt.div
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.death <- wt.death
mut.death[n - 1] <- 0.05

#' Steady state and clonal fraction
steadyState <- numericSteadyState(
    wt.ic = wt.ic, mut.ic = mut.ic,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death
)
clonal.fraction <- getClonalFraction(steadyState)
#' Plot the data
plotFractionParams(
    fraction.df = clonal.fraction,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death,
    show.var = c("death"),
    colour.vec = rep.int("#FA2D1C", n), alpha.vec = seq.int(0.2, 1, length.out = n),
    show.tree = T, y.pos = c(1,1,1,1,1,1), tree.height = 0.2, y.lim = 0.5
)
#savePlot(filename="4C-right.pdf", width = 60, height.mult = 1)

#' Plot for increased self-renewal in a single compartment
#' Mutant cell parameters
mut.div <- wt.div
mut.sr <- wt.sr
mut.sr[3] <- 0.49
mut.adj <- wt.adj
mut.death <- wt.death
#' Steady state and clonal fraction
steadyState <- numericSteadyState(
    wt.ic = wt.ic, mut.ic = mut.ic,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death
)
clonal.fraction <- getClonalFraction(steadyState)
#' Plot
plotFractionParams(
    fraction.df = clonal.fraction,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death,
    show.var = c("sr"),
    colour.vec = rep.int("#FA2D1C", n), alpha.vec = seq.int(0.2, 1, length.out = n),
    show.tree = T, y.pos = c(1,1,1,1,1,1), tree.height = 0.2, y.lim = 0.5
)
#savePlot(filename="4D-left.pdf", width = 60, height.mult = 1)

#' Plot for greater self-renewal in an amplifying compartment
#' Mutant cell parameters
mut.div <- wt.div
mut.sr <- wt.sr
mut.sr[n - 1] <- 0.49
mut.adj <- wt.adj
mut.death <- wt.death
#' Steady state and clonal fraction
steadyState <- numericSteadyState(
    wt.ic = wt.ic, mut.ic = mut.ic,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death
)
clonal.fraction <- getClonalFraction(steadyState)
#' Plot
plotFractionParams(
    fraction.df = clonal.fraction,
    wt.div = wt.div, mut.div = mut.div,
    wt.sr = wt.sr, mut.sr = mut.sr,
    wt.adj = wt.adj, mut.adj = mut.adj,
    wt.death = wt.death, mut.death = mut.death,
    show.var = c("sr"), 
    colour.vec = rep.int("#FA2D1C", n), alpha.vec = seq.int(0.2, 1, length.out = n),
    show.tree = T, y.pos = c(1,1,1,1,1,1), tree.height = 0.2, y.lim = 0.5
)
#savePlot(filename="4D-right.pdf", width = 60, height.mult = 1)
