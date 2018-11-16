#' Fig-4F.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulateModel.R")
source("R/plotFunctions.R")

#' Branched population structure -- how do changes in the mutant's differentiation
#' bias affect the clonal burden.

#' Parameter defintions for this section
n <- 6
adj <- createAdjacency(n, entries =
                           list(
                               c(1,2,1),
                               c(2,3,0.5), c(2,5,0.5),
                               c(3,4,1),
                               c(5,6,1)
                           )
)
#' Initial conditions
wt.ic <- c(0.9, rep.int(0, n - 1))
mut.ic <- c(0.1, rep.int(0, n - 1))
#' Wild-type cell parameters
wt.div <- seq(from = 1.0, to = 3, length.out = n)
wt.div[c(4,6)] <- wt.div[1]
wt.sr <- c(0.5, rep.int(0.4, n - 2), 0.0)
wt.adj <- adj
wt.death <- c(0, rep.int(0.1, n - 1))
#' Mutant cell parameters
#' Here we just modify the branching fraction
mut.div <- wt.div
mut.sr <- wt.sr
mut.adj <- wt.adj
mut.adj[2,3] <- 1.5 * wt.adj[2,3]
mut.adj[2,5] <- 1 - mut.adj[2,3]
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
    show.var = c("branch"),
    colour.vec = ifelse(seq_len(n) <= 4, "#FDB12B", "#FA2D1C"), alpha.vec = c(0.4, 0.55, 0.7, 1.0, 0.7, 1.0),
    show.tree = T, y.pos = c(1,1,0,0,1,1), tree.height = 0.3, y.lim = -0.5
)
#savePlot(filename="4E-right.pdf", width = 60, height.mult = 1.3)

