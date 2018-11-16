#' Fig-S6.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulateModel.R")
source("R/plotFunctions.R")

#' Replicating specific experimental results.

#' Parameter definitions for this section
name.vec <- c(`1` = "HSC", `2` = "CMP", `3` = "MEP", `4` = "MkP", `5` = "Meg", `6` = "PLT",
              `7` = "EP", `8` = "RET", `9` = "CFU-G", `10` = "GRA")
n <- length(name.vec)
adj <- createAdjacency(n, entries =
                           list(
                               c(1, 2, 1.0), #c(1, 4, 0.1),
                               c(2, 3, 0.5), c(2, 9, 0.5),
                               c(3, 4, 0.5), c(3, 7, 0.5),
                               c(4, 5, 1.0), c(5, 6, 1.0),
                               c(7, 8, 1.0),
                               c(9, 10, 1.0)
                           )
)
y.pos <- c(1,1,0,-1,-1,-1,0,0,1,1)
#colour.vec <- ifelse(seq_len(n) <= 6, "#FDB12B", ifelse(seq_len(n) <= 8, "#FA2D1C", "#919191"))
colour.vec <- ifelse(seq_len(n) <= 3, "#FDB12B", ifelse(seq_len(n) <= 6, "#B97E42", ifelse(seq_len(n) <= 8, "#FA2D1C", "#919191")))

#' Initial conditions
wt.ic <- c(0.99, rep.int(0, n - 1))
mut.ic <- c(0.01, rep.int(0, n - 1))
#' Wild-type cell parameters
# wt.div <- seq(from = 1.0, to = 3, length.out = n)
# wt.div[c(6,9,11)] <- wt.div[1]
# wt.sr <- seq(from = 0.5, to = 0.0, length.out = n)
wt.div <- c(1.0, 1.4, 1.8, 2.2, 2.6, 2.0, 2.2, 1.0, 2.8, 2.0)
wt.sr <- c(0.5, 0.3, 0.3, 0.3, 0.3, 0.0, 0.3, 0.0, 0.3, 0.0)
wt.adj <- adj
wt.death <- c(0, rep.int(0.2, n - 1))
wt.death[7] <- 0.6 # Higher due to lack of EPO


# Gradual platelet expansion (P433) ----

#' Mutant cell parameters
mut.div <- wt.div
mut.div[1] <- 1.8 * wt.div[1]
mut.div[2] <- 1.8 * wt.div[2]
mut.div[3] <- 1.8 * wt.div[3]
mut.div[4] <- 1.8 * wt.div[4]
mut.div[5] <- 1.8 * wt.div[5]
mut.div[7] <- 2.0 * wt.div[7]
mut.div[9] <- 1.8 * wt.div[9]

mut.sr <- wt.sr
mut.sr[4] <- 0.49
mut.sr[5] <- 0.49
mut.sr[7] <- 0.46

mut.adj <- wt.adj
mut.adj[2,3] <- 0.7
mut.adj[2,9] <- 1 - mut.adj[2,3]
mut.adj[3,4] <- 0.95
mut.adj[3,7] <- 1 - mut.adj[3,4]

mut.death <- wt.death
mut.death[5] <- 0.1
mut.death[7] <- 0.1
#mut.death[8] <- 0

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
    colour.vec = colour.vec,
    #alpha.vec = ifelse(seq_len(n) %in% c(2,5,8), 0.6, 1),
    show.tree = T, y.pos = y.pos,
    compartment.names = name.vec, lab.offset = -0.35, y.lim = -1.75
)
#savePlot(filename = "S6-A.pdf", width = 200, height.mult = 0.7)



# Red cell only (P328) ----

#' Mutant cell parameters
mut.div <- wt.div
mut.div[1] <- 1.8 * wt.div[1]
mut.div[2] <- 1.8 * wt.div[2]
mut.div[3] <- 1.8 * wt.div[3]
mut.div[4] <- 1.8 * wt.div[4]
mut.div[5] <- 1.8 * wt.div[5]
mut.div[7] <- 2.0 * wt.div[7]
mut.div[9] <- 1.8 * wt.div[9]

mut.sr <- wt.sr
mut.sr[5] <- 0.49
mut.sr[7] <- 0.46

mut.adj <- wt.adj
mut.adj[2,3] <- 0.7
mut.adj[2,9] <- 1 - mut.adj[2,3]
mut.adj[3,4] <- 0.4
mut.adj[3,7] <- 1 - mut.adj[3,4]

mut.death <- wt.death
mut.death[5] <- 0.1
mut.death[7] <- 0.1
#mut.death[8] <- 0

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
    colour.vec = colour.vec,
    #alpha.vec = ifelse(seq_len(n) %in% c(2,5,8), 0.6, 1),
    show.tree = T, y.pos = y.pos,
    compartment.names = name.vec, lab.offset = -0.35, y.lim = -1.75
)
#savePlot(filename = "S6-B.pdf", width = 200, height.mult = 0.7)



# Gradual platelet and late erythroid (P393) ----
 
#' Mutant cell parameters
mut.div <- wt.div
mut.div[1] <- 1.8 * wt.div[1]
mut.div[2] <- 1.8 * wt.div[2]
mut.div[3] <- 1.8 * wt.div[3]
mut.div[4] <- 1.8 * wt.div[4]
mut.div[5] <- 1.8 * wt.div[5]
mut.div[7] <- 2.0 * wt.div[7]
mut.div[9] <- 1.8 * wt.div[9]

mut.sr <- wt.sr
mut.sr[4] <- 0.49
mut.sr[5] <- 0.49
mut.sr[7] <- 0.46

mut.adj <- wt.adj

#mut.adj[2,3] <- 1.3 * wt.adj[2,3]
#mut.adj[2,9] <- 1 - mut.adj[2,3]
mut.adj[3,4] <- 1.4 * wt.adj[3,4]
mut.adj[3,7] <- 1 - mut.adj[3,4]

mut.death <- wt.death
mut.death[4] <- 0.1
mut.death[5] <- 0.1
mut.death[7] <- 0.1
#mut.death[8] <- 0

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
    colour.vec = colour.vec,
    #alpha.vec = ifelse(seq_len(n) %in% c(2,5,8), 0.6, 1),
    show.tree = T, y.pos = y.pos,
    compartment.names = name.vec, lab.offset = -0.35, y.lim = -1.75
)
#savePlot(filename = "S6-C.pdf", width = 200, height.mult = 0.7)
