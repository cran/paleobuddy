## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(dpi=300)

## ----results = "hide", message = FALSE----------------------------------------
# importing the package functions
library(paleobuddy)

## -----------------------------------------------------------------------------
# we set a seed so the results are reproducible
set.seed(1)

# set the necessary parameters
# initial number of species
n0 <- 1

# speciation rate - approx. 1 speciation event every 4my
# we are trying to create a big phylogeny so phytools can function better
lambda <- 0.25

# extinction rate - approx. 1 extinction event every 10my
mu <- 0.15

# maximum simulation time - species that die after this are considered extant
tMax <- 50

# run the simulation
sim <- bd.sim(n0, lambda, mu, tMax)

# take a look at the way the result is organized
sim

## -----------------------------------------------------------------------------
# draw simulation
draw.sim(sim, showLabel = FALSE)

## -----------------------------------------------------------------------------
# there are currently not many customization options for phylogenies
phy <- make.phylo(sim)

# plot it with APE - hide tip labels since there are a lot so it looks cluttered
ape::plot.phylo(phy, show.tip.label = FALSE)
ape::axisPhylo()

# plot the molecular phylogeny
ape::plot.phylo(ape::drop.fossil(phy), show.tip.label = FALSE)
ape::axisPhylo()

## -----------------------------------------------------------------------------
# set a seed
set.seed(3)

# create simulation
# note nExtant, defining we want 200 or more extant species at the end
sim <- bd.sim(n0, lambda, mu, tMax, nExtant = c(200, Inf))

# check the number of extant species
paste0("Number of species alive at the end of the simulation: ", 
       sum(sim$EXTANT))

## -----------------------------------------------------------------------------
# might look a bit cluttered
ape::plot.phylo(ape::drop.fossil(make.phylo(sim)), show.tip.label = FALSE)
ape::axisPhylo()

## -----------------------------------------------------------------------------
# we set a seed so the results are reproducible
set.seed(5)

# set the necessary parameters
# initial number of species
n0 <- 1

# speciation rate - it can be any function of time!
lambda <- function(t) {
  0.1 + 0.001*t
}

# extinction rate - also can be any function of time
mu <- function(t) {
  0.03 * exp(-0.01*t)
}

# maximum simulation time - species that die after this are considered extant
tMax <- 50

# run the simulation
sim <- bd.sim(n0, lambda, mu, tMax)

# check the resulting clade out
ape::plot.phylo(make.phylo(sim), show.tip.label = FALSE)
ape::axisPhylo()

## -----------------------------------------------------------------------------
# again set a seed
set.seed(1)

# set the sampling rate
# using a simple case - there will be on average T occurrences,
# per species, where T is the species duration
rho <- 1

# bins - used to represent the uncertainty in fossil occurrence times
bins <- seq(tMax, 0, -1)
# this is a simple 1my bin vector, but one could use the GSA timescale
# or something random, etc

# run the sampling simulation only for the first 10 species for brevity's sake
# returnAll = TRUE makes it so the occurrences are returned as binned as well
# (e.g. an occurrence at time 42.34 is returned as between 42 and 41)
fossils <- suppressMessages(sample.clade(sim = sim, rho = rho, 
                                      tMax = tMax, S = 1:10, 
                                      bins = bins, returnAll = TRUE))
# suppressing messages - the message is to inform the user how
# many speciesleft no fossils. In this case, it was 0

# take a look at how the output is organized
head(fossils)

## -----------------------------------------------------------------------------
# take only first 5 species with head
simHead <- head(sim, 10)

# draw longevities with fossil time points
suppressMessages(draw.sim(simHead, fossils = fossils))

## -----------------------------------------------------------------------------
# draw longevities with fossil time ranges
suppressMessages(draw.sim(simHead, fossils = fossils[, -3]))

## -----------------------------------------------------------------------------
# make a copy
pFossils <- fossils

# change the extant column
pFossils["Extant"][pFossils["Extant"] == FALSE] = "extant"
pFossils["Extant"][pFossils["Extant"] == TRUE] = "extinct"

# change column names
colnames(pFossils) <- c("Species", "Status", "min_age", "max_age")

# check it out
head(pFossils)

## -----------------------------------------------------------------------------
per.capita <- function(faBins, laBins, bins) {
  # create vectors to hold species that were born before and die in interval i,
  # species who were born in i and die later,
  # and species who were born before and die later
  NbL <- NFt <- Nbt <- rep(0, length(bins))
  
  # for each interval
  for (i in 1:length(bins)) {
    # number of species that were already around before i and are not seen again
    NbL[i] <- sum(faBins > bins[i] & laBins == bins[i])
    
    # number of species that were first seen in i and are seen later
    NFt[i] <- sum(faBins == bins[i] & laBins < bins[i])
    
    # number of species that were first seen before i and are seen after i
    Nbt[i] <- sum(faBins > bins[i] & laBins < bins[i])
  }
  
  # calculate the total rates
  p <- log((NFt + Nbt) / Nbt)
  q <- log((NbL + Nbt) / Nbt)
  return(list(p = p, q = q))
}

## -----------------------------------------------------------------------------
# get the species names
ids <- unique(fossils$Species)

# get the first appearance bins - the first time in bins where the fossil was seen (lower bound)
faBins <- unlist(lapply(ids, function(i) max(fossils$MaxT[fossils$Species == i])))

# get the last appearance bins - last time in bins where the fossil was seen (upper bound)
laBins <- unlist(lapply(ids, function(i) min(fossils$MinT[fossils$Species == i])))

# create the bins vector we have been using
bins <- seq(tMax, 0, -0.1)
# note this has a high resolution, the actual stratigraphic ranges are much coarser

# get the estimates
pc <- per.capita(faBins, laBins, bins)

## -----------------------------------------------------------------------------
# set a seed
set.seed(2)

# parameters to set things up
n0 <- 1
tMax <- 20

# speciation can be dependent on a time-series variable as well as time
lambda <- function(t, env) {
  0.01 * t + 0.01*exp(0.01*env)
}

# let us use the package's temperature data
data(temp)
# this could instead  be data(co2), the other environmental
# data frame supplied by paleobuddy

# we can make extinction be age-dependent by creating a shape parameter
mu <- 10
mShape <- 0.5
# this will make it so species durations are distributed
# as a Weibull with scale 10 and shape 2

# run the simulation
# we pass the shape and environmental parameters
sim <- suppressMessages(bd.sim(n0, lambda, mu, tMax, 
                               mShape = mShape, envL = temp))
# note that lShape and envM also exist
# the defaults for all of these customization options is NULL

# check out the phylogeny
ape::plot.phylo(make.phylo(sim), show.tip.label = FALSE)
ape::axisPhylo()

## -----------------------------------------------------------------------------
# speciation may be a step function, presented
# as a rates vector and a shift times vector
lList <- c(0.1, 0.2, 0.05)
lShifts <- c(0, 10, 15)
# lShifts could also be c(tMax, tMax - 10, tMax - 15) for identical results

# make.rate is the function bd.sim calls to create a 
# ratem but we will use it here to see how it looks
lambda <- make.rate(lList, tMax = 20, rateShifts = lShifts)
t <- seq(0, 20, 0.1)
plot(t, rev(lambda(t)), type = 'l', main = "Step function speciation rate",
     xlab = "Time (Mya)", ylab = "Speciation rate", xlim = c(20, 0))

## -----------------------------------------------------------------------------
# it is not possible to combine the lambda(t, env) and lList methods listed above
# but we can create a step function dependent on environmental data with ifelse
mu_t <- function(t, env) {
  ifelse(t < 10, env,
         ifelse(t < 15, env * 2, env / 2))
}

# pass it to make.rate with environmental data
mu <- make.rate(mu_t, tMax = tMax, envRate = temp)
plot(t, rev(mu(t)), type = 'l', main = "Environmental step function extinction rate",
     xlab = "Time (Mya)", ylab = "Rate", xlim = c(20, 0))

## ----eval=FALSE---------------------------------------------------------------
#  # set seed again
#  set.seed(1)
#  
#  # age-dependent speciation with a step function
#  lList <- c(10, 5, 12)
#  lShifts <- c(0, 10, 15)
#  lShape <- 2
#  
#  # age-dependent extinction with a step function of an environmental variable
#  q <- function(t, env) {
#    ifelse(t < 10, 2*env,
#           ifelse(t < 15, 1.4*env, env / 2))
#  }
#  
#  # note shape can be time-dependent as well, though
#  # we advise for variation not to be too abrupt due
#  # to computational issues
#  mShape <- function(t) {
#    return(1.2 + 0.025*t)
#  }
#  
#  # run the simulation
#  sim <- suppressMessages(bd.sim(n0, lList, q, tMax, lShape = lShape,
#                                 mShape = mShape,
#                                 envM = temp, lShifts = lShifts))
#  
#  # check out the phylogeny
#  ape::plot.phylo(make.phylo(sim), show.tip.label = FALSE)

## -----------------------------------------------------------------------------
# as an example, we will use a PERT distribution, 
# a hat-shaped distribution used in PyRate

# preservation function
dPERT <- function(t, s, e, sp, a = 3, b = 3, log = FALSE) {

  # check if it is a valid PERT
  if (e >= s) {
    message("There is no PERT with e >= s")
    return(rep(NaN, times = length(t)))
  }

  # find the valid and invalid times
  id1 <- which(t <= e | t >= s)
  id2 <- which(!(t <= e | t >= s))
  t <- t[id2]

  # initialize result vector
  res <- vector()

  # if user wants a log function
  if (log) {
    # invalid times get -Inf
    res[id1] <- -Inf

    # valid times calculated with log
    res[id2] <- log(((s - t) ^ 2)*((-e + t) ^ 2)/((s - e) ^ 5*beta(a,b)))
  }
  
  # otherwise
  else{
    res[id1] <- 0

    res[id2] <- ((s - t) ^ 2)*((-e + t) ^ 2)/((s - e) ^ 5*beta(a,b))
  }

  return(res)
}

# set seed
set.seed(1)

# generate a quick simulation
sim <- bd.sim(n0, lambda = 0.1, mu = 0.05, tMax = 20)

# sample for the first 10 species
fossils <- suppressMessages(sample.clade(sim = sim, rho = 3, 
                                      tMax = 20, adFun = dPERT))
# here we return true times of fossil occurrences

# draw longevities with fossil occurrences
draw.sim(sim, fossils = fossils)

## -----------------------------------------------------------------------------
# set a seed 
set.seed(1)

# simple simulation, starting with more than one species
sim <- bd.sim(n0 = 2, lambda = 0.1, mu = 0, tMax = 20, nFinal = c(20, Inf))

# separate the lineages
clades <- find.lineages(sim)

# plot each phylogeny

# clade 1
ape::plot.phylo(make.phylo(clades$clade_1$sim), show.tip.label = FALSE)

# clade 2
ape::plot.phylo(make.phylo(clades$clade_2$sim), show.tip.label = FALSE)

