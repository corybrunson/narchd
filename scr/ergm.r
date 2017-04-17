### exponential random graph models

## setup

rm(list = ls())
# working directory
setwd("~/Documents/CQM/network analysis of clinical data/")
# packages
pkgs <- c("igraph", "ergm")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    stopifnot(require(pkg, character.only = TRUE))
  }
}

# load and tidy data
source("scr/setup/tidy.r")

# read graph object
load("data/scientometrics/cite_graph.rdata")

# full citation graph
cite_ntwk <- as.network(as_adj(cite_graph))
plot(cite_ntwk)

# predictor variables
## common authorship
article_authors <- incl_dat$authors %>% strsplit(split = "\\|")
common_authors <- sapply(1:length(article_authors), function(i) {
  sapply(1:length(article_authors), function(j) {
    length(intersect(article_authors[[i]], article_authors[[j]]))
  })
})
## manual grouping
network::set.vertex.attribute(cite_ntwk, "group", incl_dat$review_group)
## type of publication
network::set.vertex.attribute(cite_ntwk, "venue", V(cite_graph)$venue)
## delay between publication dates
n <- vcount(cite_graph)
debut_val <- as.numeric(as.Date(V(cite_graph)$debut))
delay <- matrix(rep(debut_val, times = n), ncol = n) -
  matrix(rep(debut_val, each = n), ncol = n)
### transform for symmetry and non-negativity, based on inspection
hist(delay[which(as.matrix(cite_ntwk) == 1)])
table(cut(delay[which(as.matrix(cite_ntwk) == 1)],
          breaks = seq(-7000, 7000, 1000))) /
  table(cut(delay,
            breaks = seq(-7000, 7000, 1000)))
table(cut(rnorm(n = 1000, mean = 1461, sd = 1461),
          breaks = seq(-7000, 7000, 1000)))
### caution: linear recency term makes "sleeping beauties" nigh impossible
delay_abs <- abs(delay - 1461)
delay_norm <- exp(-(delay - 1461)^2 / (2 * 1461^2))

# routine goodness-of-fit test
routine.gof <- function(mod) {
  res <- gof(mod,
             GOF = ~ idegree + odegree + triadcensus +
               espartners + dspartners + distance,
             verbose = TRUE, interval = 1000)
  par(mfrow = c(2, 3))
  plot(res, plotlogodds = TRUE)
  par(mfrow = c(1, 1))
  res
}

# models
## frequency
mod1 <- ergm(cite_ntwk ~ edges)
plot(simulate(mod1))
gof1 <- routine.gof(mod1)
summary(mod1)
## frequency + venue
mod2 <- ergm(cite_ntwk ~ edges +
               nodeifactor("venue", base = 1))
plot(simulate(mod2))
gof2 <- routine.gof(mod2)
summary(mod2)
## frequency + venue + common authorship
mod3a <- ergm(cite_ntwk ~ edges +
                nodeifactor("venue", base = 1) +
                edgecov(common_authors))
plot(simulate(mod3a))
gof3a <- routine.gof(mod3a)
summary(mod3a)
## frequency + venue + group
mod3b <- ergm(cite_ntwk ~ edges +
                nodeifactor("venue", base = 1) +
                nodematch(attrname = "group", diff = TRUE))
plot(simulate(mod3b))
gof3b <- routine.gof(mod3b)
summary(mod3b)
## frequency + venue + recency
theta3 <- c(mod2$coef, edgecov.delay_norm = -1)
mod3 <- ergm(cite_ntwk ~ edges +
               nodeifactor("venue", base = 1) +
               #edgecov(x = delay),
               #edgecov(x = delay_abs),
               edgecov(x = delay_norm),
             verbose = TRUE, control = control.ergm(init = theta3))
plot(simulate(mod3))
gof3 <- routine.gof(mod3)
summary(mod3)
## frequency + venue + recency + impacts
theta4 <- c(mod3$coef, gwidegree = -1)
mod4 <- ergm(cite_ntwk ~ edges +
               nodeifactor("venue", base = 1) +
               edgecov(x = delay_norm) +
               gwidegree(decay = 0, fixed = FALSE, cutoff = 7),
             verbose = TRUE, control = control.ergm(init = theta4))
plot(simulate(mod4))
gof4 <- routine.gof(mod4)
summary(mod4)
## frequency + venue + recency + impacts + shared partners
### dyadwise
theta5a <- c(enformulate.curved(mod4)$theta, gwdsp = 1)
mod5a <- ergm(cite_ntwk ~ edges +
               nodeifactor("venue", base = 1) +
               edgecov(x = delay_norm) +
               gwidegree(decay = 0, fixed = FALSE, cutoff = 7) +
               gwdsp(alpha = 0, fixed = FALSE, cutoff = 3),
             verbose = TRUE,
             control = control.ergm(init = theta5a))
plot(simulate(mod5a))
gof5a <- routine.gof(mod5a)
summary(mod5a)
### edgewise
theta5b <- c(enformulate.curved(mod4)$theta, gwesp = 1)
mod5b <- ergm(cite_ntwk ~ edges +
               nodeifactor("venue", base = 1) +
               edgecov(x = delay_norm) +
               gwidegree(decay = 0, fixed = FALSE, cutoff = 7) +
               gwesp(alpha = 0, fixed = FALSE, cutoff = 3),
             verbose = TRUE,
             control = control.ergm(init = theta5b))
plot(simulate(mod5b))
gof5b <- routine.gof(mod5b)
summary(mod5b)
### both dyadwise and edgewise
theta5c <- c(enformulate.curved(mod4)$theta, gwdsp = 1, gwesp = 1)
mod5c <- ergm(cite_ntwk ~ edges +
               nodeifactor("venue", base = 1) +
               edgecov(x = delay_norm) +
               gwidegree(decay = 0, fixed = FALSE, cutoff = 7) +
               gwdsp(alpha = 0, fixed = FALSE, cutoff = 3) +
               gwesp(alpha = 0, fixed = FALSE, cutoff = 3),
             verbose = TRUE,
             control = control.ergm(init = theta5c))
plot(simulate(mod5c))
gof5c <- routine.gof(mod5c)
summary(mod5c)
### best fit
mod5 <- mod5b
gof5 <- gof5b
## frequency + common authors + venue + recency + impacts + shared partners
mod6a <- ergm(cite_ntwk ~ edges +
                edgecov(x = common_authors) +
                nodeifactor("venue", base = 1) +
                edgecov(x = delay_norm) +
                gwidegree(decay = 0, fixed = FALSE, cutoff = 7) +
                #gwdsp(alpha = 0, fixed = FALSE, cutoff = 3) +
                gwesp(alpha = 0, fixed = FALSE, cutoff = 3),
              verbose = TRUE)
plot(simulate(mod6a))
gof6a <- routine.gof(mod6a)
summary(mod6a)
## freq + groups + venue + recency + impacts + shared partners
mod6b <- ergm(cite_ntwk ~ edges +
                nodematch(attrname = "group", diff = TRUE) +
                nodeifactor("venue", base = 1) +
                edgecov(x = delay_norm) +
                gwidegree(decay = 0, fixed = FALSE, cutoff = 7) +
                #gwdsp(alpha = 0, fixed = FALSE, cutoff = 3) +
                gwesp(alpha = 0, fixed = FALSE, cutoff = 3),
              verbose = TRUE)
plot(simulate(mod6b))
gof6b <- routine.gof(mod6b)
summary(mod6b)

# restrict to journal articles
cite_graph <- delete_vertices(cite_graph,
                              which(V(cite_graph)$venue != "Journal"))
cite_ntwk <- as.network(as_adj(cite_graph))
plot(cite_ntwk)

## frequency
mod1 <- ergm(cite_ntwk ~ edges)
summary(mod1)
plot(simulate(mod1))
gof1 <- gof(mod1,
            GOF = ~ idegree + odegree + triadcensus +
              espartners + dspartners + distance,
            verbose = TRUE, interval = 1000)
par(mfrow = c(2, 3))
plot(gof1, plotlogodds = TRUE)
par(mfrow = c(1, 1))
if (FALSE) {
  ## frequency + recency
  network::set.vertex.attribute(cite_ntwk, "debut",
                                as.numeric(as.Date(V(cite_graph)$debut)))
  el <- as_edgelist(cite_graph, names = FALSE)
  ### publication dates
  debut_val <- as.numeric(as.Date(V(cite_graph)$debut))
  hist(debut_val)
  hist(max(debut_val) + 365 - debut_val)
  hist(log(max(debut_val) + 365 - debut_val))
  now <- max(debut_val) + 365
  ### differences in publication dates
  debut_diff <- as.numeric(as.Date(V(cite_graph)$debut[el[, 1]])) -
    as.numeric(as.Date(V(cite_graph)$debut[el[, 2]]))
  hist(debut_diff)
  hist(log(365 + debut_diff))
  hist(debut_diff ^ .5)
  ### models
  theta2a <- c(mod1$coef, nodecov.debut = .1)
  mod2a <- ergm(cite_ntwk ~ edges +
                  nodecov("debut", transform = function(x) log(now - x)),
                verbose = TRUE,
                control = control.ergm(init = theta2a))
  summary(mod2a)
  plot(simulate(mod2a))
  gof2a <- gof(mod2a,
               GOF = ~ idegree + odegree + triadcensus +
                 espartners + dspartners + distance,
               verbose = TRUE, interval = 1000)
  par(mfrow = c(2, 3))
  plot(gof2a, plotlogodds = TRUE)
  par(mfrow = c(1, 1))
  theta2b <- c(mod1$coef, absdiff0.5.debut = .1)
  mod2b <- ergm(cite_ntwk ~ edges +
                  absdiff("debut", pow = .5),
                verbose = TRUE,
                control = control.ergm(init = theta2b))
  summary(mod2b)
  plot(simulate(mod2b))
  gof2b <- gof(mod2b,
               GOF = ~ idegree + odegree + triadcensus +
                 espartners + dspartners + distance,
               verbose = TRUE, interval = 1000)
  par(mfrow = c(2, 3))
  plot(gof2b, plotlogodds = TRUE)
  par(mfrow = c(1, 1))
  theta2c <- c(mod1$coef, nodecov.debut = .1, absdiff0.5.debut = .1)
  mod2c <- ergm(cite_ntwk ~ edges +
                  nodecov("debut", transform = function(x) log(now - x)) +
                  absdiff("debut", pow = .5),
                verbose = TRUE,
                control = control.ergm(init = theta2c))
  summary(mod2c)
  plot(simulate(mod2c))
  gof2c <- gof(mod2c,
               GOF = ~ idegree + odegree + triadcensus +
                 espartners + dspartners + distance,
               verbose = TRUE, interval = 1000)
  par(mfrow = c(2, 3))
  plot(gof2c, plotlogodds = TRUE)
  par(mfrow = c(1, 1))
  mod2 <- mod2b
}
## frequency + recency
n <- vcount(cite_graph)
debut_val <- as.numeric(as.Date(V(cite_graph)$debut))
### matrix of differences
delay <- matrix(rep(debut_val, times = n), ncol = n) -
  matrix(rep(debut_val, each = n), ncol = n)
### models
theta2 <- c(mod1$coef, edgecov.delay = 0)
mod2 <- ergm(cite_ntwk ~ edges + edgecov(x = delay),
             verbose = TRUE, control = control.ergm(init = theta2))
summary(mod2)
plot(simulate(mod2))
gof2 <- gof(mod2,
            GOF = ~ idegree + odegree + triadcensus +
              espartners + dspartners + distance,
            verbose = TRUE, interval = 1000)
par(mfrow = c(2, 3))
plot(gof2, plotlogodds = TRUE)
par(mfrow = c(1, 1))
## frequency + recency + sender + receiver (NEED TO USE LARGEST COMPONENT)
mod3 <- ergm(cite_ntwk ~ edges + edgecov(x = delay) +
               sender + receiver,
             verbose = TRUE)
summary(mod3)
plot(simulate(mod3))
gof3 <- gof(mod3,
            GOF = ~ idegree + odegree + triadcensus +
              espartners + dspartners + distance,
            verbose = TRUE, interval = 1000)
par(mfrow = c(2, 3))
plot(gof3, plotlogodds = TRUE)
par(mfrow = c(1, 1))
## frequency + recency + transitivity
theta4a <- c(mod3$coef, transitiveties = .1)
mod4a <- ergm(cite_ntwk ~ edges +
                absdiff("debut", pow = .5) +
                sender + receiver +
                transitiveties,
              verbose = TRUE,
              control = control.ergm(init = theta4a, MCMLE.density.guard = 50))
summary(mod4a)
theta4b <- c(mod2$coef, colSums(sna::triad.census(cite_ntwk))[-1])
mod4b <- ergm(cite_ntwk ~ edges + absdiff("debut", pow = .5) + triadcensus,
              verbose = TRUE,
              control = control.ergm(init = theta4b, MCMLE.density.guard = 50))
summary(mod4b)


cite_subgraph <- delete_vertices(cite_graph, which(degree(cite_graph) == 0))
cite_comp <- components(cite_subgraph)$membership
wh_lc <- which(cite_comp == names(which.max(table(cite_comp))))
cite_lc <- induced_subgraph(cite_subgraph, wh_lc)

cite_lc_ntwk <- as.network(as_adj(cite_lc))
plot(cite_lc_ntwk)
