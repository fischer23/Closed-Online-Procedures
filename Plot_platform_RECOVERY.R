# Used to calculate Table 2 in the paper "The Online Closure Principle"

rm(list = ls())
library(ggplot2)
library(MASS)
library(patchwork)
library(Matrix)
library(mvtnorm)

### load the procedures
source("OCP_Procedures.R")

# Test is not finished for Sotrovimab. It is only included to determine \alpha_{13}.

names <- c(
  "Dexamethasone", "Lopinavir-ritonavir", "Hydroxychloroquine", "Azithromycin",
  "Tocilizumab", "Convalescent plasma", "Casirivimab-Imdevimab", "Aspirin",
  "Colchicine", "Baricitinib", "High-dose steroids", "Empagliflozin", "Sotrovimab"
)

lags <- c(0, 1, 2, 3, 4, 5, 3, 3, 3, 3, 1, 2, 2)

p <- c(0.0003, 0.58, 0.1, 0.99, 0.007, 0.34, 0.001, 0.35, 0.63, 0.026, 0.0012, 0.64, 1)

n <- 13

### We run the experiment for different gamma
n_q <- 3

alpha_ind_Spending <- matrix(0, ncol = n_q, nrow = n)
rej_Spending <- matrix(0, ncol = n_q, nrow = n)
n_rej_Spending <- rep(0, n_q)

alpha_ind_closed_ADDIS <- matrix(0, ncol = n_q, nrow = n)
rej_closed_ADDIS <- matrix(0, ncol = n_q, nrow = n)
n_rej_closed_ADDIS <- rep(0, n_q)

alpha_ind_Graph_online <- matrix(0, ncol = n_q, nrow = n)
rej_Graph_online <- matrix(0, ncol = n_q, nrow = n)
n_rej_Graph_online <- rep(0, n_q)

alpha_ind_closed_Alpha <- matrix(0, ncol = n_q, nrow = n)
rej_closed_Alpha <- matrix(0, ncol = n_q, nrow = n)
n_rej_closed_Alpha <- rep(0, n_q)

alpha_ind_Alpha <- matrix(0, ncol = n_q, nrow = n)
rej_Alpha <- matrix(0, ncol = n_q, nrow = n)
n_rej_Alpha <- rep(0, n_q)



count <- 1

for (q in c(0.6, 0.7, 0.8)) {
  ### Initialise Hyperparameters
  gamma <- q^(1:n) * (1 - q) / q
  tau <- 0.8
  lambda <- 0.16
  w <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  w[w == 0] <- 1
  w <- matrix(gamma[w], n, n)
  w[upper.tri(w) == 0] <- 0
  alpha <- 0.05


  ## ADDIS-Spending
  alpha_ind_Spending[, count] <- ADDIS_Spending(alpha, gamma, tau, lambda, lags, p, n)
  rej_Spending[, count] <- alpha_ind_Spending[, count] >= p
  n_rej_Spending[count] <- sum(rej_Spending[, count])

  ## Closed ADDIS-Spending
  alpha_ind_closed_ADDIS[, count] <- closed_ADDIS_Spending(alpha, gamma, tau, lambda, lags, p, n)
  rej_closed_ADDIS[, count] <- alpha_ind_closed_ADDIS[, count] >= p
  n_rej_closed_ADDIS[count] <- sum(rej_closed_ADDIS[, count])

  ## Online-Graph
  alpha_ind_Graph_online[, count] <- Online_Graph(alpha, gamma, w, p, n)
  rej_Graph_online[, count] <- alpha_ind_Graph_online[, count] >= p
  n_rej_Graph_online[count] <- sum(rej_Graph_online[, count])

  # Closed Alpha-Spending
  alpha_ind_closed_Alpha[, count] <- closed_Alpha_Spending(alpha, gamma, p, n)
  rej_closed_Alpha[, count] <- alpha_ind_closed_Alpha[, count] >= p
  n_rej_closed_Alpha[count] <- sum(rej_closed_Alpha[, count])

  # Alpha-Spending
  alpha_ind_Alpha[, count] <- alpha_Spending(alpha, gamma, p, n)
  rej_Alpha[, count] <- alpha_ind_Alpha[, count] >= p
  n_rej_Alpha[count] <- sum(rej_Alpha[, count])

  count <- count + 1
}
