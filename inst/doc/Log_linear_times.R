## ----include = FALSE----------------------------------------------------------
set.seed(20241017)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
K <- 10^5

## ----setup--------------------------------------------------------------------
library(data.table)
library(nhppp)

## ----population---------------------------------------------------------------
pop <- setDT(
  list(
    id = 1:K,
    alpha = stats::rnorm(n = K, mean = -4, sd = 0.5),
    beta = truncnorm::rtruncnorm(n = K, mean = 0.03, sd = 0.003, a = 0, b = Inf),
    T0 = rep(40, K),
    T1 = stats::runif(n = K, min = 50, max = 100)
  )
)
setindex(pop, id)
pop

## ----function-lambda-0--------------------------------------------------------
l <- function(t, alpha = pop$alpha, beta = pop$beta, ...) exp(alpha + beta * t)

## -----------------------------------------------------------------------------
l(t = 45, alpha = -4, beta = 0.03)

## -----------------------------------------------------------------------------
l(t = 45:50, alpha = -4, beta = 0.03)

l(t = matrix(45:50, ncol = 1), alpha = -4, beta = 0.03)

## -----------------------------------------------------------------------------
l(t = 45, alpha = pop$alpha[1:5], beta = pop$beta[1:5])

## -----------------------------------------------------------------------------
l(t = matrix(c(45, 50, 45, 47.4, 30), ncol = 1), alpha = pop$alpha[1:5], beta = pop$beta[1:5])

## -----------------------------------------------------------------------------
t_mat <- matrix(c(
  45, 50, 45, 47.4, 30,
  45.1, 50.1, 45.5, 47.8, 38,
  48, 52.7, 60.1, 70.1, 99.9
), ncol = 3, byrow = TRUE)
t_mat

l(t = t_mat, alpha = pop$alpha[1:5], beta = pop$beta[1:5])

## -----------------------------------------------------------------------------
t_mat <- matrix(rep(c(45, 50, 55), each = K), ncol = 3, byrow = TRUE)
l(t = t_mat) |> head()

## ----Lambda-------------------------------------------------------------------
L <- function(t, alpha = pop$alpha, beta = pop$beta, ...) {
  (exp(alpha + beta * t) - exp(alpha)) / beta
}

## ----Lambda-inv---------------------------------------------------------------
Li <- function(z, alpha = pop$alpha, beta = pop$beta, ...) {
  (log(beta * z + exp(alpha)) - alpha) / beta
}

## ----method-1-----------------------------------------------------------------
tictoc::tic("Method 1 (nonvectorized)")
t_nonvec_special_case <- rep(NA, K)
for (k in 1:K) {
  t1 <- nhppp::draw_sc_loglinear(
    intercept = pop$alpha[k],
    slope = pop$beta[k],
    t_min = pop$T0[k],
    t_max = pop$T1[k],
    atmost1 = TRUE
  )
  if (length(t1) != 0) {
    t_nonvec_special_case[k] <- t1
  }
}
tictoc::toc(log = TRUE)
pop[, t_nonvec_special_case := t_nonvec_special_case]

## ----method-2-----------------------------------------------------------------
tictoc::tic("Method 2 (vectorized, thinning)")
M <- 5
break_points <- seq.int(from = 40, to = 100, length.out = M + 1)
breaks_mat <- matrix(rep(break_points, each = K), nrow = K)

lmaj_mat <- nhppp::get_step_majorizer(
  fun = l,
  breaks = breaks_mat,
  is_monotone = TRUE
)

pop[
  ,
  t_thinning := nhppp::vdraw_intensity(
    lambda = l,
    lambda_maj_matrix = lmaj_mat,
    rate_matrix_t_min = 40,
    rate_matrix_t_max = 100,
    t_min = pop$T0,
    t_max = pop$T1,
    atmost1 = TRUE,
    atmostB = 5
  )
]
tictoc::toc(log = TRUE) # timer end

## ----method-3-----------------------------------------------------------------
tictoc::tic("Method 3 (inversion)")
pop[
  ,
  t_inversion := nhppp::vdraw_cumulative_intensity(
    Lambda = L,
    Lambda_inv = Li,
    t_min = pop$T0,
    t_max = pop$T1,
    atmost1 = TRUE
  )
]
tictoc::toc(log = TRUE) # timer end

## ----qq-plots, fig.alt="QQ plots comparing simulated times with the three methods. The QQ plots indicate excellent agreement."----
qqplot(pop$t_nonvec_special_case, pop$t_thinning)
qqplot(pop$t_nonvec_special_case, pop$t_inversion)

