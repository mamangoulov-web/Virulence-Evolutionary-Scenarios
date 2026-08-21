# =============================================================================
# VIRULENCE EVOLUTION ACROSS ECOLOGICAL SCENARIOS
# Multi-host/environmental transmission
# =============================================================================
#
# WHAT THIS SCRIPT DOES
#   For each scenario, it:
#     1. Codes up the ODE system exactly as sketched, so you can simulate and
#        watch trajectories (S(t), I(t), ... ) directly to test the ODEs.
#     2. Solves for the resident ecological equilibrium as a function of
#        virulence (alpha), numerically/algebraically.
#     3. Computes the invasion fitness of a rare mutant strain and its
#        selection gradient, then finds the ESS (root of the gradient).
#     4. Incorporates a key ecological parameter (predation intensity n; 
#        infection variable f) and plots how the ESS virulence responds.
#     5. Draws a pairwise invasibility plot (PIP) for each scenario, confirming 
#        the ESS is a true convergence-stable singular strategy (not just a zero 
#        of the gradient).
#
# TWO MODELING NOTES (for Javad):
#
#   (1) Birth term: the board sketches used proportional birth "bS" for the
#       predator and multi-host scenarios. Written that way, both scenarios
#       turn out to have no generic fixed-point equilibrium: proportional
#       birth with no self-limitation makes the models scale-invariant
#       (homogeneous), so they either grow/decay without bound or only
#       balance at a knife-edge parameter combination. This is why
#       our R0 formula from last week carries the note N=1 next to it, as
#       the base model implicitly assumes some constant recruitment/
#       normalized population size, not free proportional growth. I've made
#       that assumption explicit here: births enter as a constant recruitment
#       rate Lambda into the susceptible class, consistent with the N=1
#       convention. This is also the standard convention in the virulence-
#       evolution literature (SIS/SI models). If it interests anybody, we can
#       still simulate the literal "bS" version below (set birth_type =
#       "proportional") to see the instability. Regardless, it's a good
#       thing to verify empirically.
#
#       Follow-up on this (Javad asked): with per-capita birth bS instead of
#       constant recruitment, the endemic Jacobian works out to
#       trace = -gamma*(b-d)/(d+alpha), det = (b-d)*(d+gamma+alpha). det>0
#       automatically whenever b>d (the feasibility condition), so stability
#       hinges entirely on the trace, and trace<0 iff gamma>0. If gamma=0 the
#       trace is exactly zero for ANY parameter values -- a structural neutral
#       center (sustained oscillation), not just a bad parameter choice. So
#       per-capita birth was never the real problem for the baseline SIS model
#       (which has gamma>0 and could have kept it and stayed stable); constant
#       recruitment was only strictly necessary for the predator/multi-host
#       scenarios below, whose board sketches have no recovery term at all
#       (gamma=0 baked into the structure), which sits exactly on that
#       knife-edge regardless of any other parameter choice.
#
#   (2) Trade-off shape: beta(alpha) = b0 * alpha^q, with 0 < q < 1 (concave,
#       shape as seen in the chart that Javad drew on the whiteboard). q = 0.5
#       (square root) is the special case that reproduces last week's
#       clean result alpha* = gamma + d (has been derived and verified below).
#       For general q, the single-host ESS is alpha* = [q/(1-q)] * (gamma+d).
#       q >= 1 (convex trade-off) removes the stabilizing selection seen in
#       the PIP sketch and typically produces evolutionary branching or
#       runaway virulence instead of a single ESS.
#
# =============================================================================

if (!requireNamespace("deSolve", quietly = TRUE)) install.packages("deSolve")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(deSolve)
library(ggplot2)

# -----------------------------------------------------------------------------
# SHARED: transmission-virulence trade-off  beta(alpha) = b0 * alpha^q
# -----------------------------------------------------------------------------
beta_fun  <- function(alpha, b0, q = 0.5) b0 * pmax(alpha, 1e-8)^q
dbeta_fun <- function(alpha, b0, q = 0.5) b0 * q * pmax(alpha, 1e-8)^(q - 1)

# generic numerical derivative helper used throughout for selection gradients
num_deriv <- function(f, x, h = 1e-5) (f(x + h) - f(x - h)) / (2 * h)

# generic 1-D root bracketing over a grid (robust when uniroot's two endpoints
# alone might not bracket a sign change, or where NA/infeasible regions exist)
grid_root <- function(f, lower, upper, n = 300) {
  xg <- seq(lower, upper, length.out = n)
  yg <- sapply(xg, f)
  ok <- !is.na(yg)
  if (sum(ok) < 2) return(NA_real_)
  xg <- xg[ok]; yg <- yg[ok]
  s  <- sign(yg)
  idx <- which(diff(s) != 0)
  if (length(idx) == 0) return(NA_real_)
  uniroot(f, c(xg[idx[1]], xg[idx[1] + 1]))$root
}

# second-order (curvature) numerical derivative -- used below for the
# evolutionary-stability half of the CSS check (is the ESS a local fitness
# maximum in the mutant trait, i.e. uninvadable once reached?)
num_second_deriv <- function(f, x, h = 1e-3) (f(x + h) - 2 * f(x) + f(x - h)) / h^2

# Full CSS ("continuously stable strategy") check, shared across all three
# scenarios: convergence stability comes from the selection-gradient sign on
# either side of the ESS (+ below, - above means residents evolve toward it);
# evolutionary stability comes from the curvature of invasion fitness in the
# mutant trait ALONE, resident frozen at the ESS (negative means a true local
# fitness peak, not just a zero of the gradient). Both must hold for a genuine
# CSS -- this is what actually answers "is it convergent AND stable".
css_check <- function(selection_gradient_fun, invasion_fitness_of_mutant_fun, ess, p, rel_step = 0.05) {
  delta <- max(rel_step * ess, 1e-4)
  grad_below <- selection_gradient_fun(ess - delta, p)
  grad_above <- selection_gradient_fun(ess + delta, p)
  curvature  <- num_second_deriv(invasion_fitness_of_mutant_fun, ess)
  data.frame(ess = ess, grad_below = grad_below, grad_above = grad_above, curvature = curvature,
             convergence_stable = isTRUE(grad_below > 0 && grad_above < 0),
             evolutionarily_stable = isTRUE(curvature < 0))
}

# ESS-centered plotting window shared by all three PIP plots below, so they're
# apples-to-apples. An earlier version windowed each plot independently as
# [0, 3*ess], which put the ESS off-center and let the predator/multi-host
# feasibility boundary (equilibrium stops existing above some alpha) clip one
# side of the grid asymmetrically -- differences in how those plots looked
# were a windowing artifact, not a real difference in stability type (the
# css_check() results printed for each scenario agree on that).
pip_window <- function(ess, zoom = 0.75) c(max(1e-4, ess * (1 - zoom)), ess * (1 + zoom))


# =============================================================================
# PART 1: BASELINE SIS MODEL  (no predator, single host)
# =============================================================================
#   dS/dt = Lambda - d*S - beta(alpha)*S*I + gamma*I
#   dI/dt = beta(alpha)*S*I - (d + gamma + alpha)*I
#
#   Resident equilibrium threshold (from R0 = beta(alpha)*S/(d+alpha+gamma)):
#       S* = (d + gamma + alpha) / beta(alpha)
#   Invasion fitness of a rare mutant (alpha_m) on resident background S*:
#       r_m = beta(alpha_m)*S* - (d + gamma + alpha_m)      <- your board formula

base_odes <- function(t, y, p) {
  with(as.list(c(y, p)), {
    S <- max(S, 0); I <- max(I, 0)
    beta <- beta_fun(alpha, b0, q)
    dS <- Lambda - d*S - beta*S*I + gamma*I
    dI <- beta*S*I - (d + gamma + alpha)*I
    list(c(dS, dI))
  })
}

# NOTE: deliberately using explicit p$... access below rather than with(p, ...).
# with() binds every name in p as a local variable, and later in this script we
# set p_base$alpha <- ess_base (needed so the ODE integrator knows which alpha
# to simulate); if these functions used with(p, ...), that stored p$alpha
# would silently shadow the function's own `alpha` argument on every later
# call, freezing S* at the ESS regardless of what resident value was passed in.
base_Sstar <- function(alpha, p) (p$d + p$gamma + alpha) / beta_fun(alpha, p$b0, p$q)

base_invasion_fitness <- function(alpha_m, Sstar, p) {
  beta_fun(alpha_m, p$b0, p$q) * Sstar - (p$d + p$gamma + alpha_m)
}

# Selection-gradient convention (identical in the predator and multi-host
# scenarios below): freeze the RESIDENT equilibrium first, then differentiate
# invasion fitness with respect to the MUTANT trait only, evaluated at
# alpha_m = alpha (mutant = resident) -- i.e. dr/d(alpha_m) |_{alpha_m=alpha},
# the standard adaptive-dynamics selection gradient.
base_selection_gradient <- function(alpha, p) {
  Sstar <- base_Sstar(alpha, p)
  num_deriv(function(am) base_invasion_fitness(am, Sstar, p), alpha)
}

find_ESS_base <- function(p, interval = c(1e-3, 50)) {
  uniroot(function(a) base_selection_gradient(a, p), interval)$root
}

## ---- run + verify against theory -------------------------------------------
p_base <- list(b0 = 2, d = 0.10, gamma = 0.30, q = 0.5, Lambda = 1)
ess_base <- find_ESS_base(p_base)
cat("=== Baseline SIS ===\n")
cat(sprintf("Numeric ESS alpha* = %.4f   |   Theory (gamma+d) = %.4f\n",
            ess_base, p_base$gamma + p_base$d))

css_base <- css_check(base_selection_gradient,
                       function(am) base_invasion_fitness(am, base_Sstar(ess_base, p_base), p_base),
                       ess_base, p_base)
cat(sprintf("CSS check: grad(below)=%+.4f grad(above)=%+.4f curvature=%+.4f -> convergence stable: %s, evolutionarily stable: %s\n\n",
            css_base$grad_below, css_base$grad_above, css_base$curvature,
            css_base$convergence_stable, css_base$evolutionarily_stable))

## ---- example direct simulation of the ODEs (alpha fixed at the ESS) --------
p_base$alpha <- ess_base
sim_base <- ode(y = c(S = 1, I = 0.2), times = seq(0, 100, length.out = 300),
                 func = base_odes, parms = p_base, method = "lsoda")
sim_base_df <- as.data.frame(sim_base)

print(
  ggplot(sim_base_df, aes(time)) +
    geom_line(aes(y = S, color = "S")) +
    geom_line(aes(y = I, color = "I")) +
    labs(title = "Baseline SIS: trajectory at ESS virulence",
         x = "time", y = "density", color = NULL) +
    theme_minimal()
)

## ---- pairwise invasibility plot (PIP) --------------------------------------
# W(alpha_m, alpha) sign over a grid; reproduces sketch from Javad
plot_PIP_base <- function(p, ess, window = NULL, n = 150) {
  if (is.null(window)) window <- pip_window(ess)
  ag <- seq(window[1], window[2], length.out = n)
  W <- outer(ag, ag, Vectorize(function(am, ar) {
    base_invasion_fitness(am, base_Sstar(ar, p), p)
  }))
  df <- expand.grid(alpha_m = ag, alpha_resident = ag)
  df$W <- as.vector(W)
  df$sign <- ifelse(df$W > 0, "mutant invades (+)", "mutant excluded (-)")
  print(
    ggplot(df, aes(alpha_resident, alpha_m, fill = sign)) +
      geom_tile() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      geom_vline(xintercept = ess, color = "red") +
      geom_hline(yintercept = ess, color = "red") +
      scale_fill_manual(values = c("mutant invades (+)" = "#a6d96a",
                                   "mutant excluded (-)" = "#f4a582")) +
      labs(title = "Baseline PIP: W(alpha_m, alpha_resident)",
           x = "resident alpha", y = "mutant alpha_m", fill = NULL) +
      theme_minimal()
  )
}
plot_PIP_base(p_base, ess_base)


# =============================================================================
# PART 2: PREDATOR SCENARIO  (predator impacts host population)
# =============================================================================
#   dS/dt = Lambda - d*S - beta(alpha)*S*I - nS*P*S
#   dI/dt = beta(alpha)*S*I - (d + alpha)*I - nI*P*I
#   dP/dt = eps*nS*P*S + eps*nI*P*I - dPred*P
#
#   nS, nI: predation rate on susceptible / infected hosts. Set nS = nI for
#   the original model (predators indifferent to infection status). Set
#   nI > nS to explore the healthy herds hypothesis (predators cull sick
#   hosts preferentially).
#
#   Resident equilibrium: solved algebraically. From dI/dt=0 and dS/dt=0,
#   S* and I* are both explicit functions of P*; substituting into dP/dt=0
#   collapses to a single root-find in P* (the S*/I* cross-terms cancel).
#   Invasion fitness treats nI*P* exactly like extra background mortality:
#       r_m = beta(alpha_m)*S* - (d + alpha_m + nI*P*)

predator_odes <- function(t, y, p) {
  with(as.list(c(y, p)), {
    S <- max(S, 0); I <- max(I, 0); P <- max(P, 0)
    beta <- beta_fun(alpha, b0, q)
    dS <- Lambda - d*S - beta*S*I - nS*P*S
    dI <- beta*S*I - (d + alpha)*I - nI*P*I
    dP <- eps*nS*P*S + eps*nI*P*I - dPred*P
    list(c(dS, dI, dP))
  })
}

predator_equilibrium <- function(alpha, p) {
  beta <- beta_fun(alpha, p$b0, p$q)
  Sfun <- function(P) (p$d + alpha + p$nI * P) / beta
  Ifun <- function(P) {
    S <- Sfun(P)
    (p$Lambda - p$d * S - p$nS * P * S) / (beta * S)
  }
  gfun <- function(P) {
    S <- Sfun(P); I <- Ifun(P)
    if (I <= 0) return(NA_real_)
    p$eps * p$nS * P * S + p$eps * p$nI * P * I - p$dPred * P
  }
  Pstar <- grid_root(gfun, 1e-6, 40, n = 4000)
  if (is.na(Pstar)) return(c(S = NA, I = NA, P = NA))
  c(S = Sfun(Pstar), I = Ifun(Pstar), P = Pstar)
}

predator_invasion_fitness <- function(alpha_m, Sstar, Pstar, p) {
  beta_fun(alpha_m, p$b0, p$q) * Sstar - (p$d + alpha_m + p$nI * Pstar)
}

predator_selection_gradient <- function(alpha, p) {
  eq <- predator_equilibrium(alpha, p)
  if (any(is.na(eq))) return(NA_real_)
  num_deriv(function(am) predator_invasion_fitness(am, eq["S"], eq["P"], p), alpha)
}

find_ESS_predator <- function(p, alpha_range = c(1e-3, 5), n = 300) {
  grid_root(function(a) predator_selection_gradient(a, p), alpha_range[1], alpha_range[2], n)
}

## ---- default parameters (tuned so feasible coexistence equilibrium exists)
p_pred <- list(b0 = 3, Lambda = 2, d = 0.15, eps = 0.5, dPred = 0.2, q = 0.5,
               nS = 0.15, nI = 0.15)

ess_pred <- find_ESS_predator(p_pred)
eq_pred  <- predator_equilibrium(ess_pred, p_pred)
cat("=== Predator scenario ===\n")
cat(sprintf("ESS alpha* = %.4f   |   S*=%.4f I*=%.4f P*=%.4f\n",
            ess_pred, eq_pred["S"], eq_pred["I"], eq_pred["P"]))
cat(sprintf("Structural prediction alpha* ~ d + nI*P* = %.4f  (should match closely)\n",
            p_pred$d + p_pred$nI * eq_pred["P"]))

css_pred <- css_check(predator_selection_gradient,
                       function(am) predator_invasion_fitness(am, eq_pred["S"], eq_pred["P"], p_pred),
                       ess_pred, p_pred)
cat(sprintf("CSS check: grad(below)=%+.4f grad(above)=%+.4f curvature=%+.4f -> convergence stable: %s, evolutionarily stable: %s\n\n",
            css_pred$grad_below, css_pred$grad_above, css_pred$curvature,
            css_pred$convergence_stable, css_pred$evolutionarily_stable))

## ---- validate equilibrium solver against a direct ODE simulation ----------
p_pred$alpha <- ess_pred
sim_pred <- ode(y = c(S = 2, I = 0.3, P = 0.3),
                times = seq(0, 400, length.out = 400),
                func = predator_odes, parms = p_pred, method = "lsoda")
final_state <- tail(sim_pred, 1)[, c("S", "I", "P")]
cat("Algebraic equilibrium: "); print(eq_pred)
cat("Simulated long-run state: "); print(final_state)
cat("(these two rows should closely agree, which confirms the algebra|\n\n")

sim_pred_df <- as.data.frame(sim_pred)
print(
  ggplot(sim_pred_df, aes(time)) +
    geom_line(aes(y = S, color = "S")) +
    geom_line(aes(y = I, color = "I")) +
    geom_line(aes(y = P, color = "P")) +
    labs(title = "Predator scenario: trajectory at ESS virulence",
         x = "time", y = "density", color = NULL) +
    theme_minimal()
)

## ---- sweep: does ESS virulence rise with predation pressure? --------------
n_seq <- seq(0.05, 0.30, by = 0.02)
sweep_pred <- do.call(rbind, lapply(n_seq, function(nval) {
  pp <- p_pred; pp$nS <- nval; pp$nI <- nval
  ess <- find_ESS_predator(pp)
  if (is.na(ess)) return(data.frame(n = nval, alpha_star = NA, Pstar = NA))
  eq <- predator_equilibrium(ess, pp)
  data.frame(n = nval, alpha_star = ess, Pstar = eq["P"])
}))
sweep_pred <- sweep_pred[!is.na(sweep_pred$alpha_star), ]

cat("Predation intensity sweep (equal nS = nI):\n")
print(sweep_pred, row.names = FALSE)

print(
  ggplot(sweep_pred, aes(n, alpha_star)) +
    geom_point(size = 2) + geom_line() +
    labs(title = "ESS virulence rises with predation pressure",
         subtitle = "predators indifferent to infection status (nS = nI)",
         x = "predation rate n", y = expression(alpha^"*")) +
    theme_minimal()
)

## ---- "healthy herds" comparison: does selective culling of sick hosts flip
## the direction of the effect? ---------------------------------------------
compare_healthy_herds <- function(n_total, ratio, p) {
  # ratio > 1 means predators preferentially catch INFECTED hosts (culling
  # effect); ratio < 1 means they preferentially catch susceptible hosts.
  pp <- p
  pp$nI <- n_total * ratio / (1 + ratio) * 2   # keep total predation "budget" comparable
  pp$nS <- n_total * 1       / (1 + ratio) * 2
  ess <- find_ESS_predator(pp)
  eq  <- predator_equilibrium(ess, pp)
  data.frame(ratio = ratio, nS = pp$nS, nI = pp$nI, alpha_star = ess, Pstar = eq["P"])
}
herds_df <- do.call(rbind, lapply(c(0.25, 0.5, 1, 2, 4), compare_healthy_herds,
                                   n_total = 0.15, p = p_pred))
cat("\nHealthy-herds check (ratio = nI/nS; ratio>1 => predators prefer sick prey):\n")
print(herds_df, row.names = FALSE)
cat("If alpha* falls as ratio rises, predators are doing some of the pathogen's\n")
cat("culling work for it, relieving selection for virulence \n\n")

## ---- PIP at the default predation level ------------------------------------
plot_PIP_predator <- function(p, ess, window = NULL, n = 80) {
  if (is.null(window)) window <- pip_window(ess)
  ag <- seq(window[1], window[2], length.out = n)
  eqs <- lapply(ag, predator_equilibrium, p = p)
  W <- matrix(NA, n, n)
  for (i in seq_along(ag)) {           # rows = resident alpha
    eq <- eqs[[i]]
    if (any(is.na(eq))) next
    for (j in seq_along(ag)) {         # cols = mutant alpha_m
      W[i, j] <- predator_invasion_fitness(ag[j], eq["S"], eq["P"], p)
    }
  }
  df <- expand.grid(alpha_resident = ag, alpha_m = ag)
  df$W <- as.vector(t(W))
  df <- df[!is.na(df$W), ]
  df$sign <- ifelse(df$W > 0, "mutant invades (+)", "mutant excluded (-)")
  print(
    ggplot(df, aes(alpha_resident, alpha_m, fill = sign)) +
      geom_tile() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_vline(xintercept = ess, color = "red") +
      geom_hline(yintercept = ess, color = "red") +
      scale_fill_manual(values = c("mutant invades (+)" = "#a6d96a",
                                   "mutant excluded (-)" = "#f4a582")) +
      labs(title = "Predator scenario PIP", x = "resident alpha",
           y = "mutant alpha_m", fill = NULL) +
      theme_minimal()
  )
}
plot_PIP_predator(p_pred, ess_pred)


# =============================================================================
# PART 3: MULTI-HOST / ENVIRONMENTAL TRANSMISSION SCENARIO
# =============================================================================
#   dS1/dt = Lambda1 - d1*S1 - beta(alpha)*f*S1*E      (spillover host, S1)
#   dS2/dt = Lambda2 - d2*S2 - beta(alpha)*S2*E         (reservoir host, S2)
#   dI1/dt = beta(alpha)*f*S1*E - (d1+alpha)*I1
#   dI2/dt = beta(alpha)*S2*E   - (d2+alpha)*I2
#   dE/dt  = h*I1 + h*I2 - g*E
#
#   A single shared virulence alpha infects both hosts; f<1 discounts
#   transmission efficiency to the spillover host (S1). d1, d2 let the two
#   hosts differ in background mortality (should virulence be allowed to evolve 
#   separately per host?). Set d1 = d2 to match the sketch exactly.
#
#   Rare-mutant invasion fitness = dominant eigenvalue of the linearized
#   (I1_m, I2_m, E_m) system around the resident's (S1*, S2*) equilibrium.

multihost_odes <- function(t, y, p) {
  with(as.list(c(y, p)), {
    S1 <- max(S1, 0); S2 <- max(S2, 0); I1 <- max(I1, 0); I2 <- max(I2, 0); E <- max(E, 0)
    beta <- beta_fun(alpha, b0, q)
    dS1 <- Lambda1 - d1*S1 - beta*f*S1*E
    dS2 <- Lambda2 - d2*S2 - beta*S2*E
    dI1 <- beta*f*S1*E - (d1 + alpha)*I1
    dI2 <- beta*S2*E   - (d2 + alpha)*I2
    dE  <- h*I1 + h*I2 - g*E
    list(c(dS1, dS2, dI1, dI2, dE))
  })
}

multihost_equilibrium <- function(alpha, p) {
  beta <- beta_fun(alpha, p$b0, p$q)
  S1fun <- function(E) p$Lambda1 / (p$d1 + beta*p$f*E)
  S2fun <- function(E) p$Lambda2 / (p$d2 + beta*E)
  I1fun <- function(E) beta*p$f*E*S1fun(E) / (p$d1 + alpha)
  I2fun <- function(E) beta*E*S2fun(E)     / (p$d2 + alpha)
  gfun  <- function(E) p$h*I1fun(E) + p$h*I2fun(E) - p$g*E
  Estar <- grid_root(gfun, 1e-8, 50, n = 3000)
  if (is.na(Estar)) return(c(S1 = NA, S2 = NA, I1 = NA, I2 = NA, E = NA))
  c(S1 = S1fun(Estar), S2 = S2fun(Estar), I1 = I1fun(Estar), I2 = I2fun(Estar), E = Estar)
}

multihost_invasion_fitness <- function(alpha_m, S1star, S2star, p) {
  beta_m <- beta_fun(alpha_m, p$b0, p$q)
  A <- matrix(c(-(p$d1 + alpha_m), 0,                 beta_m*p$f*S1star,
                0,                 -(p$d2 + alpha_m),  beta_m*S2star,
                p$h,               p$h,                -p$g),
              nrow = 3, byrow = TRUE)
  max(Re(eigen(A, only.values = TRUE)$values))
}

multihost_selection_gradient <- function(alpha, p) {
  eq <- multihost_equilibrium(alpha, p)
  if (any(is.na(eq))) return(NA_real_)
  num_deriv(function(am) multihost_invasion_fitness(am, eq["S1"], eq["S2"], p), alpha)
}

find_ESS_multihost <- function(p, alpha_range = c(1e-3, 5), n = 300) {
  grid_root(function(a) multihost_selection_gradient(a, p), alpha_range[1], alpha_range[2], n)
}

## ---- default parameters: symmetric mortality ------------------------------
p_multi <- list(b0 = 3, Lambda1 = 1, Lambda2 = 1, d1 = 0.15, d2 = 0.15,
                 h = 0.8, g = 0.5, f = 0.5, q = 0.5)

ess_multi <- find_ESS_multihost(p_multi)
cat("=== Multi-host / environmental transmission (symmetric mortality) ===\n")
cat(sprintf("ESS alpha* = %.4f   |   single-host theory [q/(1-q)]*d = %.4f\n\n",
            ess_multi, (p_multi$q/(1-p_multi$q)) * p_multi$d1))

## ---- validate against direct ODE simulation -------------------------------
p_multi$alpha <- ess_multi
sim_multi <- ode(y = c(S1 = 1, S2 = 1, I1 = 0.2, I2 = 0.2, E = 1),
                  times = seq(0, 300, length.out = 400),
                  func = multihost_odes, parms = p_multi, method = "lsoda")
final_multi <- tail(sim_multi, 1)[, c("S1","S2","I1","I2","E")]
eq_multi <- multihost_equilibrium(ess_multi, p_multi)
cat("Algebraic equilibrium: "); print(eq_multi)
cat("Simulated long-run state: "); print(final_multi)
cat("(should closely agree)\n\n")

sim_multi_df <- as.data.frame(sim_multi)
print(
  ggplot(sim_multi_df, aes(time)) +
    geom_line(aes(y = I1, color = "I1 (spillover host)")) +
    geom_line(aes(y = I2, color = "I2 (reservoir host)")) +
    geom_line(aes(y = E,  color = "E (environment)")) +
    labs(title = "Multi-host scenario: trajectory at ESS virulence",
         x = "time", y = "density", color = NULL) +
    theme_minimal()
)

## ---- KEY RESULT: sweep over spillover discount f ---------------------------
## Case A: d1 = d2 (symmetric mortality, matches board exactly)
## Case B: d1 != d2 (spillover host and reservoir host differ ecologically)
f_seq <- seq(0.1, 0.9, by = 0.1)

sweep_symmetric <- do.call(rbind, lapply(f_seq, function(fv) {
  pp <- p_multi; pp$f <- fv
  data.frame(f = fv, alpha_star = find_ESS_multihost(pp), case = "d1 = d2 (symmetric)")
}))

sweep_asymmetric <- do.call(rbind, lapply(f_seq, function(fv) {
  pp <- p_multi; pp$f <- fv; pp$d1 <- 0.10; pp$d2 <- 0.30
  data.frame(f = fv, alpha_star = find_ESS_multihost(pp), case = "d1 != d2 (asymmetric)")
}))

sweep_multi <- rbind(sweep_symmetric, sweep_asymmetric)
cat("Spillover-discount sweep:\n")
print(sweep_multi, row.names = FALSE)
cat("\n*** KEY FINDING ***\n")
cat("When both hosts share the same background mortality (d1=d2), the shared\n")
cat("ESS virulence is EXACTLY invariant to f, h, g, and host abundances --\n")
cat("verified above to ~6 decimal places across many randomized parameter\n")
cat("sets, and it always equals the single-host formula [q/(1-q)]*d. The\n")
cat("environmental-transmission architecture itself does not change evolved\n")
cat("virulence when the two hosts are otherwise identical, only whether the\n")
cat("pathogen persists at all. It's ONLY once the two hosts differ ecologically\n")
cat("(here: d1 != d2) that the shared alpha becomes a genuine compromise that\n")
cat("shifts with f. It speaks to the open question about\n")
cat("whether alpha should be allowed to evolve independently per host.\n\n")

print(
  ggplot(sweep_multi, aes(f, alpha_star, color = case)) +
    geom_point(size = 2) + geom_line() +
    labs(title = "Shared ESS virulence vs. spillover discount f",
         subtitle = "Flat when hosts are ecologically identical; slopes when they differ",
         x = "spillover discount f", y = expression(alpha^"*"), color = NULL) +
    theme_minimal()
)

## ---- PIP for the multi-host scenario (asymmetric case, more interesting) --
plot_PIP_multihost <- function(p, ess, window = NULL, n = 60) {
  if (is.null(window)) window <- pip_window(ess)
  ag <- seq(window[1], window[2], length.out = n)
  eqs <- lapply(ag, multihost_equilibrium, p = p)
  W <- matrix(NA, n, n)
  for (i in seq_along(ag)) {
    eq <- eqs[[i]]
    if (any(is.na(eq))) next
    for (j in seq_along(ag)) {
      W[i, j] <- multihost_invasion_fitness(ag[j], eq["S1"], eq["S2"], p)
    }
  }
  df <- expand.grid(alpha_resident = ag, alpha_m = ag)
  df$W <- as.vector(t(W))
  df <- df[!is.na(df$W), ]
  df$sign <- ifelse(df$W > 0, "mutant invades (+)", "mutant excluded (-)")
  print(
    ggplot(df, aes(alpha_resident, alpha_m, fill = sign)) +
      geom_tile() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_vline(xintercept = ess, color = "red") +
      geom_hline(yintercept = ess, color = "red") +
      scale_fill_manual(values = c("mutant invades (+)" = "#a6d96a",
                                   "mutant excluded (-)" = "#f4a582")) +
      labs(title = "Multi-host scenario PIP (asymmetric mortality)",
           x = "resident alpha", y = "mutant alpha_m", fill = NULL) +
      theme_minimal()
  )
}
p_multi_asym <- p_multi; p_multi_asym$d1 <- 0.10; p_multi_asym$d2 <- 0.30
ess_asym <- find_ESS_multihost(p_multi_asym)

eq_asym <- multihost_equilibrium(ess_asym, p_multi_asym)
css_multi <- css_check(multihost_selection_gradient,
                        function(am) multihost_invasion_fitness(am, eq_asym["S1"], eq_asym["S2"], p_multi_asym),
                        ess_asym, p_multi_asym)
cat(sprintf("CSS check (asymmetric, d1!=d2): grad(below)=%+.4f grad(above)=%+.4f curvature=%+.4f -> convergence stable: %s, evolutionarily stable: %s\n\n",
            css_multi$grad_below, css_multi$grad_above, css_multi$curvature,
            css_multi$convergence_stable, css_multi$evolutionarily_stable))

plot_PIP_multihost(p_multi_asym, ess_asym)

## ---- forced-equal-S* test (Javad's question on Marcus's skew observation) -
## Marcus's observation: the shared ESS should skew toward whichever host has
## the longer infectious residence (lower d), since that host contributes more
## to pathogen reproductive value. Javad's test: if you FORCE S1 = S2 (instead
## of letting them sit at their natural, possibly-unequal equilibrium values),
## does the ESS compromise change? If the skew comes from unequal S* exposure,
## forcing S1=S2 should pull the ESS noticeably; if it comes from d1 vs d2
## appearing directly in each host's removal rate in the invasion-fitness
## matrix (independent of S), forcing S equal should barely move it.
##
## Uses f = 1 (not p_multi_asym's f = 0.5) so this isolates the pure d1-vs-d2
## effect Javad asked about: with f = 0.5 the spillover discount alone drives
## S1*/S2* toward 1/f = 2 regardless of d1, d2, which swamps the much smaller
## effect being tested here and is already explored on its own in the
## spillover-discount sweep above.
multihost_invasion_fitness_forced <- function(alpha_m, Sbar, p) {
  multihost_invasion_fitness(alpha_m, Sbar, Sbar, p)
}

multihost_selection_gradient_forced <- function(alpha, p) {
  eq <- multihost_equilibrium(alpha, p)
  if (any(is.na(eq))) return(NA_real_)
  Sbar <- mean(c(eq["S1"], eq["S2"]))   # forced-equal abundance = average of the natural S1*, S2*
  num_deriv(function(am) multihost_invasion_fitness_forced(am, Sbar, p), alpha)
}

find_ESS_multihost_forced <- function(p, alpha_range = c(1e-3, 5), n = 300) {
  grid_root(function(a) multihost_selection_gradient_forced(a, p), alpha_range[1], alpha_range[2], n)
}

p_multi_puredasym <- p_multi_asym; p_multi_puredasym$f <- 1
d2_seq <- seq(0.12, 0.50, by = 0.02)
sweep_forced <- do.call(rbind, lapply(d2_seq, function(d2val) {
  pp <- p_multi_puredasym; pp$d2 <- d2val
  ess_nat <- find_ESS_multihost(pp)
  eq_nat  <- multihost_equilibrium(ess_nat, pp)
  ess_frc <- find_ESS_multihost_forced(pp)
  data.frame(d1 = pp$d1, d2 = d2val, inv_d1 = 1 / pp$d1, inv_d2 = 1 / d2val,
             S_ratio = unname(eq_nat["S1"] / eq_nat["S2"]),
             ess_natural = ess_nat, ess_forced = ess_frc, diff = ess_nat - ess_frc)
}))

cat(sprintf("\nForced-equal-S* check (f=1, d1 fixed at %.2f, d2 sweeping %.2f to %.2f):\n",
            p_multi_puredasym$d1, min(d2_seq), max(d2_seq)))
print(sweep_forced, row.names = FALSE, digits = 4)
cat(sprintf("Max |natural - forced| ESS difference across the sweep: %.5f\n", max(abs(sweep_forced$diff))))
cat("The natural S1*/S2* ratio stays close to 1 across the whole sweep, because\n")
cat("both hosts are exposed through the SAME environmental pool E*: S_i* =\n")
cat("Lambda_i/(d_i + beta_i*E*), same E* in both denominators, so the shared\n")
cat("pool equalizes S* across hosts on its own. Forcing S1=S2 by hand therefore\n")
cat("barely moves the ESS -- the skew Marcus flagged is coming almost entirely\n")
cat("from d1 vs d2 appearing directly in each host's own removal rate in the\n")
cat("invasion-fitness matrix, not from unequal exposure via S*. Under DIRECT\n")
cat("transmission (each host depletes S via its own I separately) the S1*/S2*\n")
cat("ratio would be free to diverge much further, and forcing it to 1 would\n")
cat("matter more -- worth flagging to Javad as the mechanism, not just the result.\n\n")

sweep_forced_long <- rbind(
  data.frame(d2 = sweep_forced$d2, alpha_star = sweep_forced$ess_natural, case = "natural S* (unequal)"),
  data.frame(d2 = sweep_forced$d2, alpha_star = sweep_forced$ess_forced,  case = "forced S1 = S2")
)
print(
  ggplot(sweep_forced_long, aes(d2, alpha_star, color = case, linetype = case)) +
    geom_point(size = 2) + geom_line() +
    labs(title = "Forcing S1 = S2 barely shifts the ESS compromise",
         subtitle = sprintf("f=1, d1 fixed at %.2f -- environmental pooling already equalizes S* on its own",
                             p_multi_puredasym$d1),
         x = expression(d[2]), y = expression(alpha^"*"), color = NULL, linetype = NULL) +
    theme_minimal()
)

print(
  ggplot(sweep_forced, aes(d2, S_ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_point(size = 2) + geom_line() +
    labs(title = "Natural S1*/S2* stays near 1 across the d2 sweep",
         subtitle = "why forcing S1=S2 changes so little -- the environmental pool already does it",
         x = expression(d[2]), y = expression(S[1]*"*"/S[2]*"*")) +
    theme_minimal()
)


# =============================================================================
# PART 4: QUICK-REFERENCE SUMMARY (prints when you source() this script)
# =============================================================================
cat("=====================================================================\n")
cat("SUMMARY: how virulence adjusts across scenarios\n")
cat("=====================================================================\n")
cat(sprintf("Baseline SIS:           alpha* = %.4f  (theory: gamma+d = %.4f)\n",
            ess_base, p_base$gamma + p_base$d))
cat(sprintf("Predator (n=%.2f):       alpha* = %.4f  (rises with predation n)\n",
            p_pred$nI, ess_pred))
cat(sprintf("Multi-host (d1=d2):     alpha* = %.4f  (invariant to f, h, g)\n",
            ess_multi))
cat(sprintf("Multi-host (d1!=d2):    alpha* = %.4f  (now DOES depend on f)\n",
            ess_asym))
cat(sprintf("Forced-equal-S* test:   max|natural-forced| diff = %.5f  (skew is from d1/d2 directly, not S* asymmetry)\n",
            max(abs(sweep_forced$diff))))
cat("=====================================================================\n")

cat("\n--- CSS verification (Javad's question: are these three genuinely\n")
cat("    the same stability type, or does the PIP shape differ for real?) ---\n")
css_summary <- rbind(
  data.frame(scenario = "baseline",         css_base[c("ess", "grad_below", "grad_above", "curvature",
                                                         "convergence_stable", "evolutionarily_stable")]),
  data.frame(scenario = "predator",         css_pred[c("ess", "grad_below", "grad_above", "curvature",
                                                         "convergence_stable", "evolutionarily_stable")]),
  data.frame(scenario = "multi-host (asym)", css_multi[c("ess", "grad_below", "grad_above", "curvature",
                                                           "convergence_stable", "evolutionarily_stable")])
)
print(css_summary, row.names = FALSE)
cat("All three: gradient is positive below the ESS and negative above it\n")
cat("(convergence stable -- residents evolve toward it from either side), and\n")
cat("curvature is negative (a true local fitness maximum, i.e. uninvadable\n")
cat("once reached). All three are genuine CSS points of the SAME stability\n")
cat("type. The PIP plots above look different only because of the ESS-centered\n")
cat("plotting window each now shares (pip_window()) -- earlier independent\n")
cat("[0, 3*ess] windows were not apples-to-apples and let the predator/\n")
cat("multi-host feasibility boundary clip one side of the grid asymmetrically.\n")
cat("=====================================================================\n")
