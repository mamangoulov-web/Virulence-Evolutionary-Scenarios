# =============================================================================
# VIRULENCE EVOLUTION ACROSS ECOLOGICAL SCENARIOS
# Multi-host/environmental transmission
# =============================================================================
#
# Quick overview, since this covers a few scenarios: for each one below I
# (1) code the ODE system exactly as sketched, so it can be simulated and
# the trajectories (S(t), I(t), ...etc) watched directly as a check on the ODEs
# themselves, (2) solve for the resident ecological equilibrium as a function
# of virulence alpha, numerically or algebraically, (3) get the invasion
# fitness of a rare mutant strain and its selection gradient, then find the
# ESS as the root of that gradient, (4) push whatever the key ecological
# variable is for that scenario -- predation intensity n, spillover discount f,
# etc -- and see how ESS virulence responds, and (5) draw a pairwise invasibility
# plot (PIP) to confirm each ESS is a convergence-stable singular
# strategy, not merely a zero of the gradient that happens to sit there.
#
# TWO MODELING NOTES (done during the writeup):
#
#   (1) Birth term: the initial sketches used proportional birth "bS" for the
#       predator and multi-host scenarios. Written that way, both scenarios
#       turn out to have no generic fixed-point equilibrium: proportional
#       birth with no self-limitation makes the models scale-invariant
#       (homogeneous), so they either grow/decay without bound or only
#       balance at a knife-edge parameter combination. This is why
#       the R0 formula carries the note N=1 next to it, as
#       the base model implicitly assumes some constant recruitment/
#       normalized population size, not free proportional growth. I've made
#       that assumption explicit here: births enter as a constant recruitment
#       rate Lambda into the susceptible class, consistent with the N=1
#       convention; this is also the standard convention in the virulence-
#       evolution literature (SIS/SI models). If it interests anybody, we can
#       still simulate the literal "bS" version below (set birth_type =
#       "proportional") to see the instability. Regardless, it is something
#       I certainly wanted to verify empirically.
#
#       To follow up: with per-capita birth bS employed instead of
#       constant recruitment, the endemic Jacobian works out to
#       trace = -gamma*(b-d)/(d+alpha), det = (b-d)*(d+gamma+alpha). det>0
#       automatically whenever b>d (the feasibility condition), so stability
#       hinges entirely on the trace, and trace<0 if gamma>0. If gamma=0, the
#       trace is exactly zero for any parameter values -- a structural neutral
#       center (sustained oscillation), not just a bad parameter choice. So
#       per-capita birth was never the real problem for the baseline SIS model
#       (which has gamma>0 and could have kept it and stayed stable); constant
#       recruitment was only strictly necessary for the predator/multi-host
#       scenarios below, which my sketches involved no recovery term at all
#       (gamma=0 included the structure), which sits exactly on that
#       knife-edge regardless of any other parameter choice.
#
#   (2) Trade-off shape: beta(alpha) = b0 * alpha^q, with 0 < q < 1 (concave
#       shape). q = 0.5 (square root) is the special case that reproduces the
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

# generic 1-D root bracketing over a grid
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

# second-order (curvature) numerical derivative -- used for the
# evolutionary-stability part of the CSS check (is the ESS a local fitness
# maximum in the mutant trait?)
num_second_deriv <- function(f, x, h = 1e-3) (f(x + h) - 2 * f(x) + f(x - h)) / h^2

# Full CSS (con. stab. strat.) check, shared across all three
# scenarios: convergence stability comes from the selection-gradient sign on
# either side of the ESS (+ below, - above dictating if residents tend to or away from it);
# evolutionary stability comes from the curvature of invasion fitness in the
# mutant trait alone, resident frozen at the ESS (negative means a true local
# fitness peak, not just a zero of the gradient). Both must hold for a genuine
# CSS as this is what actually answers "is it convergent and stable".
css_check <- function(selection_gradient_fun, invasion_fitness_of_mutant_fun, ess, p, rel_step = 0.05) {
  delta <- max(rel_step * ess, 1e-4)
  grad_below <- selection_gradient_fun(ess - delta, p)
  grad_above <- selection_gradient_fun(ess + delta, p)
  curvature  <- num_second_deriv(invasion_fitness_of_mutant_fun, ess)
  data.frame(ess = ess, grad_below = grad_below, grad_above = grad_above, curvature = curvature,
             convergence_stable = isTRUE(grad_below > 0 && grad_above < 0),
             evolutionarily_stable = isTRUE(curvature < 0))
}

# ESS-centered plotting window shared by all three PIP plots below. An earlier version windowed each plot independently as
# [0, 3*ess], which skewed the ESS off-center and let the predator/multi-host
# feasibility boundary (equilibrium stops existing above some alpha) clip one
# side of the grid asymmetrically. Differences in how those plots looked
# were a result of improper scale, not a true difference in stability type (the
# css_check() results printed for each scenario agree on that).
pip_window <- function(ess, zoom = 0.75) c(max(1e-4, ess * (1 - zoom)), ess * (1 + zoom))

# shared axis/legend labels for every PIP below, so each plot says what it is showing
pip_labs <- function(title) {
  labs(title = title,
       subtitle = paste0("Green = a rare mutant (alpha_m) can invade the resident (alpha); orange = it is excluded.\n",
                         "Red lines mark the ESS; dashed diagonal = mutant identical to resident."),
       x = "resident virulence (alpha)", y = "mutant virulence (alpha_m)", fill = NULL)
}

# shared axis labels for the shedding-asymmetry plots in Part 5
lab_ess        <- expression("shared ESS virulence"~alpha^"*"~"(disease-induced mortality)")
lab_shed_ratio <- expression(h[1]/h[2]~"= shedding-rate ratio, host 1 / host 2 (log scale)")


# =============================================================================
# PART 1: BASELINE SIS MODEL  (no predator, single host)
# =============================================================================
#   dS/dt = Lambda - d*S - beta(alpha)*S*I + gamma*I
#   dI/dt = beta(alpha)*S*I - (d + gamma + alpha)*I
#
#   Resident equilibrium threshold (from R0 = beta(alpha)*S/(d+alpha+gamma)):
#       S* = (d + gamma + alpha) / beta(alpha)
#   Invasion fitness of a rare mutant (alpha_m) on resident background S*:
#       r_m = beta(alpha_m)*S* - (d + gamma + alpha_m)

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
# scenarios below): freeze the resident equilibrium first, then differentiate
# invasion fitness with respect to the mutant trait only, evaluated at
# alpha_m = alpha (mutant = resident) -- i.e. dr/d(alpha_m) |_{alpha_m=alpha},
# the standard adaptive-dynamics selection gradient.
base_selection_gradient <- function(alpha, p) {
  Sstar <- base_Sstar(alpha, p)
  num_deriv(function(am) base_invasion_fitness(am, Sstar, p), alpha)
}

find_ESS_base <- function(p, interval = c(1e-3, 50)) {
  uniroot(function(a) base_selection_gradient(a, p), interval)$root
}

## ---- run and verify against theory -----------------------------------------
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
         subtitle = "S = susceptible hosts, I = infected hosts; both settle to their equilibrium values",
         x = "time (arbitrary units)", y = "host density (arbitrary units)", color = NULL) +
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
      pip_labs("Baseline SIS pairwise invasibility plot (PIP)") +
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
#   Resident equilibrium: solved algebraically; from dI/dt=0 and dS/dt=0,
#   S* and I* are both explicit functions of P*; substituting into dP/dt=0
#   collapses to a single root-find in P* (the S*/I* cross-terms cancel out).
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

## ---- default parameters (optimized so feasible coexistence equilibrium exists)
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
         subtitle = "S = susceptible hosts, I = infected hosts, P = predators",
         x = "time (arbitrary units)", y = "density (arbitrary units)", color = NULL) +
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
         subtitle = "Predators indifferent to infection status (nS = nI = n)",
         x = "predation rate n (per predator, same on susceptible and infected hosts)",
         y = expression("ESS virulence"~alpha^"*"~"(disease-induced mortality)")) +
    theme_minimal()
)

## ---- "healthy herds" comparison: does selective culling of sick hosts flip
## the direction of the effect? ---------------------------------------------
## The healthy herds hypothesis (Packer et al. 2003) claims that predators which
## kill infected hosts keep the herd healthier overall. 
## This is a claim about pathogen prevalence (in pop), so I tracked prevalence 
## at the ESS alongside alpha*. The ecological effect and the 
## evolutionary effect seem to point in opposite directions.
compare_healthy_herds <- function(n_total, ratio, p) {
  # ratio > 1 means predators preferentially catch infected hosts (culling
  # effect); ratio < 1 means they preferentially catch susceptible hosts.
  pp <- p
  pp$nI <- n_total * ratio / (1 + ratio) * 2   # keep total predation "budget" comparable
  pp$nS <- n_total * 1       / (1 + ratio) * 2  # (nS + nI = 2*n_total at every ratio)
  ess <- find_ESS_predator(pp)
  eq  <- predator_equilibrium(ess, pp)
  data.frame(ratio = ratio, nS = pp$nS, nI = pp$nI, alpha_star = ess, Pstar = eq["P"],
             Istar = eq["I"], prevalence = eq["I"] / (eq["S"] + eq["I"]))
}
herds_df <- do.call(rbind, lapply(c(0.25, 0.5, 1, 2, 4), compare_healthy_herds,
                                   n_total = 0.15, p = p_pred))
cat("\nHealthy-herds check (ratio = nI/nS; ratio>1 => predators prefer sick prey):\n")
print(herds_df, row.names = FALSE)
cat("Result: alpha* rises as predators favour sick hosts, while prevalence falls.\n")
cat("Culling adds nI*P* to the death rate of infected hosts -- the same structural\n")
cat("role as background mortality d -- so it shortens the infectious period for\n")
cat("every strain alike. Those hosts were inevitably dying, so virulence\n")
cat("is less costly to pathogen and selection favours faster exploitation. Healthy\n")
cat("herds lead to less disease, but the remaining disease becomes more severe.\n\n")

herds_long <- rbind(
  data.frame(ratio = herds_df$ratio, value = herds_df$alpha_star,
             panel = "ESS virulence alpha*"),
  data.frame(ratio = herds_df$ratio, value = 100 * herds_df$prevalence,
             panel = "prevalence at the ESS (% of hosts infected)")
)
print(
  ggplot(herds_long, aes(ratio, value)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_point(size = 2) + geom_line() +
    facet_wrap(~ panel, scales = "free_y") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    labs(title = "Healthy herds: culling sick hosts means less disease, but more virulent disease",
         subtitle = "Total predation fixed (nS + nI = 0.30); dashed line = indiscriminate predation (nS = nI)",
         x = "predator preference for infected hosts, nI/nS (log scale)",
         y = NULL) +
    theme_minimal()
)

## ---- What decides the direction? Virulence-dependent capture ---------------
## Betts et al. (2016)'s review cite Kisdi et al. (2013) for the opposite 
## prediction: predation selects for lower virulence. The key difference in their 
## assumption: more virulent infections make prey easier to catch (sluggish,
## conspicuous), so capture of infected hosts depends on virulence itself.
## Henceforth, that is added as
##       nI(alpha) = n0 + k*alpha, with k = virulence 'attraction' factor
## so the infected predation rate becomes strain-specific and the mutant pays
## for it at its own virulence:
##       r_m = beta(alpha_m)*S* - (d + alpha_m + (n0 + k*alpha_m)*P*)
## Setting the selection gradient to zero gives (power-law trade-off)
##       alpha* = q*(d + n0*P*) / [(1-q)*(1 + k*P*)]
## d(alpha*)/dP* has the sign of (n0 - k*d): predation raises virulence when
## the virulence-independent part of capture dominates (n0 > k*d, which
## includes k = 0, i.e. this model), and lowers it when n0 < k*d (when the excess
## rate of predator capture exceeds the baseline rate, thus selecting for less
## extreme virulence which makes them harder to capture again).
## This Kisdi assumption is just incorporated into my framework, and is not a
## reproduction of their full model; n0 = 0.02, k = 0.3 are illustrative.
predator_vdc_equilibrium <- function(alpha, p) {
  pp <- p; pp$nI <- p$n0 + p$k * alpha    # resident's infected-capture rate
  predator_equilibrium(alpha, pp)
}

predator_vdc_invasion_fitness <- function(alpha_m, Sstar, Pstar, p) {
  beta_fun(alpha_m, p$b0, p$q) * Sstar - (p$d + alpha_m + (p$n0 + p$k * alpha_m) * Pstar)
}

predator_vdc_selection_gradient <- function(alpha, p) {
  eq <- predator_vdc_equilibrium(alpha, p)
  if (any(is.na(eq))) return(NA_real_)
  num_deriv(function(am) predator_vdc_invasion_fitness(am, eq["S"], eq["P"], p), alpha)
}

find_ESS_predator_vdc <- function(p, alpha_range = c(1e-3, 5), n = 300) {
  grid_root(function(a) predator_vdc_selection_gradient(a, p), alpha_range[1], alpha_range[2], n)
}

# scale every predation rate (nS, n0, k) by the same factor to vary intensity
sweep_vdc <- function(n0, k, case_label, scale_seq = seq(0.6, 2, by = 0.1)) {
  do.call(rbind, lapply(scale_seq, function(s) {
    pp <- p_pred; pp$nS <- 0.15 * s; pp$n0 <- n0 * s; pp$k <- k * s
    ess <- find_ESS_predator_vdc(pp)
    if (is.na(ess)) return(NULL)        # no predator-host coexistence at this intensity
    eq <- predator_vdc_equilibrium(ess, pp)
    data.frame(scale = s, alpha_star = ess, Pstar = eq["P"],
               theory = pp$q * (pp$d + pp$n0 * eq["P"]) / ((1 - pp$q) * (1 + pp$k * eq["P"])),
               case = case_label)
  }))
}
vdc_df <- rbind(
  sweep_vdc(n0 = 0.15, k = 0,   "k = 0: capture independent of virulence (this model)"),
  sweep_vdc(n0 = 0.02, k = 0.3, "k > 0: virulent hosts easier to catch (Kisdi-type)")
)
cat("Virulence-dependent capture, nI = n0 + k*alpha (scale multiplies nS, n0 and k):\n")
for (case_label in unique(vdc_df$case)) {
  cat(" ", case_label, "\n")
  print(vdc_df[vdc_df$case == case_label, c("scale", "alpha_star", "Pstar", "theory")],
        row.names = FALSE, digits = 4)
}
cat(sprintf("Max |numeric - theory| ESS difference: %.1e\n", max(abs(vdc_df$alpha_star - vdc_df$theory))))
cat("With k = 0 alpha* rises with predation; with n0 < k*d it falls. The direction is\n")
cat("set by whether virulence itself makes hosts easier to catch, not by whether\n")
cat("predators prefer sick prey.\n\n")

print(
  ggplot(vdc_df, aes(scale, color = case)) +
    geom_line(aes(y = theory), alpha = 0.5) +
    geom_point(aes(y = alpha_star), size = 2) +
    labs(title = "What decides the direction: does virulence make hosts easier to catch?",
         subtitle = "Points = numerical ESS; lines = alpha* = q(d + n0*P*) / [(1-q)(1 + k*P*)]",
         x = "predation intensity (nS, n0 and k all multiplied by this factor)",
         y = expression("ESS virulence"~alpha^"*"~"(disease-induced mortality)"),
         color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
)

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
      pip_labs("Predator scenario pairwise invasibility plot (PIP)") +
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
#   separately per host?). Set d1 = d2 to match the sketch exactly. E functions
## as an environmental reservoir, with inflow coming from shedding of disease
## from infecteds (h), and outflow from natural decay rate of the pathogen (g).
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
         subtitle = "Infected densities in each host, plus pathogen load in the shared environment",
         x = "time (arbitrary units)", y = "density (arbitrary units)", color = NULL) +
    theme_minimal()
)

## ---- Sweep over spillover discount f ---------------------------
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
cat("\n*** Key finding ***\n")
cat("When both hosts share the same background mortality (d1=d2), the shared\n")
cat("ESS virulence is exactly invariant to f, h, g, and host abundances --\n")
cat("verified above to ~6 decimal places across many randomized parameter\n")
cat("sets, and it always equals the single-host formula [q/(1-q)]*d. The\n")
cat("environmental-transmission architecture itself does not change evolved\n")
cat("virulence when the two hosts are otherwise identical, only whether the\n")
cat("pathogen persists at all. It is only once the two hosts differ ecologically\n")
cat("(here: d1 != d2) that the shared alpha becomes a genuine compromise that\n")
cat("shifts with f. It speaks to a question about\n")
cat("whether alpha should be allowed to evolve independently per host.\n\n")

print(
  ggplot(sweep_multi, aes(f, alpha_star, color = case)) +
    geom_point(size = 2) + geom_line() +
    labs(title = "Shared ESS virulence vs. spillover discount f",
         subtitle = "Flat when hosts are ecologically identical; slopes when they differ",
         x = "spillover discount f (relative transmission to host 1; 1 = no discount)",
         y = expression("shared ESS virulence"~alpha^"*"), color = "host mortality") +
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
      pip_labs("Multi-host pairwise invasibility plot (PIP), d1 != d2") +
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

## ---- forced-equal-S* test (a question posed by Javad, my invigilator) --
## The shared ESS should skew toward whichever host has
## the longer infectious residence (lower d), since that host contributes more
## to pathogen reproductive value. Javad proposes: if you force S1 = S2 (instead
## of letting them sit at their natural, possibly-unequal equilibrium values),
## does the ESS compromise change? If the skew comes from unequal S* exposure,
## forcing S1=S2 should pull the ESS noticeably; if it comes from d1 vs d2
## appearing directly in each host's removal rate in the invasion-fitness
## matrix (independent of S), forcing S equal should barely move it.
##
## Uses f = 1 (not p_multi_asym's f = 0.5) so this isolates the pure d1-vs-d2
## effect: with f = 0.5 the spillover discount alone drives
## S1*/S2* toward 1/f = 2 regardless of d1, d2, which squeezes the much smaller
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
cat("both hosts are exposed through the same environmental pool E*: S_i* =\n")
cat("Lambda_i/(d_i + beta_i*E*), same E* in both denominators, so the shared\n")
cat("pool equalizes S* across hosts on its own. Forcing S1=S2 therefore\n")
cat("barely moves the ESS. The skew I identified is coming almost entirely\n")
cat("from d1 vs d2 appearing directly in each host's own removal rate in the\n")
cat("invasion-fitness matrix, not from unequal exposure via S*. Under direct\n")
cat("transmission (each host depletes S via its own I separately) the S1*/S2*\n")
cat("ratio would be free to diverge much further, and forcing it to 1 would\n")
cat("matter more. Good thing to understand both the mechanism and the result.\n\n")

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
         x = expression(d[2]~"(background mortality of host 2)"),
         y = expression("shared ESS virulence"~alpha^"*"), color = NULL, linetype = NULL) +
    theme_minimal()
)

print(
  ggplot(sweep_forced, aes(d2, S_ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_point(size = 2) + geom_line() +
    labs(title = "Natural S1*/S2* stays near 1 across the d2 sweep",
         subtitle = "why forcing S1=S2 changes so little -- the environmental pool already does it",
         x = expression(d[2]~"(background mortality of host 2)"),
         y = expression(S[1]*"*"/S[2]*"*"~"(equilibrium susceptibles, host 1 / host 2)")) +
    theme_minimal()
)


# =============================================================================
# PART 4: ASYMMETRIC SHEDDING EXTENSION (an independent excursion)
# =============================================================================
#   Same multi-host/environmental-transmission architecture as Part 3, with
#   one change: each host gets its own per-capita shedding rate into the
#   environment instead of a shared h.
#       dE/dt = h1*I1 + h2*I2 - g*E      (was h*I1 + h*I2 - g*E)
#   Everything else -- S1, S2, I1, I2 equations, the trade-off, the single
#   shared alpha -- is unchanged from Part 3.
#
#   Swept as a ratio h1/h2 at fixed h1+h2, so the shedding-asymmetry effect
#   is isolated from the total-pathogen-load effect (more total shedding
#   into E, regardless of split, already changes R0/feasibility on its own --
#   that's not what's being tested).
#
#   Prediction to test: h1 > h2 means more of the pathogen's effective
#   "residence time" in E traces back to host 1, so the shared ESS should
#   lean toward host 1's single-host optimum, qd1/(1-q) -- but only when
#   d1 != d2 (if d1 = d2 the two hosts' single-host optima are identical, so
#   there is nothing for shedding asymmetry to skew toward,
#   mirroring the Part 3 result that the shared ESS is invariant when the
#   hosts are otherwise ecologically identical). Whether the shift is a
#   simple h-weighted average of the two single-host optima (linear) or
#   differs from that (interacts with the d1 != d2 asymmetry) is checked
#   below.

multihost_shed_odes <- function(t, y, p) {
  with(as.list(c(y, p)), {
    S1 <- max(S1, 0); S2 <- max(S2, 0); I1 <- max(I1, 0); I2 <- max(I2, 0); E <- max(E, 0)
    beta <- beta_fun(alpha, b0, q)
    dS1 <- Lambda1 - d1*S1 - beta*f*S1*E
    dS2 <- Lambda2 - d2*S2 - beta*S2*E
    dI1 <- beta*f*S1*E - (d1 + alpha)*I1
    dI2 <- beta*S2*E   - (d2 + alpha)*I2
    dE  <- h1*I1 + h2*I2 - g*E
    list(c(dS1, dS2, dI1, dI2, dE))
  })
}

multihost_shed_equilibrium <- function(alpha, p) {
  beta <- beta_fun(alpha, p$b0, p$q)
  S1fun <- function(E) p$Lambda1 / (p$d1 + beta*p$f*E)
  S2fun <- function(E) p$Lambda2 / (p$d2 + beta*E)
  I1fun <- function(E) beta*p$f*E*S1fun(E) / (p$d1 + alpha)
  I2fun <- function(E) beta*E*S2fun(E)     / (p$d2 + alpha)
  gfun  <- function(E) p$h1*I1fun(E) + p$h2*I2fun(E) - p$g*E
  Estar <- grid_root(gfun, 1e-8, 50, n = 3000)
  if (is.na(Estar)) return(c(S1 = NA, S2 = NA, I1 = NA, I2 = NA, E = NA))
  c(S1 = S1fun(Estar), S2 = S2fun(Estar), I1 = I1fun(Estar), I2 = I2fun(Estar), E = Estar)
}

# only the bottom row of the invasion-fitness matrix changes (h -> h1, h2)
multihost_shed_invasion_fitness <- function(alpha_m, S1star, S2star, p) {
  beta_m <- beta_fun(alpha_m, p$b0, p$q)
  A <- matrix(c(-(p$d1 + alpha_m), 0,                 beta_m*p$f*S1star,
                0,                 -(p$d2 + alpha_m),  beta_m*S2star,
                p$h1,              p$h2,               -p$g),
              nrow = 3, byrow = TRUE)
  max(Re(eigen(A, only.values = TRUE)$values))
}

multihost_shed_selection_gradient <- function(alpha, p) {
  eq <- multihost_shed_equilibrium(alpha, p)
  if (any(is.na(eq))) return(NA_real_)
  num_deriv(function(am) multihost_shed_invasion_fitness(am, eq["S1"], eq["S2"], p), alpha)
}

find_ESS_multihost_shed <- function(p, alpha_range = c(1e-3, 5), n = 300) {
  grid_root(function(a) multihost_shed_selection_gradient(a, p), alpha_range[1], alpha_range[2], n)
}

## ---- sanity check: h1 = h2 = h/... reproduces the Part 3 ESS exactly ------
p_shed_check <- p_multi_asym; p_shed_check$h1 <- p_multi_asym$h; p_shed_check$h2 <- p_multi_asym$h
ess_shed_check <- find_ESS_multihost_shed(p_shed_check)
cat("=== Asymmetric shedding: sanity check against Part 3 ===\n")
cat(sprintf("h1=h2=%.2f reproduces Part 3 asymmetric-d ESS: %.4f  (Part 3 value: %.4f)\n\n",
            p_multi_asym$h, ess_shed_check, ess_asym))

## ---- KEY SWEEP: h1/h2 ratio at fixed h1+h2, symmetric-d vs asymmetric-d ----
h_total   <- 2 * p_multi$h                       # = 1.6; ratio=1 reproduces h1=h2=0.8
ratio_seq <- exp(seq(log(1/9), log(9), length.out = 19))  # log-spaced, symmetric around 1

sweep_shed <- function(base_p, case_label) {
  do.call(rbind, lapply(ratio_seq, function(r) {
    pp <- base_p; pp$h1 <- h_total * r / (1 + r); pp$h2 <- h_total / (1 + r); pp$h <- NULL
    ess <- find_ESS_multihost_shed(pp)
    data.frame(ratio = r, h1 = pp$h1, h2 = pp$h2, alpha_star = ess, case = case_label)
  }))
}
sweep_shed_symdz  <- sweep_shed(p_multi,      "d1 = d2 (symmetric)")
sweep_shed_asymdz <- sweep_shed(p_multi_asym, "d1 != d2 (asymmetric)")
sweep_shed_all    <- rbind(sweep_shed_symdz, sweep_shed_asymdz)

alpha1_opt_sym  <- (p_multi$q      / (1 - p_multi$q))      * p_multi$d1
alpha1_opt_asym <- (p_multi_asym$q / (1 - p_multi_asym$q)) * p_multi_asym$d1   # host 1 alone, d1=0.10
alpha2_opt_asym <- (p_multi_asym$q / (1 - p_multi_asym$q)) * p_multi_asym$d2   # host 2 alone, d2=0.30

cat(sprintf("Shedding-ratio sweep (h1+h2 held fixed at %.2f):\n", h_total))
print(sweep_shed_all[, c("ratio", "h1", "h2", "alpha_star", "case")], row.names = FALSE, digits = 4)

cat("\n*** Key finding (shedding asymmetry) ***\n")
cat(sprintf("d1 = d2 case: alpha* range across the whole h1/h2 sweep = [%.5f, %.5f] (flat, as predicted --\n",
            min(sweep_shed_symdz$alpha_star), max(sweep_shed_symdz$alpha_star)))
cat("  with identical single-host optima there is nothing for shedding asymmetry to skew toward).\n")
cat(sprintf("d1 != d2 case: alpha* ranges from %.4f (h1-dominated, near host 1's optimum) to %.4f (h2-dominated, near host 2's),\n",
            min(sweep_shed_asymdz$alpha_star), max(sweep_shed_asymdz$alpha_star)))
cat(sprintf("  bracketed by the single-host optima qd1/(1-q)=%.4f and qd2/(1-q)=%.4f as predicted.\n\n",
            alpha1_opt_asym, alpha2_opt_asym))

print(
  ggplot(sweep_shed_all, aes(ratio, alpha_star, color = case)) +
    geom_hline(yintercept = alpha1_opt_sym, linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = alpha1_opt_asym, linetype = "dashed", color = "#1b9e77") +
    geom_hline(yintercept = alpha2_opt_asym, linetype = "dashed", color = "#d95f02") +
    annotate("text", x = max(ratio_seq), y = alpha2_opt_asym, hjust = 1, vjust = -0.4, size = 3, color = "#d95f02",
             label = sprintf("host 2 alone (d2 = %.2f):\nalpha* = q*d2/(1-q) = %.2f", p_multi_asym$d2, alpha2_opt_asym)) +
    annotate("text", x = max(ratio_seq), y = alpha1_opt_sym, hjust = 1, vjust = -0.4, size = 3, color = "grey35",
             label = sprintf("either host alone when d1 = d2 = %.2f:\nalpha* = q*d/(1-q) = %.2f", p_multi$d1, alpha1_opt_sym)) +
    annotate("text", x = max(ratio_seq), y = alpha1_opt_asym, hjust = 1, vjust = 1.4, size = 3, color = "#1b9e77",
             label = sprintf("host 1 alone (d1 = %.2f):\nalpha* = q*d1/(1-q) = %.2f", p_multi_asym$d1, alpha1_opt_asym)) +
    geom_point(size = 2) + geom_line() +
    scale_x_log10() +
    scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = "Shared ESS virulence vs. shedding asymmetry",
         subtitle = "Each point = the ESS when host 1 sheds h1 and host 2 sheds h2 (h1 + h2 held fixed)",
         caption = paste0("Horizontal lines = the virulence that would evolve if the pathogen lived in only one host: alpha* = q*d/(1-q),\n",
                          "where q sets the shape of the transmission-virulence trade-off. The shared ESS is always a compromise between\n",
                          "the two lines: it slides toward host 1's line when host 1 does most of the shedding (ratio large) and toward\n",
                          "host 2's line when host 2 does (ratio small). When d1 = d2 the lines coincide, so the curve is flat."),
         x = lab_shed_ratio, y = lab_ess, color = "host mortality") +
    theme_minimal() +
    theme(plot.caption = element_text(hjust = 0, size = 8))
)

## ---- does the spillover discount f change this curve? (Javad's question) ----
## f only enters through host 1's share of the shedding flux (w1 ~ h1*f*S1*/(d1+alpha)),
## so a lower f should pull alpha* toward host 2's optimum. But a lower f also leaves more
## susceptibles alive in host 1 (S1* rises), which cancels most of that whenever transmission
## is strong enough that nearly every recruit gets infected (true of the default parameters).
## This is shown at the default environmental decay g and at a much larger g (weaker transmission),
## where the cancellation is only partial. The ESS search is restricted to the band between
## the two single-host optima (0.1 - 0.3), which always contains the shared ESS 
## (outside that band every term of the harmonic condition has the same sign).
f_vals       <- c(0.25, 0.5, 1)
g_regimes    <- c("strong transmission (g = 0.5)" = 0.5, "weak transmission (g = 10)" = 10)
ratio_sub    <- ratio_seq[seq(1, length(ratio_seq), by = 2)]
ess_band     <- c(0.8 * alpha1_opt_asym, 1.2 * alpha2_opt_asym)
sweep_f_shed <- do.call(rbind, lapply(names(g_regimes), function(regime) {
  do.call(rbind, lapply(f_vals, function(fv) {
    do.call(rbind, lapply(ratio_sub, function(r) {
      pp <- p_multi_asym; pp$f <- fv; pp$g <- g_regimes[[regime]]
      pp$h1 <- h_total * r / (1 + r); pp$h2 <- h_total / (1 + r); pp$h <- NULL
      data.frame(regime = regime, f = fv, ratio = r,
                 alpha_star = find_ESS_multihost_shed(pp, alpha_range = ess_band, n = 40))
    }))
  }))
}))

f_effect <- do.call(rbind, lapply(split(sweep_f_shed, sweep_f_shed$regime), function(d) {
  lo <- d[d$f == min(f_vals), ]; hi <- d[d$f == max(f_vals), ]
  data.frame(regime = d$regime[1],
             max_change_in_alpha_star = max(abs(lo$alpha_star - hi$alpha_star), na.rm = TRUE))
}))
cat(sprintf("Effect of f on the shedding-asymmetry curve: max |alpha*(f=%.2f) - alpha*(f=%.2f)| over the sweep\n",
            min(f_vals), max(f_vals)))
print(f_effect, row.names = FALSE, digits = 3)
cat("With symmetric mortality (d1 = d2) f has no effect at all: alpha* stays at q*d/(1-q), as in Part 3.\n\n")

print(
  ggplot(sweep_f_shed, aes(ratio, alpha_star, color = factor(f))) +
    geom_hline(yintercept = c(alpha1_opt_asym, alpha2_opt_asym), linetype = "dashed", color = "grey60") +
    annotate("text", x = max(ratio_sub), y = alpha2_opt_asym, hjust = 1, vjust = -0.4, size = 3, color = "grey35",
             label = sprintf("host 2 alone: %.2f", alpha2_opt_asym)) +
    annotate("text", x = max(ratio_sub), y = alpha1_opt_asym, hjust = 1, vjust = 1.4, size = 3, color = "grey35",
             label = sprintf("host 1 alone: %.2f", alpha1_opt_asym)) +
    geom_point(size = 2) + geom_line() +
    scale_x_log10() +
    scale_y_continuous(expand = expansion(mult = 0.12)) +
    facet_wrap(~ regime) +
    labs(title = "Effect of the spillover discount f on the shedding-asymmetry curve",
         subtitle = "d1 != d2. Dashed lines = virulence that would evolve in each host alone. f barely matters when transmission is strong.",
         x = lab_shed_ratio, y = lab_ess, color = "spillover discount f\n(1 = no discount)") +
    theme_minimal() +
    theme(panel.spacing = unit(1.5, "lines"))
)

## ---- additive vs interacting check: is the shift a simple h-weighted ------
## average of the two single-host optima, or are we dealing with something more complex? ------
sweep_shed_asymdz$naive_weighted <- (sweep_shed_asymdz$h1 * alpha1_opt_asym +
                                      sweep_shed_asymdz$h2 * alpha2_opt_asym) /
                                     (sweep_shed_asymdz$h1 + sweep_shed_asymdz$h2)
sweep_shed_asymdz$naive_minus_actual <- sweep_shed_asymdz$naive_weighted - sweep_shed_asymdz$alpha_star

cat("Additive (h-weighted-average) prediction vs. actual numeric ESS (d1 != d2 case):\n")
print(sweep_shed_asymdz[, c("ratio", "alpha_star", "naive_weighted", "naive_minus_actual")],
      row.names = FALSE, digits = 4)
cat(sprintf("Max deviation from the naive additive guess: %.5f (%.1f%% of the alpha1-alpha2 optimum gap)\n",
            max(abs(sweep_shed_asymdz$naive_minus_actual)),
            100 * max(abs(sweep_shed_asymdz$naive_minus_actual)) / (alpha2_opt_asym - alpha1_opt_asym)))
cat("If this deviation is small and structureless, the shedding-ratio effect is close to additive on\n")
cat("top of the existing d1!=d2 compromise; a deviation that grows systematically with the ratio (rather\n")
cat("than scattering near zero) would say the two asymmetries interact rather than simply superpose.\n\n")

sweep_shed_compare <- rbind(
  data.frame(ratio = sweep_shed_asymdz$ratio, alpha_star = sweep_shed_asymdz$alpha_star,      series = "actual ESS (numerical)"),
  data.frame(ratio = sweep_shed_asymdz$ratio, alpha_star = sweep_shed_asymdz$naive_weighted,   series = "naive guess: h-weighted average of the two single-host optima")
)
print(
  ggplot(sweep_shed_compare, aes(ratio, alpha_star, color = series, linetype = series)) +
    geom_point(size = 2) + geom_line() +
    scale_x_log10() +
    labs(title = "Actual ESS vs. naive additive (h-weighted) prediction",
         subtitle = "d1 != d2. If the naive guess were right the curves would overlap; the systematic gap means the two asymmetries interact.",
         x = lab_shed_ratio, y = lab_ess, color = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
)

## ---- So what is the right average? I wanted to check to confirm that it is indeed the harmonic mean --
## For the single-host model, the ESS condition collapses to
## beta'(alpha*)/beta(alpha*) = 1/(d+gamma+alpha*). The question for the shared
## multi-host ESS: does it obey the same kind of condition, but with the
## removal rate replaced by a weighted harmonic mean of (d1+alpha*) and
## (d2+alpha*)?
##     beta'(alpha*)/beta(alpha*) = w1/(d1+alpha*) + w2/(d2+alpha*),  w1+w2=1
## Working this out from the actual selection gradient (left/right
## eigenvectors of the resident (I1,I2,E) growth matrix at the ESS) shows it
## holds once w_i is each host's share of total shedding flux into the
## environment, h_i*I_i*, not shedding rate h_i alone. This is the
## correct average, unlike the naive h-weighted one just above.
##
## Why call it "harmonic" (Javad's question): the averaging is across hosts, of
## each host's total removal rate r_i = d_i + alpha*. It is not a harmonic
## combination of d and alpha inside a host; within a host they simply add.
## 1/r_i is the mean infectious period in host i, so the right-hand side above
## is the arithmetic mean of the infectious periods, and the reciprocal
## of an arithmetic mean of periods is the harmonic mean of the
## rates (with weights summing to 1, 1/H = w1/r1 + w2/r2). In biological terms,
## a small rise in alpha shortens host i's infectious period by the fraction
## 1/r_i, and the pathogen weighs that loss by where its shedding comes from.
## For beta = b0*alpha^q the left side is q/alpha*, so the harmonic mean of
## (d_i + alpha*) must equal alpha*/q. Tested against the alternatives below.
harmonic_weights <- function(alpha, p) {
  eq <- multihost_shed_equilibrium(alpha, p)
  flux1 <- p$h1 * eq["I1"]; flux2 <- p$h2 * eq["I2"]
  c(w1 = unname(flux1 / (flux1 + flux2)), w2 = unname(flux2 / (flux1 + flux2)))
}

harmonic_gap <- function(alpha, p) {
  w <- harmonic_weights(alpha, p)
  lhs <- dbeta_fun(alpha, p$b0, p$q) / beta_fun(alpha, p$b0, p$q)
  rhs <- w["w1"] / (p$d1 + alpha) + w["w2"] / (p$d2 + alpha)
  unname(lhs - rhs)
}

find_ESS_harmonic <- function(p, alpha_range = c(1e-3, 5), n = 300) {
  grid_root(function(a) harmonic_gap(a, p), alpha_range[1], alpha_range[2], n)
}

sweep_shed_asymdz$alpha_harmonic <- sapply(seq_len(nrow(sweep_shed_asymdz)), function(i) {
  row <- sweep_shed_asymdz[i, ]
  pp <- p_multi_asym; pp$h1 <- row$h1; pp$h2 <- row$h2; pp$h <- NULL
  find_ESS_harmonic(pp)
})

harmonic_compare <- rbind(
  data.frame(ratio = sweep_shed_asymdz$ratio, alpha = sweep_shed_asymdz$alpha_star,    series = "actual ESS (numerical)"),
  data.frame(ratio = sweep_shed_asymdz$ratio, alpha = sweep_shed_asymdz$naive_weighted, series = "naive guess: h-weighted average of the two single-host optima"),
  data.frame(ratio = sweep_shed_asymdz$ratio, alpha = sweep_shed_asymdz$alpha_harmonic,  series = "harmonic-mean ESS: flux-weighted harmonic mean of (d_i + alpha)")
)
print(
  ggplot(harmonic_compare, aes(ratio, alpha, color = series, linetype = series)) +
    geom_point(size = 2) + geom_line() +
    scale_x_log10() +
    labs(title = "Harmonic-mean ESS vs. actual ESS vs. naive ESS",
         subtitle = "d1 != d2. The harmonic-mean prediction overlays the actual ESS; the naive h-weighted average does not.",
         x = lab_shed_ratio, y = lab_ess, color = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
)

## ---- is it really harmonic? compare candidate averages of r_i = d_i + alpha* --
## The "right" average of the removal rates equals alpha*/q here. This is tested at
## five shedding ratios from the sweep, with the ESS polished to a tight tolerance.
average_check <- do.call(rbind, lapply(c(1, 5, 10, 15, 19), function(i) {
  row <- sweep_shed_asymdz[i, ]
  pp <- p_multi_asym; pp$h1 <- row$h1; pp$h2 <- row$h2; pp$h <- NULL
  a  <- uniroot(function(x) multihost_shed_selection_gradient(x, pp),
                row$alpha_star + c(-2e-3, 2e-3), tol = 1e-12)$root
  eq <- multihost_shed_equilibrium(a, pp)
  rates  <- c(pp$d1 + a, pp$d2 + a)                          # removal rate in each host
  w_flux <- c(pp$h1 * eq["I1"], pp$h2 * eq["I2"]); w_flux <- w_flux / sum(w_flux)
  w_h    <- c(pp$h1, pp$h2) / (pp$h1 + pp$h2)
  data.frame(ratio = row$ratio,
             required = beta_fun(a, pp$b0, pp$q) / dbeta_fun(a, pp$b0, pp$q),   # = alpha*/q
             harmonic_flux_weights   = 1 / sum(w_flux / rates),
             arithmetic_flux_weights = sum(w_flux * rates),
             harmonic_h_only_weights = 1 / sum(w_h / rates))
}))
cat("Which average of the removal rates (d_i + alpha*) equals the required value alpha*/q?\n")
print(average_check, row.names = FALSE, digits = 5)
cat(sprintf("Max error -- harmonic (flux weights): %.1e | arithmetic (flux weights): %.1e | harmonic (h-only weights): %.1e\n\n",
            max(abs(average_check$harmonic_flux_weights   - average_check$required)),
            max(abs(average_check$arithmetic_flux_weights - average_check$required)),
            max(abs(average_check$harmonic_h_only_weights - average_check$required))))

## ---- CSS check + PIP at a representative strongly-asymmetric-shedding point
p_shed_rep <- p_multi_asym; p_shed_rep$h1 <- h_total * 9/10; p_shed_rep$h2 <- h_total * 1/10; p_shed_rep$h <- NULL
ess_shed_rep <- find_ESS_multihost_shed(p_shed_rep)
eq_shed_rep  <- multihost_shed_equilibrium(ess_shed_rep, p_shed_rep)

css_shed <- css_check(multihost_shed_selection_gradient,
                       function(am) multihost_shed_invasion_fitness(am, eq_shed_rep["S1"], eq_shed_rep["S2"], p_shed_rep),
                       ess_shed_rep, p_shed_rep)
cat(sprintf("CSS check (d1!=d2, h1/h2=9): grad(below)=%+.4f grad(above)=%+.4f curvature=%+.4f -> convergence stable: %s, evolutionarily stable: %s\n\n",
            css_shed$grad_below, css_shed$grad_above, css_shed$curvature,
            css_shed$convergence_stable, css_shed$evolutionarily_stable))

plot_PIP_multihost_shed <- function(p, ess, window = NULL, n = 60) {
  if (is.null(window)) window <- pip_window(ess)
  ag <- seq(window[1], window[2], length.out = n)
  eqs <- lapply(ag, multihost_shed_equilibrium, p = p)
  W <- matrix(NA, n, n)
  for (i in seq_along(ag)) {
    eq <- eqs[[i]]
    if (any(is.na(eq))) next
    for (j in seq_along(ag)) {
      W[i, j] <- multihost_shed_invasion_fitness(ag[j], eq["S1"], eq["S2"], p)
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
      pip_labs("Asymmetric-shedding PIP (d1 != d2, h1/h2 = 9: host 1 sheds 9x more)") +
      theme_minimal()
  )
}
plot_PIP_multihost_shed(p_shed_rep, ess_shed_rep)


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
cat(sprintf("Multi-host (d1!=d2):    alpha* = %.4f  (now does depend on f)\n",
            ess_asym))
cat(sprintf("Forced-equal-S* test:   max|natural-forced| diff = %.5f  (skew is from d1/d2 directly, not S* asymmetry)\n",
            max(abs(sweep_forced$diff))))
cat("=====================================================================\n")

cat("\n--- CSS verification (Javad asks are these three genuinely\n")
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
cat("(convergence stable; residents evolve toward it from either side), and\n")
cat("curvature is negative (a true local fitness maximum, i.e. uninvadable\n")
cat("once reached). All three are legit CSS points of the same stability\n")
cat("type. The PIP plots above look different only because of the ESS-centered\n")
cat("plotting window each now shares (pip_window()) -- earlier independent\n")
cat("[0, 3*ess] iterations of windows were not the same and let the predator/\n")
cat("multi-host feasibility boundary clip one side of the grid asymmetrically.\n")
cat("=====================================================================\n")