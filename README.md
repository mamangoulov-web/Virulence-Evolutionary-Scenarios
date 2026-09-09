# Virulence Evolution Across Ecological Scenarios

This is a project from a research collaboration I've been doing at SFU with Javad, a PhD student, and his supervisor Dr. Ben Ashby, in the field of host-parasite coevolution and disease dynamics. I'm not enrolled there; I'm a graduate from the University of Toronto, but Dr. Ailene MacPherson, also of SFU, introduced me to the group based on my background, and I've been working through some virulence-evolution models with Javad ever since, first by hand and then checking everything numerically here. This is my first project, so there's certainly room for improvement, but I think it represents a solid collection of my findings thus far.

My background is pretty unique: I specialized in biology and minored in chemistry at UofT — so I've got a strong analytical background and understanding of the reasons behind all of this — but took a liking for epidemiology and computational biology towards the back half of my fourth year in Toronto. After venturing into the veterinary field, I realised that my passion for the industry wasn't as keen as I had hoped, and pivoted to a more statistical, theoretical world, where I've so far greatly enjoyed my career. I hope to begin a Master's or PhD in comp-bio or epidemiology in 2028, and this research will help kick-start my journey into this field.

## Background

Everything started as a basic derivation: take a baseline SIS model with demography, solve for the endemic equilibrium, write down the invasion fitness of a rare mutant strain, and find the ESS (evolutionarily stable) virulence via the standard adaptive-dynamics selection-gradient method. I derived all of that by hand first, then coded it up (Part 1 below) to check the numbers against the algebra.

After that, we sketched two extensions on a whiteboard which became the foundation of everything below:

- what happens to ESS virulence when there's a predator hunting both susceptible and infected hosts (Part 2)
- what happens when the pathogen has to jump between two different host species through a shared environmental reservoir — think a spillover host and a wildlife reservoir (Part 3)

Later I added a third extension on my own (Part 5): letting the two hosts shed the pathogen into the environment at different rates instead of a shared rate, to see whether that skews the shared ESS toward whichever host is shedding more.

None of this is published or "done"; it's an active back-and-forth with my collaborators, so treat the open questions section below as truly open, not rhetorical.

## What's in the script

Everything lives in one file, [`virulence_evolution_scenarios.R`](virulence_evolution_scenarios.R), split into five parts:

1. **Baseline SIS model** — no predator, one host. Verifies the by-hand derivation that ESS virulence alpha\* = gamma + d, for a square-root transmission-virulence trade-off.
2. **Predator scenario** — adds a predator compartment that hunts both S and I hosts. Predation acts structurally like extra background mortality, so ESS virulence should rise with predation pressure. Also investigates the "healthy herds" concept — does it change if predators preferentially cull sick hosts?
3. **Multi-host / shared environmental transmission** — a spillover host and a reservoir host both feed into (and get infected from) a shared environmental pool. Sweeps how the shared ESS depends on transmission efficiency to the spillover host, and on whether the two hosts share the same background mortality.
4. **Summary + CSS check** — prints a wrap-up comparing all three scenarios, and formally verifies each ESS is a continuously stable strategy, not just a root of the selection gradient.
5. **Asymmetric shedding** — extends Part 3 so each host sheds into the environment at its own rate rather than a shared one. Tests whether the resulting ESS shift is just a weighted average of the two single-host optima, or something more complicated.

For every scenario, the script:

- codes up the actual ODEs so you can simulate them and watch the trajectories directly
- solves for the resident ecological equilibrium as a function of virulence, algebraically where possible and by root-finding where it isn't
- gets the invasion fitness of a rare mutant and the resulting selection gradient, then finds the ESS
- draws a pairwise invasibility plot (PIP) confirming the ESS is convergence- and evolutionarily stable

## Key findings so far

- Predation pressure pushes ESS virulence up, which matches the prediction as predation is structurally just extra mortality.
- In the shared-environment scenario, if the two hosts have identical background mortality, the shared ESS is completely invariant to transmission efficiency, shedding rate, and environmental decay — it only starts to matter once the hosts are ecologically different.
- When the hosts do differ, the shared ESS becomes a compromise between their two single-host optima, skewed toward whichever host contributes more to the environmental reservoir (through unequal exposure, or unequal shedding). That skew is not a simple weighted average; there is an interaction between the two asymmetries (by up to ~27% of the gap between the single-host optima).

## Open questions

- Should virulence actually be allowed to evolve separately per host, instead of being forced into one shared trait? Every multi-host scenario here assumes a single shared alpha.

## Running it

You'll need R with the `deSolve` and `ggplot2` packages (the script installs them automatically if they're missing). Then:

```r
source("virulence_evolution_scenarios.R")
```

Bear with me here; the code does take a few minutes to run so don't assume anything is broken.

## License

[MIT](LICENSE)
