# Linear Algebra Honors Project: Rocketry

Math 7 (Linear Algebra) Honors Project, Las Positas College, Fall 2026
Student: Talha Ahmed Shaik · Mentor: Prof. Bhairav Singh

## Overview

This project recovers the transonic drag curve, Cd vs. Mach (M 0.8 to 1.2), of the Zenith 2 high-power rocket using least squares. It then uses numerical linear algebra to measure how trustworthy that recovered curve is.

Coast-phase flight data produces far more equations than unknowns, so the problem is an overdetermined system Ax = b. Least squares gives a best fit; the condition number of A says how much noise in the data is amplified in the answer. The central question is how the shape of the trajectory changes the condition number of the design matrix.

## Parts

1. **Theory:** matrix norms, condition number κ(A) = σ_max / σ_min, perturbation bounds, and why κ(AᵀA) = κ(A)² favors QR/SVD over the normal equations.
2. **Application:** build the least-squares problem from synthetic coast-phase data, compute κ for different trajectory shapes, apply weighted least squares, estimate uncertainty on Cd(M), and compare against RasAero II.

## Data

Zenith 2 has not flown yet, so synthetic trajectories are the primary dataset. Real flight telemetry may be added later if the rocket flies.

Workflow: known Cd(M) → simulate trajectory → sample and add sensor noise → recover Cd(M) → compare recovered vs. known.

## Repository layout

| Folder | Contents |
|---|---|
| `sim/` | Synthetic trajectory generation |
| `analysis/` | Least-squares drag extraction and conditioning analysis |
| `data/` | Generated datasets and reference exports (e.g. RasAero II CSVs) |
| `docs/` | Report, poster, notes, and references |

## Timeline

| Date | Milestone |
|---|---|
| Sept 11 | Proposal submitted |
| Oct 2 | Literature search; least-squares problem built and tested on simulated data |
| Oct 23 | Conditioning analysis and baseline drag recovery on synthetic data |
| Nov 13 | Weighted least squares and uncertainty quantification |
| Nov 25 | Rough draft to Prof. Singh |
| Dec 4 | Final report, code, and poster |

## Tools

MATLAB
