# Finite empirical P and BH

For m=18,710, q=0.10, B=1,000, p_j=j/1,001. A possible BH crossing at bin j requires cumulative rank k≥ceil(187100*j/1001). Only j≤100 can be rejected. The first bin requires k=187, not 188. All 1,001 bins and feasibility flags are in P_RESOLUTION_BH_ANALYTICAL_AUDIT.tsv.

This is a cumulative-bin condition: 187 minimum-P genes are sufficient for a first-bin crossing, but are not necessary for a crossing at a larger P bin. A global-null screen could cross at any feasible bin. No overall BH probability follows solely from the minimum-bin calculation.

Under exchangeable continuous null ranks, E(N_min)=18710/1001=18.69130869. If gene ranks additionally were independent discrete uniform variables, Pr(N_min≥187)=1.73654e-116. This is an illustrative first-bin binomial probability, not the fitted-model probability of any BH discovery and not a claim that genes are independent. Non-continuous scores or ties can be more conservative.

Observed: 1,388/18,710 genes (7.41849%) occupy the floor; all 1,388 are among the 3,999 BH records. The largest rejected P is 21/1,001=0.02097902098 at rank 3,999. Tail score differences within the minimum-P bin cannot be resolved by P ranking. Consult OBSERVED_P_GRID_DIAGNOSTICS.tsv for HC and effect-qualified counts at every bin and score ranges; BH_REJECTION_BOUNDARY.tsv retains nearby tied records.

The actual global-null minimum-bin count ranges 10–29. With upward common-grid rounding, B also gives 0/100 screens with a discovery. Thus resolution materially reduces the opportunity for sparse extreme-tail calls. The exact portion attributable to resolution, aggregation and calibration is unresolved; these diagnostics do not establish a universal probability of 0/100.

The bounded first-B analysis used archived draws only. Observed BH counts: B250=3,852, B500=3,923, B1000=3,999. Spike discoveries: B250 all zero, B500 mean451.15, B1000 mean460.083. No additional simulation budget was selected and no formal P/FDR was replaced. Larger-B convergence remains untested.
