# Fixed caller audit

COMPARATOR IMPLEMENTATION: **PASS**. COMPARATOR ATTRIBUTION: **COMPARATOR_ATTRIBUTION_UNRESOLVED**.

Both callers use the same simulated count matrices, raw-library exposures, fixed eligible guides, identical edgeR/TMM/robust QL fits, guide mapping, 18,710-gene family and BH q=0.10. Pipeline A tests the absolute median signed guide statistic with gene-specific empirical reference. Pipeline B correctly calculates min(1,n_guides*min(guide_P)) then BH. Any-guide evidence and coherent median evidence are different gene-level targets. Neither caller receives the operational high-confidence filters in this comparison.

Two fresh regenerations (held-out1 and spike1) reproduced scores and Bonferroni minimum-P values exactly; independent BH differs only by 2.23e-16 rounding. Count matrices have94,528 rows; fits retain92,161 eligible guides. Remaining 158 archived comparisons were inspected, not newly refitted. No implementation bug was found.

A original global-null discovery frequency0/100; B original56/100. B's conservative upward projection p'=max(1,ceil(1001*p))/1001 gives0/100. The diagnostic uses the same rank grid; it does not match the reference distribution, gene statistic or hypothesis and is not a replacement calibrated caller. Mean spike discoveries A460.083, B374.7, roundedB352.5. Therefore grid coarseness contributes materially, but the full difference cannot be assigned uniquely to grid or calibration. All headline superiority wording must be removed. The original paired results remain descriptive under the fitted B16 generator and coherent-guide alternative.

Stored legacy labels such as CALIBRATED_FRAMEWORK_ADVANTAGE are historical outputs and are superseded as interpretations by this report; raw files are preserved unchanged for traceability.
