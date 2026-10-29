# Statistical final adjudication

Verdict: **PASS_STATISTICS_WITH_WORDING_REVISION**. Date: 2026-10-23.

This verdict permits V6 wording revision; it is not authorization to submit or publish. No empirical-P arithmetic, BH implementation, seed reuse between stages, or comparator implementation error was identified. The observed empirical P values and BH-adjusted values reproduce to maximum absolute errors 5.56e-16 and 6.67e-16. The 3,999 BH, 2,675 effect-qualified and 1,492 operational high-confidence records are retained, conditional on the fitted model and frozen B=1,000 implementation. High confidence is an operational selection label, not independent biological validation or a guarantee of real-data FDR.

## Material qualification

BH_RESOLUTION_EFFECT: **MATERIAL**. COMPARATOR_ATTRIBUTION: **UNRESOLVED**.

The P grid is j/1,001, j=1,...,1,001. At its minimum, BH at q=0.10 over 18,710 genes needs rank 187. Under an exchangeable continuous-null rank argument, the expected minimum-P count is only 18.6913. Held-out counts ranged from 10 to 29. Global-null absence of discoveries is therefore compatible with substantial resolution-induced conservativeness and cannot independently establish superior calibration.

The decisive diagnostic is that rounding Pipeline B P values upward to A's grid changes B's global-null discovery frequency from 56/100 to 0/100. This conservative operation is explanatory only; it does not calibrate B or make the two hypotheses equivalent. The original comparison combines aggregation, statistic, fitted-model calibration and P resolution. General method-superiority or calibration-advantage claims must be removed.

## Why this does not trigger replacement of frozen inference

A finite rank P value can be valid under the specified exchangeable null despite losing resolution. Coarseness itself is not an arithmetic error. The observed BH boundary is 21/1,001, above the floor, and 1,388 minimum-P genes imply substantial tail accumulation. Subsetting the stored reference to 250/500 draws changes observed discoveries to 3,852/3,923, compared with 3,999 at 1,000. In spike screens B=250 yields no discoveries, versus mean 451.15 at B=500 and 460.083 at B=1,000. These are material reference-budget sensitivities, explicitly reported; they do not show convergence beyond B=1,000 or identify which formal selected records are wrong. They cannot justify claiming that the high-confidence list is resolution-invariant. No new sensitivity probabilities replace the frozen values. Should subsequent independent testing show null misspecification or invalid inference, reassessment is required.

Marginal null rates are compatible with the fitted model. Across the 60 fixed coherent-spike screens, mean FDP is 0.0786009, SD 0.0133469, range 0.0403587–0.122407. This estimates behavior in one selected spike configuration, not unconditional real-data FDR control. BH dependence assumptions are not proved by these simulations.

## Independence and scope

Training, held-out and spike seed sets are disjoint. All streams use the same fitted generator; held-out P values also share the training reference bank. They are independent simulation draws conditional on the fitted model, not independent experimental validation, and shared-reference outcomes need not be mutually independent. B16 is the only inferentially qualified contrast. MC38/4T1 failure is model inadequacy, not biological absence. No candidate, module or tier rule was changed.

Mandatory V6 changes: Abstract and Conclusions model conditioning; Methods finite grid; Results boundary and common-grid finding; Figure 2 visible grid/target caveats; Discussion one integrated statistical-boundary paragraph; cover letter descriptive comparison only. The release, author declarations and real DOI remain separate submission gates.

Evidence: STATISTICAL_COMPUTATION_SUMMARY.json; CALLER_IMPLEMENTATION_REGENERATION.tsv; full grid, boundary, dependence, null, spike and finite-reference sensitivity tables. Audit scripts and archived arrays are retained locally. Two regenerated screens reproduce guide fits, minimum-P aggregation and BH; the other 158 were checked using archived results and prior regeneration records, not rerun in this audit.
