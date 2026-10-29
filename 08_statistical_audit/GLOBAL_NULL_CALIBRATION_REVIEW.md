# Global-null and mixed-null review

Global-null calibration: **SUPPORTED_WITHIN_MODEL**. Mixed-null FDR interpretation: **SUPPORTED_WITHIN_TESTED_SPIKES**, with narrow generalizability.

Marginal calibration is distinct from BH calling behavior. Across 100 screens, mean P≤0.05 rate=0.0501378942, between-screen SD=0.0016181933 and descriptive MCSE=0.0001618193; mean P≤0.01 rate=0.0099877071, SD=0.0007796002, MCSE=0.0000779600. Finite-grid reference rates are 50/1001=0.04995005 and 10/1001=0.00999001. Pooled median=0.500; discrete-reference KS distance≈0.000572. Pooled genes are not independent experimental replicates; the KS value is descriptive, not an independent-gene test.

All 100 screens have zero BH discoveries. This agrees with sparse tail counts on the coarse grid and is not proof of unconditional FDR control. Each held-out screen has independently generated counts conditional on the fitted parameters; reference-bank reuse and common parameter estimation remain. Between-screen MCSE describes the fixed-reference simulation exercise, not all uncertainty in model fitting or training-bank estimation.

A fixed, evenly spaced subset of 1,000 alphabetical genes has pairwise signed-score correlations with mean -0.0001068 and 1st/99th percentiles -0.23165/0.23221 across 100 draws. This limited noisy diagnostic does not establish independence, PRDS or a general BH theorem. Generator and normalization dependencies remain.

Sixty coherent-spike screens use one fixed stratified 600-gene set, 200 genes per absolute log2 effect 0.5/1/2 and balanced signs. Mean BH powers are 0.255833/0.866333/0.997. Mean realized FDP=0.07860093 (SD0.01334687; range0.04035874–0.12240664), mean false calls36.25 and mean total calls460.0833. The 18,110 nominally unspiked genes have mean P≤0.05 rate0.05012976 and P≤0.01 rate0.01007086. Composition rescaling alters their relative abundance; this is examined in the archived SPIKE_COMPOSITION table, not assumed absent. Target prevalence is 600/18710=3.20684%, selected on baseline abundance and guide count. Coherent guide effects do not cover guide inefficiency, off-target effects, other prevalence or alternative generators. No broad real-data FDR claim is supported.


Archived spike composition check: comparator exposure-normalized shifts were approximately 1e-14 log2; immune columns shifted by +0.01908724408218 log2. Thus unspiked raw counts were not intentionally spiked, but immune library composition was changed; nominal unspiked-gene behavior is empirical evidence within this construction, not absence of composition effects.
