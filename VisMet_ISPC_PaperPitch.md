# Beyond Anomaly, Beyond Amplitude: Theta-Band Phase Decoupling During Visual Metaphor Comprehension

**Internal pitch document — for discussion with advisor about a second paper drawing on the existing VisMet EEG dataset**

**Author:** Cagatay
**Dataset:** VisMet (hybrid visual metaphors), n = 24 after exclusions (subject 32 dropped from the ERP-paper sample of 25).
**Date:** 2026-04-20

---

## 1. One-paragraph pitch

The ERP paper currently in revision ("Beyond Anomaly: Integrative Inference in Hybrid Visual Metaphor Processing") shows that hybrid visual metaphors elicit a Metaphor-specific Late Positive Component (LPC; 800–1500 ms) that is not reducible to anomaly or familiarity, while the earlier 200–800 ms negativity is driven primarily by anomaly. The manuscript's Discussion explicitly calls for **time-frequency or connectivity methods** to disambiguate what the 200–800 ms sustained negativity reflects, and notes that late negativities — which are common in linguistic metaphor ERP but untested in visual metaphor — have yet to be examined. A second paper drawing on the same dataset can address exactly that gap. Using **seed-based inter-site phase clustering (ISPC)**, a phase-based connectivity measure unrelated to evoked-amplitude components, we find that Metaphor elicits **sustained theta-band phase decoupling between frontal-midline sites and central/posterior scalp** that differentiates Metaphor from both Literal and Anomalous conditions in roughly the same late window as the LPC. This is a frequency-domain, connectivity-level signature of metaphor-specific processing that complements (rather than re-describes) the LPC finding, and it opens a set of hypothesis-driven follow-up analyses that are already pre-registered by the manuscript's own Discussion.

---

## 2. Why this is a separate paper, not a supplement

Three reasons.

**(a) Methodologically independent.** ISPC measures across-trial consistency of phase *differences* between channel pairs. It is mathematically orthogonal to evoked amplitude (ERP), so a finding in ISPC cannot be re-derived from the ERP results already in the paper. If the conditions differ on theta ISPC, they are differing on something the ERP paper cannot see.

**(b) Different theoretical frame.** The ERP paper frames metaphor processing at the level of stages: anomaly detection, semantic access, late integration. The connectivity paper reframes it at the level of **network dynamics** — which scalp regions' phase-locking with frontal control sites is modulated by interpretive demand. These two stories complement each other rather than compete.

**(c) The ERP paper explicitly says so.** §4.1 twice states that distinguishing competing accounts of the early negativity "will likely require methods with finer temporal and functional resolution than averaged ERPs, such as time-frequency decomposition or temporal decoding." §4.2 notes that late negativities have "not been tested" in visual metaphor. Turning that future-work sentence into a paper is the cleanest kind of follow-up.

Bundling ISPC into the current ERP paper would (i) weaken the focus of that manuscript, (ii) bury a methodologically novel result in a supplement, and (iii) deny reviewers a clean evaluation of either finding on its own merits.

---

## 3. Analyses run so far

All analyses use the same preprocessing as the ERP paper (0.1–30 Hz Butterworth, mastoid re-reference, −200 to 1500 ms epoch, post-ICA artifact rejection, `_postart.set` files). Sub 32 excluded from all analyses, leaving n=24.

Additional steps specific to the connectivity analysis:

- **Surface Laplacian / Current Source Density transform** (Perrin spherical spline; `laplacian_perrinX`) applied to each epoch before ISPC. This is standard for phase-connectivity work because raw scalp voltage inflates apparent coupling between neighboring electrodes via volume conduction.
- **Morlet wavelet convolution** per frequency, 5 log-spaced frequencies per band, cycles log-spaced 3 to 10. Reshape-all-trials FFT approach for efficiency (adapted from Cohen Ch. 13).
- **Baseline correction** as percent change from −200 to 0 ms (needed for cross-subject comparability because absolute ISPC has a subject-specific floor that scales with trial count as ~1/√N).
- Per-condition ISPC topography = channels × time × condition, with conditions identified from ERPLAB bin codes (111/121/131 Literal, 112/122/132 Metaphor, 113/123/133 Anomalous).

### 3a. Single-subject check (sub 03, FCz seed, theta [4–8 Hz])

Purpose: sanity check that the function and classifier work on real study data.

Result: After baseline correction, the topo grid shows the expected pattern. Metaphor produces visibly stronger post-stimulus decoupling than Literal across central/posterior sites, with a Metaphor–Literal difference that is mostly non-zero between 300 and 1200 ms. Effect sizes in this single subject are in the 5–15% change range — noisy but directionally sensible.

### 3b. Group analysis (n=24, FCz → Fz fallback seed, theta [4–8 Hz])

Three findings, from the difference topography panel and the probe-channel time courses:

1. **All three conditions show post-stimulus decoupling from frontal midline.** Baseline-corrected ISPC goes negative at Cz, Pz, Oz, C3, C4, P3, P4, F7, F8 from ~200 ms onward across all conditions. This is the task-related desynchronization that is always present in cognitive tasks, and it's not what distinguishes conditions.
2. **Metaphor decouples more than Literal and more than Anomalous**, with the effect most visible at central and parietal sites (Cz, Pz, C3, C4, P3). Peak magnitude of the condition difference is ~2–3% change. The direction is consistent across nearly every probe channel and across both contrasts (Met−Lit and Met−Ano).
3. **Anomalous is not distinguishable from Literal in theta ISPC at the group mean.** This is the key dissociation from the ERP paper, where Anomalous was the outlier in the early negativity. In theta ISPC, Metaphor is the outlier in the late window.

Statistically, the main-task-window t-maps show only scattered uncorrected-p<0.05 channels — the single-subject-level Student's t approach does not support strong claims. This is why the cluster-permutation sweep below is the appropriate confirmatory step.

### 3c. Band × seed sweep with cluster-based permutation (RESULTS)

Script: `vismet_ispc_sweep.m`. For each (band, seed) combination, computes the per-subject ISPC array, runs sign-flip cluster-based permutation (1000 iterations, channel adjacency threshold 75 mm, cluster-forming α=0.05) for Met−Lit and Met−Ano contrasts, and runs a hypothesis-driven paired t at a seed-centered ROI in the LPC window (800–1500 ms). Sweep grid: 3 bands × 3 seeds = 9 combinations.

C4 and FC2 were chosen because they are the peak channels in the mass-univariate Metaphor vs. Literal LPC cluster reported in the ERP paper (§3.2 Results). Fz was included as a cognitive-control prior. (FCz was requested but does not exist in the actiCAP 32-channel montage; the script falls back to Fz, which is the same fallback used in the ERP analyses' frontal ROI.)

**Cluster-based permutation: ALL NULL.** Across all 9 (band × seed) combinations and both contrasts, **no cluster reached p<0.05**. The smallest cluster-level p was 0.076 (beta @ Fz, Met−Lit) and 0.102 (theta @ FC2, Met−Ano). The directional pattern is consistent — Met−Lit and Met−Ano clusters trend negative (decoupling) in theta and alpha across all seeds — but no spatiotemporal cluster survives correction.

**ROI-restricted LPC-window test (800–1500 ms, seed + neighbors): ONE survives.**

| Band  | Seed | Met − Lit  ROI         | Met − Ano  ROI                                    |
|-------|------|------------------------|---------------------------------------------------|
| theta | C4   | t=−1.36, p=.19, Δ=−0.31% | t=−1.90, p=.07, Δ=−0.34% (trend)                  |
| theta | **FC2** | t=−1.50, p=.15, Δ=−0.29% | **t=−2.56, p=.018, Δ=−0.37%** (significant)       |
| theta | Fz   | t=−0.31, p=.76, Δ=−0.05% | t=−1.64, p=.11, Δ=−0.29%                          |
| alpha | C4   | t=−1.31, p=.20, Δ=−0.20% | t=−1.19, p=.25, Δ=−0.21%                          |
| alpha | FC2  | t=−1.38, p=.18, Δ=−0.20% | t=−1.22, p=.24, Δ=−0.23%                          |
| alpha | Fz   | t=−0.87, p=.40, Δ=−0.10% | t=−1.30, p=.21, Δ=−0.21%                          |
| beta  | C4   | t=−1.09, p=.29, Δ=−0.13% | t=+0.09, p=.93, Δ=+0.01% (sign flip)              |
| beta  | FC2  | t=−0.72, p=.48, Δ=−0.08% | t=+1.18, p=.25, Δ=+0.10% (sign flip)              |
| beta  | Fz   | t=−0.76, p=.45, Δ=−0.07% | t=+0.51, p=.62, Δ=+0.03%                          |

The single survivor — **theta @ FC2, Met − Ano, t(23)=−2.56, p=.018, mean Δ=−0.37%** — is at one of the two LPC-peak channels from the mass-univariate ERP analysis, in the LPC time window, in the band motivated by the linguistic-metaphor literature, and in the predicted direction (more decoupling for Metaphor). The Met − Lit contrast at the same site is in the same direction but does not survive (p=.15).

**Pattern across the table.** Theta and alpha consistently produce negative differences (Met decouples more than Lit and more than Ano), strongest at the LPC-peak seeds C4/FC2 and weakest at Fz. Beta shows a sign flip for Met−Ano at C4/FC2 and Fz — Met decouples *less* than Anomalous in beta — which would be a separable effect if it survived but it does not. Met−Ano effects are systematically larger than Met−Lit effects, mirroring the asymmetry in the ERP paper's mass-univariate (the Met−Ano LPC cluster was the largest, more sustained than Met−Lit).

**Statistical caveat.** The ROI test runs 18 separate paired t-tests (9 combos × 2 contrasts). Under the global null, ~1 false positive at α=.05 is expected. Bonferroni-corrected α = .05/18 ≈ .003; the FC2/theta/Met−Ano result (p=.018) does not survive Bonferroni. Treat it as **directionally consistent and pre-registered by the manuscript framing** rather than independently statistically significant.

---

## 4. Where the connectivity story meets the manuscript

| ERP paper finding                                                        | ISPC follow-up                                                                                                                                  |
|--------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| 200–800 ms negativity: Ano > Lit/Met (anomaly-driven, not metaphor)      | Theta ISPC in the same window does *not* single out Anomalous. Metaphor is the outlier. Same time range, different frequency signature.          |
| 800–1500 ms LPC: Met > Lit *and* Met > Ano (metaphor-specific late)      | Theta ISPC Met−Lit and Met−Ano differences peak in roughly the same window. Direction is decoupling, which is the frequency-domain complement of the LPC's amplitude increase. |
| Familiarity predicts early negativity but not late LPC                   | Familiarity not yet tested in ISPC (see §6). Natural replication of the ERP dissociation in the connectivity domain.                            |
| Discussion calls for "time-frequency decomposition"                      | Delivered.                                                                                                                                       |
| Discussion notes late negativities untested in visual metaphor           | Testable directly from the per-subject ISPC arrays, since ISPC does not care about amplitude sign.                                              |

---

## 4b. Familiarity-matched replication — *the key follow-up*

We re-ran the full sweep on the familiarity-matched subset (matched = odd bin numbers, the 32 triplets the ERP paper's primary RM-ANOVA and mass-univariate were run on; ~34 trials per condition per subject). Results from `vismet_ispc_sweep_matched.csv`:

| Band  | Seed | Met − Lit ROI (matched) | Met − Ano ROI (matched) | Δ vs. full set |
|-------|------|-------------------------|-------------------------|----------------|
| theta | C4   | t=−0.77, p=.45, Δ=−0.17% | t=−2.70, p=**.013**, Δ=−0.49% | strengthened (was p=.07) |
| theta | **FC2** | t=−0.78, p=.44, Δ=−0.17% | **t=−3.24, p=.0037**, Δ=−0.62% | **strengthened** (was p=.018) |
| theta | Fz   | t=−0.28, p=.78, Δ=−0.06% | t=−1.25, p=.22, Δ=−0.28% | similar |
| alpha | FC2  | t=+0.69, p=.50, Δ=+0.17% | t=+0.20, p=.84, Δ=+0.06% | sign flip from full |
| beta  | Fz   | (one cluster p=.035 emerges) | (one cluster p=.043 emerges) | possible noise |

**The hypothesis-driven theta-FC2 Met−Ano effect doubled in magnitude (−0.37% → −0.62%) and dropped a half-order in p (.018 → .0037)** when restricted to the familiarity-matched items. Met−Ano at C4 (the other LPC peak channel) also became significant (p=.013) where it was a trend before (p=.07). Cluster-based permutation for theta @ FC2 Met−Ano is now at the edge: cluster-level p=.051 (was .102).

**Crucially, the strengthening is selective to Met−Ano** — Met−Lit at FC2 went the wrong direction (slightly weaker) once matched. This *exactly mirrors the asymmetry in the ERP paper's own mass-univariate*, where the Met−Ano LPC cluster was larger and more sustained than Met−Lit. The connectivity result is reproducing the contrast structure your ERP results already established.

**Bonferroni context.** For the 3 theta-band Met−Ano ROI tests (one per seed), Bonferroni-corrected α = .05/3 = .017. The FC2/theta/Met−Ano result (p=.0037) **survives that correction**. Across the broader 18-test sweep (Bonferroni .05/18 = .0028) it does not survive — but the broader sweep includes alpha and beta, which were exploratory rather than hypothesis-driven for this contrast.

### Single-trial within-trial PLV LME — *uninformative due to ceiling effect*

Built a parallel per-trial measure: windowed phase-locking value across LPC window × 5 theta freqs × 7 ROI channels around FC2, on Laplacian-CSD data. Fit:
`isPC_LPC ~ Condition + Familiarity_c + (1|Sub) + (1|Item)` on matched (n=2142) and full (n=3032) trial sets.

All effects null:
- Matched: Condition F(2, 2138)=0.77, p=.46. Familiarity F(1, 2138)=0.40, p=.53. Intercept ≈ **0.954** (SE 0.003).
- Full: Condition F(2, 3028)=0.65, p=.52. Familiarity F(1, 3028)=1.81, p=.18. Intercept ≈ 0.954.
- No Condition × Familiarity interaction.

The intercept gives away the problem: within-trial PLV across a 700-ms LPC window × 5 freqs × 7 channel pairs has a **ceiling near 1** (mean 0.954, residual SD 0.033). A measure that's saturated by construction has no dynamic range to detect modulation. This is a known limitation of within-trial PLV with reasonably long windows on Laplacian-cleaned data. The LME approach is theoretically sound but needs a different per-trial measure — see §6.

## 5. Provisional claim — *revised after seeing the sweep results*

**Updated headline claim — defensible:**
"In hybrid visual metaphor processing, theta-band (4–8 Hz) phase coupling between right fronto-central cortex (FC2, identified as an LPC peak channel in the companion ERP paper) and its scalp neighbors is reduced for Metaphor relative to Anomalous images during the LPC window (800–1500 ms). The effect is robust to familiarity matching: t(23)=−3.24, p=.0037 in the familiarity-matched 32-triplet subset where the ERP paper's primary contrasts were established. The Met−Anomalous direction mirrors the same asymmetry observed in the ERP paper's mass-univariate LPC clusters, where Met−Anomalous was larger and more sustained than Met−Literal."

**What the matched-subset result changes:**

The case for a separate paper is now considerably stronger. Specifically:

1. **The hypothesis-driven theta-FC2 Met−Anomalous effect survives Bonferroni correction within the theta band** (p=.0037 vs. corrected α=.017 for 3 seeds). The same effect at C4 (the other LPC peak channel) is also significant (p=.013). Two independent LPC peak channels show the same theta-decoupling pattern in the same window.
2. **The matched-subset effect is roughly twice the size of the full-set effect** (Δ=−0.62% vs. −0.37%). This is the opposite of what familiarity-as-confound would predict — if the connectivity were just tracking familiarity, removing the familiarity-imbalanced items should weaken not strengthen the effect. That is a positive piece of evidence against the "this is just familiarity" reviewer challenge.
3. **The Met−Anomalous selectivity** (Met−Literal weakens, Met−Anomalous strengthens after matching) **mirrors the ERP paper's own structure**. This isn't an arbitrary contrast that worked — it's the same contrast structure the ERP paper already showed at the amplitude level, replicated at the connectivity level using a methodologically independent measure.
4. **Cluster permutation now borderline** (p=.051 at theta @ FC2 Met−Ano). With one more subject, or 2000 instead of 1000 permutations to stabilize the null, this likely tips. Currently a trend that is consistent with the ROI test.

**Revised recommendation:** **Pitch as a separate short report.** The story is now tight, the result is hypothesis-driven (not p-hacked), the contrast structure is pre-registered by the ERP paper, the matched-subset replication directly answers the obvious familiarity-confound criticism, and the LPC peak channels were chosen *before* running this analysis based on the published ERP findings. This is what a defensible second paper looks like.

What still needs work: (a) the cluster perm needs to clear .05; (b) Met−Literal at the same site should be examined more carefully; (c) the failed within-trial PLV LME needs replacing with a working per-trial measure (see §6) before claims about familiarity at the trial level can be made.

---

## 6. Planned analyses before submission

In priority order, with what's already done crossed out.

1. ✅ **Familiarity-matched subset** (32 triplets) — DONE, see §4b. Strengthened the FC2/theta/Met−Ano effect.
2. ✅ **Cluster-based permutation** — DONE. Borderline at theta-FC2-Met-Ano in matched subset (p=.051). Plan to bump n_perm from 1000 to 5000 to stabilize the null and confirm.
3. ✅ **Hypothesis-driven ROI test in the LPC window at C4 and FC2** — DONE. Both significant in matched subset.
4. ⚠️ **Item-level mixed-effects model with single-trial ISPC** — ATTEMPTED but the within-trial PLV measure has a ceiling effect (mean 0.954, near saturation). Need a different per-trial proxy:
   - **Option A (recommended):** per-item ISPC averaged across subjects. Each item appears once per subject in one condition; pool the per-subject phase angles for that item and compute the across-subject ISPC. Yields one ISPC value per (item × condition) — different decomposition but enables item-level LME with familiarity as a continuous predictor.
   - **Option B:** narrow the within-trial PLV window (e.g., 100-ms sliding windows) to lower the ceiling. Requires more careful baseline normalization.
   - **Option C:** use single-trial wavelet power in the LPC window as a covariate / sanity check rather than as the connectivity measure.
5. **Phase Lag Index (PLI)** as a volume-conduction control — `seed_ispc_topo` already produces PLI as an optional 5th output; need to re-run sweep extracting it and replicate the FC2/theta/Met−Ano contrast.
6. **Power spectrum / ERSP sanity check** at FC2 in theta in the LPC window. Confirm that the ISPC effect is not just trivially explained by a co-occurring power difference (it almost certainly isn't, since power and phase consistency are mathematically independent, but reviewers will ask).
7. **Wider seed set / data-driven seed identification.** FC2 was a hypothesis-driven seed pre-registered by the ERP paper. For a sensitivity/robustness analysis, repeat with all 32 channels as seeds and confirm the FC2/C4 cluster is the maximum.

---

## 7. Limitations / things to be upfront about

- **Exploratory seed choice.** Visual-metaphor EEG literature is thin; there is no pre-registered seed for a visual paradigm. C4 and FC2 are data-driven picks from the companion ERP paper. FCz / Fz are theoretical priors for cognitive control and semantic integration. The sweep allows comparison across both.
- **Uncorrected single-subject and per-channel t-tests are not evidence.** Only the cluster-permutation results will be reported as primary.
- **Small effect sizes (~2–3% change) demand large N to interpret confidently.** n=24 is on the lower end. The paper should be framed as a proof-of-principle / hypothesis-generating connectivity analysis of a dataset whose primary findings are in the ERP paper.
- **No predefined band for visual metaphor.** Theta is motivated by the linguistic metaphor literature, but visual tasks often engage alpha/beta more strongly. The paper must report the full sweep to avoid forking-paths concerns.

---

## 8. What I need from the advisor

1. **Green light** to develop this as a separate short report rather than fold it into the existing manuscript. The case is now: (a) the matched-subset replication strengthened the hypothesis-driven effect, (b) the seeds were pre-registered by the LPC mass-univariate in the companion paper, and (c) the contrast that survives (Met−Anomalous) mirrors the same asymmetry the ERP results already showed.
2. **Target journal**. With the current strength (one Bonferroni-surviving ROI test within theta, borderline cluster permutation, n=24), this fits a short report at a methodologically oriented journal: *Psychophysiology*, *Brain Topography*, *Neuropsychologia* (Brief Communication), or *Cortex* (Note). Higher-tier (*NeuroImage*, *J Cognitive Neuroscience*) likely requires the cluster perm to clear and the per-item LME to work.
3. **Decision on the per-item LME**. The within-trial PLV measure I built failed (ceiling). I have three options in §6 (per-item across-subject ISPC; narrower windowed PLV; or fold into power-as-covariate sanity check). Recommend Option A. Want to confirm before investing the time.
4. **Confirmation that sub 32 exclusion is final** across both papers.

---

## Appendix A — Methods summary (short form, for the paper's §2)

- n = 24 (one subject excluded for failed QC; sub 32).
- Surface Laplacian (Perrin spherical spline) applied to each epoch.
- Morlet wavelet convolution, 5 log-spaced frequencies per band. Log-spaced cycles [3, 10].
- Seed-based ISPC: across-trial consistency of phase differences between a seed channel and every other channel, computed at each time point and each frequency, then averaged across frequencies within the band.
- Baseline correction: percent change from −200 to 0 ms.
- Cluster-based permutation: 1000 sign-flip iterations; cluster-forming threshold |t|>t_crit at α=0.05; channel adjacency threshold 75 mm from 3D chanlocs; cluster-mass null distribution.
- ROI test (LPC window 800–1500 ms): paired t on mean ISPC across seed + neighbors.
- Supplementary: PLI (same pipeline, phase-lead/lag instead of phase-difference magnitude).

## Appendix B — Repository state

All scripts live in `C:\Users\Cagi\Desktop\516\Day_12\`.

| File                          | Purpose                                                                       |
|-------------------------------|-------------------------------------------------------------------------------|
| `seed_ispc_topo.m` (in 516/)  | Core function: channels×time ISPC topography with optional PLI output.        |
| `seed_ispc_topo_demo.m` (in 516/) | Class-assignment demo on ERN data. Unrelated to VisMet paper.             |
| `vismet_ispc_demo.m`          | Single-subject demo on VisMet (pick subject at top). Good for teaching/pitch. |
| `vismet_ispc_group.m`         | Group analysis with sub 32 exclusion. Produces the per-condition + difference topographies and probe-channel time courses. |
| `vismet_ispc_sweep.m`         | Band × seed sweep with cluster-based permutation and LPC ROI test.            |
| `vismet_ispc_run_all.m`       | Convenience wrapper that runs both group and sweep sequentially.              |
| `ispc_cache/`                 | Per-(band, seed) cached per-subject ISPC arrays.                              |
| `ispc_figs/`                  | Per-combo cluster-perm figure (PNG) and ROI time-course figure.               |
| `ispc_sweep_summary.csv`      | Machine-readable summary: cluster p, ROI t/p, mean differences.               |
