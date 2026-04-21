# ISPC-VisMet

Seed-based inter-site phase clustering (ISPC) analysis of the VisMet hybrid visual-metaphor EEG dataset (n=24, subject 32 excluded). Companion connectivity analysis to an ERP paper currently in revision — tests whether Metaphor vs. Literal/Anomalous conditions differ in phase-based connectivity in the LPC window (800–1500 ms), independent of evoked amplitude.

See [VisMet_ISPC_PaperPitch.md](VisMet_ISPC_PaperPitch.md) for full scientific rationale and results.

## Headline result

Theta @ FC2, Metaphor − Anomalous, ROI-restricted test in the 800–1500 ms window: **t(23) = −2.56, p = .018, Δ = −0.37%**. Cluster-based permutation across the full topography did not survive correction (smallest cluster p = .076).

## Pipeline

All scripts are MATLAB. Requires EEGLAB on the path (for `pop_loadset`, `topoplot`, `laplacian_perrinX`) and the preprocessed `_postart.set` files from the parent study.

Entry point: [vismet_ispc_run_all.m](vismet_ispc_run_all.m) runs both stages.

| Stage | Script | Output |
|---|---|---|
| 1 | [vismet_ispc_group.m](vismet_ispc_group.m) | Group-level theta/FCz ISPC — per-condition topo grid, Met−Lit and Met−Ano difference maps, probe-channel time courses |
| 2 | [vismet_ispc_sweep.m](vismet_ispc_sweep.m) | 3 bands × 3 seeds (theta/alpha/beta × C4/FC2/Fz) with cluster permutation + ROI test. Writes [ispc_sweep_summary.csv](ispc_sweep_summary.csv) and figures to `ispc_figs/` |
| 3 | [vismet_ispc_matched.m](vismet_ispc_matched.m) | Trial-count-matched replication (see `ispc_figs_matched/`, [ispc_sweep_summary_matched.csv](ispc_sweep_summary_matched.csv)) |
| 4 | [vismet_ispc_lme.m](vismet_ispc_lme.m) | Single-trial LME follow-up on the theta/FC2 survivor ([SingleTrial_ISPC_LPC.csv](SingleTrial_ISPC_LPC.csv), [lme_results.txt](lme_results.txt)) |

A single-subject demo is in [vismet_ispc_demo.m](vismet_ispc_demo.m).

## Method summary

- Surface Laplacian (Perrin spherical spline) applied before ISPC to suppress volume conduction.
- Morlet wavelet convolution, 5 log-spaced frequencies per band, cycles log-spaced 3→10.
- Baseline correction: percent change from −200 to 0 ms.
- Conditions from ERPLAB bin codes: 111/121/131 Literal, 112/122/132 Metaphor, 113/123/133 Anomalous.
- Cluster-based permutation: 1000 sign-flip iterations, 75 mm channel adjacency, cluster-forming α = 0.05.
- ROI test: seed + neighbors, 800–1500 ms, paired t across subjects.

## Figures

- `ispc_figs/` — main sweep (9 band×seed combos, cluster + ROI panels).
- `ispc_figs_matched/` — trial-matched replication.
- [p300_topo_perm.png](p300_topo_perm.png) — P300 topo permutation (separate class exercise).

## Data (not in repo)

`.mat` files (raw-ish inputs and caches) and unpacked manuscript sources are gitignored. Regenerate by deleting `vismet_ispc_group_cache.mat` and `ispc_cache/` then rerunning Stage 1+2.

## Other files

- [Code_Day12.m](Code_Day12.m), [NO_PEEKING_Answers_Day12.m](NO_PEEKING_Answers_Day12.m), [tf_decomp.m](tf_decomp.m), [groupFX.m](groupFX.m), [permutationTesting.m](permutationTesting.m), [p300_topo_perm.m](p300_topo_perm.m) — class exercises / utilities used by the pipeline.
- [_manuscript_text.txt](_manuscript_text.txt) — extracted text of the companion ERP manuscript.
