%% ========================================================================
%  RUN EVERYTHING: group ISPC analysis + (band x seed) sweep with
%  cluster-based permutation tests. Excludes subject 32.
%  ========================================================================
%  Prereqs:
%     - EEGLAB already started in this MATLAB session (for pop_loadset,
%       topoplot, laplacian_perrinX).
%     - C:\Users\Cagi\Desktop\516 is on path or visible (holds
%       seed_ispc_topo.m and laplacian_perrinX from class code).
%
%  What this does (in order):
%     1) vismet_ispc_group    — theta-band FCz (fallback Fz) group analysis.
%                                Produces per-condition topo grid, Met-Lit /
%                                Met-Ano difference maps, and probe-channel
%                                time courses across 9 electrodes.
%     2) vismet_ispc_sweep    — 9 combos of (band x seed) with cluster-based
%                                permutation on Met-Lit and Met-Ano plus an
%                                ROI-restricted LPC-window test.
%                                Writes ispc_sweep_summary.csv +
%                                per-combo PNG figures in ispc_figs/.
%
%  To regenerate everything from scratch (e.g., if you change the excluded
%  subject list), DELETE the caches first:
%     - Day_12\vismet_ispc_group_cache.mat
%     - Day_12\ispc_cache\*.mat
% =========================================================================

fprintf('\n===== STAGE 1: Group analysis (theta, FCz) =====\n');
run('vismet_ispc_group.m');

fprintf('\n===== STAGE 2: Band x Seed sweep with cluster permutation =====\n');
run('vismet_ispc_sweep.m');

fprintf('\n===== ALL DONE =====\n');
fprintf('Summary CSV:   C:\\Users\\Cagi\\Desktop\\516\\Day_12\\ispc_sweep_summary.csv\n');
fprintf('Per-combo figs: C:\\Users\\Cagi\\Desktop\\516\\Day_12\\ispc_figs\\\n');
fprintf('Group figs are still in open figure windows.\n');
