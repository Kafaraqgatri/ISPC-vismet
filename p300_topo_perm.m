%% P300 topographic permutation test -- all channels
% Rare vs Frequent P3 amplitude difference, with MC-corrected thresholds
% derived from max/min across channels per iteration.

clear all; close all; clc

% Ensure EEGLab is on the path (try common locations if not already)
if ~exist('topoplot','file')
    try; eeglab nogui; catch; end
end

load('P300.mat');

P3_timewin = [300 500];
P3timeidx  = dsearchn(EEG.times', P3_timewin');

freq_epochs = [EEG.event.type] == 100;
rare_epochs = [EEG.event.type] == 101;

n_chans = EEG.nbchan;
n_freq  = sum(freq_epochs);
n_rare  = sum(rare_epochs);
n_total = n_freq + n_rare;

% observed: P3 amplitude per channel (max in window) for each condition
freq_erp_all = mean(EEG.data(:,:,freq_epochs), 3);   % chan x time
rare_erp_all = mean(EEG.data(:,:,rare_epochs), 3);

freq_P3 = mean(freq_erp_all(:, P3timeidx(1):P3timeidx(2)), 2);  % chan x 1
rare_P3 = mean(rare_erp_all(:, P3timeidx(1):P3timeidx(2)), 2);
obs_diff = rare_P3 - freq_P3;   % observed topo of rare-freq P3 diff

% permutation test
n_iter = 1000;
rng(42);

perm_diff = zeros(n_chans, n_iter);
perm_max  = zeros(n_iter,1);
perm_min  = zeros(n_iter,1);

all_idx = [find(freq_epochs) find(rare_epochs)];  % pool of trial indices

fprintf('Running %d permutations across %d channels...\n', n_iter, n_chans);
for n = 1:n_iter
    ord = all_idx(randperm(n_total));
    p_freq_idx = ord(1:n_freq);
    p_rare_idx = ord(n_freq+1:end);

    p_freq_erp = mean(EEG.data(:,:,p_freq_idx), 3);
    p_rare_erp = mean(EEG.data(:,:,p_rare_idx), 3);

    p_freq_P3 = mean(p_freq_erp(:, P3timeidx(1):P3timeidx(2)), 2);
    p_rare_P3 = mean(p_rare_erp(:, P3timeidx(1):P3timeidx(2)), 2);

    d = p_rare_P3 - p_freq_P3;
    perm_diff(:,n) = d;
    perm_max(n)    = max(d);
    perm_min(n)    = min(d);

    if mod(n,100)==0, fprintf('  iter %d/%d\n', n, n_iter); end
end

% thresholds
z_crit      = abs(norminv(0.025));
perm_mean_c = mean(perm_diff, 2);
perm_std_c  = std(perm_diff, 0, 2);
z_obs       = (obs_diff - perm_mean_c) ./ perm_std_c;

% uncorrected per-channel percentiles
upper_uncorr = prctile(perm_diff, 97.5, 2);
lower_uncorr = prctile(perm_diff,  2.5, 2);
sig_uncorr   = obs_diff > upper_uncorr | obs_diff < lower_uncorr;

% MC-corrected via max/min across channels
upper_corr = prctile(perm_max, 97.5);
lower_corr = prctile(perm_min,  2.5);
sig_corr   = obs_diff > upper_corr | obs_diff < lower_corr;

fprintf('\nMC-corrected thresholds: [%.3f, %.3f] uV\n', lower_corr, upper_corr);
fprintf('Channels significant (uncorrected): %d / %d\n', sum(sig_uncorr), n_chans);
fprintf('Channels significant (MC-corrected): %d / %d\n', sum(sig_corr), n_chans);

sig_labels = {EEG.chanlocs(sig_corr).labels};
fprintf('MC-corrected significant channels: %s\n', strjoin(sig_labels, ', '));

% -------- Plotting --------
fig = figure('Position',[100 100 1300 1000],'Color','w');
t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
title(t, 'P3 mean amplitude (300-500 ms) Rare - Freq — permutation test (1000 iter)');

clim_abs = max(abs(obs_diff));

nexttile;
topoplot(obs_diff, EEG.chanlocs, 'electrodes','labels','maplimits',[-clim_abs clim_abs]);
title('Observed Rare - Freq (\muV)');
colorbar;

nexttile;
zlim = max(4, max(abs(z_obs)));
topoplot(z_obs, EEG.chanlocs, 'electrodes','labels','maplimits',[-zlim zlim]);
title('Z-score vs permutation null');
colorbar;

nexttile;
mask_u = obs_diff; mask_u(~sig_uncorr) = 0;
topoplot(mask_u, EEG.chanlocs, 'electrodes','labels','maplimits',[-clim_abs clim_abs], ...
         'emarker2',{find(sig_uncorr),'o','w',8,2});
title(sprintf('Uncorrected p<.05 (%d chans)', sum(sig_uncorr)));
colorbar;

nexttile;
mask_c = obs_diff; mask_c(~sig_corr) = 0;
topoplot(mask_c, EEG.chanlocs, 'electrodes','labels','maplimits',[-clim_abs clim_abs], ...
         'emarker2',{find(sig_corr),'o','w',8,2});
title(sprintf('MC-corrected p<.05 (%d chans)', sum(sig_corr)));
colorbar;

saveas(fig, 'p300_topo_perm.png');
fprintf('Saved figure to p300_topo_perm.png\n');

save('p300_topo_perm_results.mat', 'obs_diff','z_obs','sig_uncorr','sig_corr', ...
     'upper_corr','lower_corr','perm_diff','-v7.3');
