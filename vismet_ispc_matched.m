%% ========================================================================
%  VISMET ISPC — MATCHED-SUBSET sweep
%  ========================================================================
%  Same band x seed sweep as vismet_ispc_sweep, but restricted to the
%  familiarity-matched 32-triplet subset (per Cagatay's preprocessing:
%  matched trials = ODD bin numbers, preceded by trigger 85; unmatched =
%  EVEN bin numbers, preceded by trigger 88).
%
%  Why: the manuscript's headline ERP finding (LPC Met>Lit, Met>Ano)
%  survives in the matched subset where familiarity is equated by design.
%  Replicating that dissociation in theta-band ISPC is the cleanest
%  follow-up.
%
%  Caches to ispc_cache_matched/ . Same exclude list (sub 32).
% =========================================================================

clear; close all; clc;

%% Paths
study_dir  = 'C:\Users\Cagi\Desktop\VisualMet';
set_dir    = fullfile(study_dir, 'EEGSets');
func_dir   = 'C:\Users\Cagi\Desktop\516';
out_dir    = 'C:\Users\Cagi\Desktop\516\Day_12';
cache_dir  = fullfile(out_dir, 'ispc_cache_matched');
fig_dir    = fullfile(out_dir, 'ispc_figs_matched');
if ~exist(cache_dir,'dir'), mkdir(cache_dir); end
if ~exist(fig_dir,'dir'),   mkdir(fig_dir);   end
addpath(func_dir);
if exist('eeglab','file') ~= 2
    error('EEGLAB not on path.');
end

%% Sweep grid (same as the all-trials sweep)
BANDS = { [4 8], [8 12], [13 25] };
BAND_NAMES = {'theta','alpha','beta'};
SEEDS = {'C4','FC2','FCz'};

num_frex        = 5;
time_range      = [-200 1500];
baseline_window = [-200 0];
baseline_flag   = 1;
use_laplacian   = true;

n_perm        = 1000;
cluster_alpha = 0.05;
neighbor_dist = 75;
lpc_window    = [800 1500];

EXCLUDE = {'32'};

% Subjects
files = dir(fullfile(set_dir,'MetVis_*_postart.set'));
sub_ids = cell(1, numel(files));
for k = 1:numel(files)
    tok = regexp(files(k).name,'MetVis_(\d+)_postart\.set','tokens','once');
    if ~isempty(tok), sub_ids{k} = tok{1}; end
end
sub_ids = sub_ids(~cellfun(@isempty, sub_ids));
sub_ids = sub_ids(~ismember(sub_ids, EXCLUDE));
fprintf('Found %d subjects (after EXCLUDE = %s).\n', numel(sub_ids), strjoin(EXCLUDE,','));

results = struct([]);

for bi = 1:numel(BANDS)
    for si = 1:numel(SEEDS)
        band_name = BAND_NAMES{bi};
        seed_req  = SEEDS{si};
        tag = sprintf('%s_%s_matched', band_name, seed_req);
        cache_file = fullfile(cache_dir, sprintf('ispc_%s.mat', tag));

        fprintf('\n========== %s (%g-%g Hz) at %s — MATCHED ==========\n', ...
            upper(band_name), BANDS{bi}(1), BANDS{bi}(2), seed_req);

        if exist(cache_file,'file')
            fprintf('Loading cache: %s\n', cache_file);
            S = load(cache_file);
            ispc_group = S.ispc_group; chanlocs = S.chanlocs;
            times_out  = S.times_out;  sub_ids_used = S.sub_ids_used;
            seed_used  = S.seed_used;  cond_names = S.cond_names;
            n_trials_mat = S.n_trials_mat;
            drop = ismember(sub_ids_used, EXCLUDE);
            if any(drop)
                ispc_group = ispc_group(~drop,:,:,:);
                sub_ids_used = sub_ids_used(~drop);
                n_trials_mat = n_trials_mat(~drop,:);
            end
        else
            [ispc_group, chanlocs, times_out, sub_ids_used, seed_used, ...
             cond_names, n_trials_mat] = ...
                local_compute_matched_ispc(sub_ids, set_dir, seed_req, ...
                BANDS{bi}, num_frex, time_range, baseline_window, ...
                baseline_flag, use_laplacian);
            save(cache_file, 'ispc_group','chanlocs','times_out', ...
                'sub_ids_used','seed_used','cond_names','n_trials_mat', '-v7.3');
            fprintf('Cached %s (%d subs, mean trials L/M/A = %.1f/%.1f/%.1f).\n', ...
                tag, size(ispc_group,1), mean(n_trials_mat,1));
        end

        adj = local_chan_adjacency(chanlocs, neighbor_dist);

        fprintf('Cluster perm Met-Lit (n=%d)...\n', n_perm);
        [t_ML, clus_ML] = local_cluster_perm( ...
            squeeze(ispc_group(:,:,:,2)), squeeze(ispc_group(:,:,:,1)), ...
            adj, cluster_alpha, n_perm);
        fprintf('Cluster perm Met-Ano (n=%d)...\n', n_perm);
        [t_MA, clus_MA] = local_cluster_perm( ...
            squeeze(ispc_group(:,:,:,2)), squeeze(ispc_group(:,:,:,3)), ...
            adj, cluster_alpha, n_perm);

        seed_idx = find(strcmpi({chanlocs.labels}, seed_used), 1);
        roi_mask = false(numel(chanlocs),1); roi_mask(seed_idx) = true;
        roi_mask(adj(seed_idx,:)) = true;
        roi_chans = {chanlocs(roi_mask).labels};
        t_idx_lpc = times_out >= lpc_window(1) & times_out <= lpc_window(2);

        roi_lit = local_nanmean(local_nanmean(ispc_group(:,roi_mask,t_idx_lpc,1),2),3);
        roi_met = local_nanmean(local_nanmean(ispc_group(:,roi_mask,t_idx_lpc,2),2),3);
        roi_ano = local_nanmean(local_nanmean(ispc_group(:,roi_mask,t_idx_lpc,3),2),3);
        roi_lit = roi_lit(:); roi_met = roi_met(:); roi_ano = roi_ano(:);
        [tML, pML] = local_paired_t_1d(roi_met, roi_lit);
        [tMA, pMA] = local_paired_t_1d(roi_met, roi_ano);

        local_plot_combo(ispc_group, chanlocs, times_out, seed_used, ...
            band_name, BANDS{bi}, t_ML, clus_ML, t_MA, clus_MA, ...
            roi_chans, lpc_window, fig_dir, tag);

        row.band = band_name; row.band_hz = BANDS{bi}; row.seed = seed_used;
        row.n_subs = size(ispc_group,1);
        row.mean_trials_LMA = mean(n_trials_mat,1);
        row.ML_min_p = min([clus_ML.p_cluster, NaN]);
        row.ML_n_sig = sum([clus_ML.p_cluster] < cluster_alpha);
        row.MA_min_p = min([clus_MA.p_cluster, NaN]);
        row.MA_n_sig = sum([clus_MA.p_cluster] < cluster_alpha);
        row.roi_chans = strjoin(roi_chans,',');
        row.ML_roi_t = tML; row.ML_roi_p = pML; row.ML_roi_meanDiff = mean(roi_met-roi_lit);
        row.MA_roi_t = tMA; row.MA_roi_p = pMA; row.MA_roi_meanDiff = mean(roi_met-roi_ano);
        results = [results, row]; %#ok<AGROW>
    end
end

fprintf('\n\n================ MATCHED-SUBSET SUMMARY ================\n');
T = struct2table(results);
disp(T(:, {'band','seed','n_subs','ML_min_p','ML_n_sig','ML_roi_t','ML_roi_p','ML_roi_meanDiff', ...
          'MA_min_p','MA_n_sig','MA_roi_t','MA_roi_p','MA_roi_meanDiff'}));
writetable(T, fullfile(out_dir, 'ispc_sweep_summary_matched.csv'));


% =========================================================================
%  Compute per-subject ISPC restricted to MATCHED trials (odd bin numbers)
% =========================================================================
function [ispc_group, chanlocs_out, times_out, sub_ids_used, seed_used, cond_names, n_trials_mat] = ...
    local_compute_matched_ispc(sub_ids, set_dir, seed_req, band, num_frex, ...
    time_range, baseline_window, baseline_flag, use_laplacian)

cond_names  = {'Literal','Metaphor','Anomalous'};
valid_codes = [111 121 131 112 122 132 113 123 133];
n_conds     = 3;

ispc_group = []; chanlocs_out = []; times_out = []; sub_ids_used = {}; seed_used = '';
n_trials_mat = nan(numel(sub_ids), n_conds);

for s = 1:numel(sub_ids)
    sid = sub_ids{s};
    fname = sprintf('MetVis_%s_postart.set', sid);
    fpath = fullfile(set_dir, fname);
    if ~exist(fpath,'file'), continue; end
    fprintf('  [%2d/%d] sub %s ', s, numel(sub_ids), sid);
    EEG = pop_loadset('filename', fname, 'filepath', set_dir);

    sc = seed_req;
    if ~any(strcmpi({EEG.chanlocs.labels}, sc))
        if any(strcmpi({EEG.chanlocs.labels},'Fz')), sc = 'Fz';
        else, fprintf('-- no seed/Fz, skip\n'); continue; end
    end

    % Classify trials AND determine matched/unmatched from bin parity
    cond_of_epoch  = nan(1, EEG.trials);
    is_matched     = false(1, EEG.trials);
    for ep = 1:EEG.trials
        % Get stim-code condition from time-lock event
        types = EEG.epoch(ep).eventtype;
        lats  = EEG.epoch(ep).eventlatency;
        if ~iscell(types), types = {types}; end
        if ~iscell(lats),  lats  = {lats};  end
        numlats = cell2mat(lats);
        numtypes = nan(1, numel(types));
        for k = 1:numel(types)
            numtypes(k) = local_extract_code(types{k}, valid_codes);
        end
        keep = ~isnan(numtypes);
        if ~any(keep), continue; end
        numtypes_kept = numtypes(keep); numlats_kept = numlats(keep);
        [~, mi] = min(abs(numlats_kept));
        switch mod(numtypes_kept(mi),10)
            case 1, cond_of_epoch(ep) = 1;
            case 2, cond_of_epoch(ep) = 2;
            case 3, cond_of_epoch(ep) = 3;
        end

        % Get bin number for the time-lock event from eventbini.
        % eventbini may be a cell array (one per event in epoch).
        bin_num = NaN;
        if isfield(EEG.epoch(ep), 'eventbini')
            bn = EEG.epoch(ep).eventbini;
            if ~iscell(bn), bn = {bn}; end
            % find event index aligned to time-lock (the same event we just used)
            % map: event in 'kept' positions is at original idx find(keep)
            kept_pos = find(keep);
            tl_event_in_epoch = kept_pos(mi);
            if numel(bn) >= tl_event_in_epoch
                v = bn{tl_event_in_epoch};
                if iscell(v), v = v{1}; end
                if isnumeric(v) && ~isempty(v)
                    bin_num = v(1);
                end
            end
        end
        % matched = odd bin number per the BDF convention
        if ~isnan(bin_num) && bin_num > 0
            is_matched(ep) = mod(bin_num,2) == 1;
        end
    end

    % Restrict to MATCHED trials only
    keep_mask = is_matched;
    cond_of_epoch_m = cond_of_epoch;
    cond_of_epoch_m(~keep_mask) = NaN;
    n_per = [sum(cond_of_epoch_m==1), sum(cond_of_epoch_m==2), sum(cond_of_epoch_m==3)];
    n_trials_mat(s,:) = n_per;
    if min(n_per) < 8, fprintf('-- matched trials L/M/A=%d/%d/%d, skip\n', n_per); continue; end

    if use_laplacian
        X = [EEG.chanlocs.X]; Y = [EEG.chanlocs.Y]; Z = [EEG.chanlocs.Z];
        data_for_ispc = zeros(size(EEG.data));
        for tri = 1:EEG.trials
            [sl, ~, ~] = laplacian_perrinX(EEG.data(:,:,tri), X, Y, Z);
            data_for_ispc(:,:,tri) = sl;
        end
    else
        data_for_ispc = EEG.data;
    end

    ispc_subj = [];
    for c = 1:n_conds
        idx = (cond_of_epoch_m == c);
        data_c = data_for_ispc(:,:,idx);
        [ispc_c, ~, t_out] = seed_ispc_topo( ...
            data_c, EEG.srate, EEG.times, EEG.chanlocs, sc, ...
            band, num_frex, time_range, baseline_window, ...
            baseline_flag, 0, 0);
        if isempty(ispc_subj)
            ispc_subj = nan(size(ispc_c,1), size(ispc_c,2), n_conds);
        end
        ispc_subj(:,:,c) = ispc_c;
    end

    if isempty(ispc_group)
        ispc_group   = nan(numel(sub_ids), size(ispc_subj,1), size(ispc_subj,2), n_conds);
        chanlocs_out = EEG.chanlocs;
        times_out    = t_out;
        seed_used    = sc;
    end
    if size(ispc_subj,1) ~= size(ispc_group,2) || size(ispc_subj,2) ~= size(ispc_group,3)
        fprintf('-- size mismatch, skip\n'); continue;
    end
    ispc_group(s,:,:,:) = ispc_subj;
    sub_ids_used{end+1} = sid; %#ok<AGROW>
    fprintf('matched L/M/A=%d/%d/%d done.\n', n_per);
end

keep = squeeze(any(any(any(~isnan(ispc_group),2),3),4));
ispc_group = ispc_group(keep,:,:,:);
n_trials_mat = n_trials_mat(keep,:);
end


% =========================================================================
%  Helper functions reused from vismet_ispc_sweep.m
% =========================================================================
function [t_obs, clus] = local_cluster_perm(A, B, adj, alpha, n_perm)
D    = A - B;
n    = sum(~isnan(D),1);
m    = local_nanmean(D, 1);
sd   = local_nanstd(D, 1);
t_obs = squeeze(m ./ (sd ./ sqrt(n)));
df    = squeeze(n) - 1;
df0   = max(df(:));
t_thresh = abs(tinv(alpha/2, df0));
[obs_pos, obs_neg] = local_find_clusters(t_obs, t_thresh, adj);
n_subs = size(D,1);
max_pos = nan(n_perm,1); max_neg = nan(n_perm,1);
for p = 1:n_perm
    flip = (rand(n_subs,1) < 0.5)*2 - 1;
    Dp = D .* flip;
    mp = local_nanmean(Dp,1); sdp = local_nanstd(Dp,1);
    tp = squeeze(mp ./ (sdp ./ sqrt(n)));
    [pp, np] = local_find_clusters(tp, t_thresh, adj);
    if isempty(pp), max_pos(p) = 0; else, max_pos(p) = max(arrayfun(@(c)c.mass,pp)); end
    if isempty(np), max_neg(p) = 0; else, max_neg(p) = max(arrayfun(@(c)c.mass,np)); end
end
clus = struct([]);
for k = 1:numel(obs_pos)
    obs_pos(k).sign = +1;
    obs_pos(k).p_cluster = (1+sum(max_pos>=obs_pos(k).mass))/(1+n_perm);
    clus = [clus, obs_pos(k)]; %#ok<AGROW>
end
for k = 1:numel(obs_neg)
    obs_neg(k).sign = -1;
    obs_neg(k).p_cluster = (1+sum(max_neg>=obs_neg(k).mass))/(1+n_perm);
    clus = [clus, obs_neg(k)]; %#ok<AGROW>
end
end

function [pos_clus, neg_clus] = local_find_clusters(t_map, thresh, adj)
n_chans = size(t_map,1); n_times = size(t_map,2);
mask_pos = t_map >  thresh;
mask_neg = t_map < -thresh;
pos_clus = local_extract(mask_pos, t_map, adj, n_chans, n_times);
neg_clus = local_extract(mask_neg, t_map, adj, n_chans, n_times);
end

function clusters = local_extract(mask, tmap, adj, n_chans, n_times)
clusters = struct('channels',{},'times',{},'mass',{});
unvisited = mask;
while any(unvisited(:))
    [ci, ti] = find(unvisited, 1, 'first');
    queue = [ci, ti];
    chans_in = []; times_in = []; mass = 0;
    while ~isempty(queue)
        c = queue(1,1); t = queue(1,2); queue(1,:) = [];
        if c<1||c>n_chans||t<1||t>n_times, continue; end
        if ~unvisited(c,t), continue; end
        unvisited(c,t) = false;
        chans_in(end+1) = c; %#ok<AGROW>
        times_in(end+1) = t; %#ok<AGROW>
        mass = mass + abs(tmap(c,t));
        queue(end+1,:) = [c, t-1]; %#ok<AGROW>
        queue(end+1,:) = [c, t+1]; %#ok<AGROW>
        nb = find(adj(c,:));
        for k = 1:numel(nb)
            queue(end+1,:) = [nb(k), t]; %#ok<AGROW>
        end
    end
    clusters(end+1).channels = chans_in; %#ok<AGROW>
    clusters(end).times = times_in;
    clusters(end).mass  = mass;
end
end

function adj = local_chan_adjacency(chanlocs, max_dist)
n = numel(chanlocs); xyz = nan(n,3);
for k = 1:n
    if ~isempty(chanlocs(k).X) && ~isempty(chanlocs(k).Y) && ~isempty(chanlocs(k).Z)
        xyz(k,:) = [chanlocs(k).X, chanlocs(k).Y, chanlocs(k).Z];
    end
end
norms = sqrt(sum(xyz.^2,2));
mn = median(norms(~isnan(norms)));
if mn < 5, xyz = xyz * (90/mn); end
D = nan(n);
for i = 1:n
    for j = 1:n
        D(i,j) = norm(xyz(i,:)-xyz(j,:));
    end
end
adj = (D > 0) & (D <= max_dist);
end

function local_plot_combo(ispc_group, chanlocs, times_out, seed_used, ...
    band_name, band_hz, t_ML, clus_ML, t_MA, clus_MA, ...
    roi_chans, lpc_window, fig_dir, tag)
n_subs = size(ispc_group,1);
n_plots = 8;
t_idx_p = round(linspace(1,length(times_out),n_plots));
seed_idx = find(strcmpi({chanlocs.labels}, seed_used),1);
cmap_div = local_diverging_cmap();
sig_ML = false(numel(chanlocs), length(times_out));
for k = 1:numel(clus_ML)
    if clus_ML(k).p_cluster < 0.05
        for j = 1:numel(clus_ML(k).channels)
            sig_ML(clus_ML(k).channels(j), clus_ML(k).times(j)) = true;
        end
    end
end
sig_MA = false(numel(chanlocs), length(times_out));
for k = 1:numel(clus_MA)
    if clus_MA(k).p_cluster < 0.05
        for j = 1:numel(clus_MA(k).channels)
            sig_MA(clus_MA(k).channels(j), clus_MA(k).times(j)) = true;
        end
    end
end
clim_t = max([abs(t_ML(:)); abs(t_MA(:))]) * 0.9;
if ~isfinite(clim_t)||clim_t==0, clim_t=3; end

fig = figure('Color','w','Position',[20 20 1500 700], ...
    'Name', sprintf('MATCHED — %s @ %s (n=%d)', band_name, seed_used, n_subs));
sgtitle(sprintf(['MATCHED  n=%d  |  seed=%s  |  %s (%g-%g Hz)  |  ', ...
    'cluster t-maps; black dots = p_{cluster}<0.05'], ...
    n_subs, seed_used, band_name, band_hz(1), band_hz(2)), ...
    'FontSize', 12, 'FontWeight','bold');

row_titles = {'Met − Lit  (t)', 'Met − Ano  (t)'};
row_t  = {t_ML, t_MA};
row_sig= {sig_ML, sig_MA};
for r = 1:2
    for p = 1:n_plots
        ax = subplot(2, n_plots, (r-1)*n_plots + p);
        ti = t_idx_p(p);
        emk = {seed_idx,'o','r',6,1.5};
        sig_here = find(row_sig{r}(:,ti));
        if ~isempty(sig_here), emk = {sig_here,'.','k',12,1}; end
        topoplot(row_t{r}(:,ti), chanlocs, ...
            'maplimits',[-clim_t clim_t],'electrodes','on', ...
            'style','both','numcontour',6,'emarker2',emk);
        colormap(ax, cmap_div);
        if r == 1, title(sprintf('%d ms', round(times_out(ti))),'FontSize',9); end
        if p == 1
            pos = get(ax,'Position');
            annotation('textbox', ...
                [0.001, pos(2)+pos(4)/2 - 0.025, pos(1)-0.005, 0.05], ...
                'String', row_titles{r}, 'FontSize', 10, 'FontWeight','bold', ...
                'EdgeColor','none','HorizontalAlignment','right','VerticalAlignment','middle');
        end
        if p == n_plots
            cb = colorbar(ax,'eastoutside'); cb.Label.String='t'; cb.FontSize=8;
        end
    end
end
exportgraphics(fig, fullfile(fig_dir, sprintf('cluster_%s.png', tag)), 'Resolution', 150);

fig2 = figure('Color','w','Position',[40 40 900 400]);
ch_idx = ismember({chanlocs.labels}, roi_chans);
m_lit = squeeze(local_nanmean(local_nanmean(ispc_group(:,ch_idx,:,1),1),2));
m_met = squeeze(local_nanmean(local_nanmean(ispc_group(:,ch_idx,:,2),1),2));
m_ano = squeeze(local_nanmean(local_nanmean(ispc_group(:,ch_idx,:,3),1),2));
plot(times_out, m_lit,'b-', times_out, m_met,'r-', times_out, m_ano,'y-', 'LineWidth', 1.6); hold on;
xline(0,'k--'); yline(0,'k:');
patch([lpc_window(1) lpc_window(2) lpc_window(2) lpc_window(1)], ...
      [min(ylim) min(ylim) max(ylim) max(ylim)], [0.85 0.85 0.85], ...
      'FaceAlpha', 0.25, 'EdgeColor','none');
xlabel('Time (ms)'); ylabel('ISPC (% change)');
title(sprintf('MATCHED — ROI mean (%s) — %s @ %s', strjoin(roi_chans,','), band_name, seed_used));
legend({'Lit','Met','Ano','stim','zero','LPC'},'Location','best');
exportgraphics(fig2, fullfile(fig_dir, sprintf('roi_%s.png',tag)), 'Resolution',150);
end

function code = local_extract_code(x, valid_codes)
    code = NaN;
    if isnumeric(x) && ~isempty(x)
        if ismember(x, valid_codes), code = x; end; return;
    end
    if ~ischar(x) || isempty(x), return; end
    s = x; if s(1)=='S', s = s(2:end); end
    if ~isempty(s) && all(isstrprop(s,'digit'))
        v = str2double(s);
        if ismember(v, valid_codes), code = v; end; return;
    end
    toks = regexp(x,'\d+','match');
    for k = 1:numel(toks)
        v = str2double(toks{k});
        if ismember(v, valid_codes), code = v; return; end
    end
end

function m = local_nanmean(X, dim)
    valid = ~isnan(X); Xz = X; Xz(~valid) = 0;
    n = sum(valid, dim); s = sum(Xz, dim);
    m = s ./ n; m(n==0) = NaN;
end

function sd = local_nanstd(X, dim)
    valid = ~isnan(X); n = sum(valid, dim);
    mu = local_nanmean(X, dim);
    Xc = X - mu; Xc(~valid) = 0;
    ss = sum(Xc.^2, dim);
    sd = sqrt(ss ./ max(n-1,1));
    sd(n < 2) = NaN;
end

function [t, p] = local_paired_t_1d(a, b)
    d = a - b; d = d(~isnan(d));
    n = numel(d); m = mean(d); sd = std(d);
    t = m / (sd / sqrt(n));
    p = 2 * (1 - tcdf(abs(t), n-1));
end

function cmap = local_diverging_cmap()
n = 128;
low  = linspace(0,1,n+1); high = linspace(1,0,n+1);
r = [low(1:end-1), ones(1,n)];
g = [low(1:end-1), high(2:end)];
b = [ones(1,n),    high(2:end)];
cmap = [r(:), g(:), b(:)];
end
