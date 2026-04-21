%% ========================================================================
%  VISUAL METAPHOR — Group ISPC sweep + cluster-based permutation tests
%  ========================================================================
%  For each (band x seed) combination:
%    1) Compute per-subject seed-ISPC (Laplacian, baseline-corrected)
%       and cache to disk so subsequent runs skip the long compute.
%    2) Run cluster-based permutation on the Met-Lit and Met-Ano contrasts
%       (sign-flip permutation, channel x time clusters, 75-mm adjacency).
%    3) Run a directional, hypothesis-driven LPC-window test (800-1500 ms)
%       at a small ROI around the seed.
%    4) Produce one summary figure per (band, seed) combo and one master
%       summary table at the end.
%
%  Defaults reflect the recommendation in the manuscript discussion:
%    Bands:  theta [4 8], alpha [8 12], beta [13 25]
%    Seeds:  C4, FC2  (LPC peak channels from the mass-univariate ERP)
%            + FCz as exploratory comparison.
%  Edit BANDS and SEEDS below to add/remove sweep points.
% =========================================================================

clear; close all; clc;

%% -----------------------------------------------------------------------
%  STEP 0: Paths and parameters
%  -----------------------------------------------------------------------
study_dir  = 'C:\Users\Cagi\Desktop\VisualMet';
set_dir    = fullfile(study_dir, 'EEGSets');
func_dir   = 'C:\Users\Cagi\Desktop\516';
out_dir    = 'C:\Users\Cagi\Desktop\516\Day_12';
cache_dir  = fullfile(out_dir, 'ispc_cache');
fig_dir    = fullfile(out_dir, 'ispc_figs');
if ~exist(cache_dir,'dir'), mkdir(cache_dir); end
if ~exist(fig_dir,'dir'),   mkdir(fig_dir);   end

addpath(func_dir);
if exist('eeglab','file') ~= 2
    error('EEGLAB is not on the MATLAB path. Start EEGLAB first.');
end

% Sweep grid (edit here to subset)
BANDS = { [4 8],  [8 12], [13 25] };
BAND_NAMES = {'theta','alpha','beta'};
SEEDS = {'C4','FC2','FCz'};

% Common analysis params
num_frex        = 5;
time_range      = [-200 1500];
baseline_window = [-200 0];
baseline_flag   = 1;
use_laplacian   = true;

% Permutation params
n_perm          = 1000;
cluster_alpha   = 0.05;     % cluster-forming threshold (two-tailed t)
neighbor_dist   = 75;       % mm, for channel adjacency
lpc_window      = [800 1500];

% Subjects to exclude (failed QC etc.)
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
fprintf('Found %d subjects (after excluding %s).\n', ...
    numel(sub_ids), strjoin(EXCLUDE,','));

%% -----------------------------------------------------------------------
%  STEP 1: Sweep — compute or load each (band, seed) cache
%  -----------------------------------------------------------------------
results = struct([]);   % collected summary per combo

for bi = 1:numel(BANDS)
    for si = 1:numel(SEEDS)
        band_name = BAND_NAMES{bi};
        seed_req  = SEEDS{si};
        tag = sprintf('%s_%s', band_name, seed_req);
        cache_file = fullfile(cache_dir, sprintf('ispc_%s.mat', tag));

        fprintf('\n========== %s (%g-%g Hz) at %s ==========\n', ...
            upper(band_name), BANDS{bi}(1), BANDS{bi}(2), seed_req);

        if exist(cache_file,'file')
            fprintf('Loading cache: %s\n', cache_file);
            S = load(cache_file);
            ispc_group   = S.ispc_group;
            chanlocs     = S.chanlocs;
            times_out    = S.times_out;
            sub_ids_used = S.sub_ids_used;
            seed_used    = S.seed_used;
            cond_names   = S.cond_names;

            % Filter out any excluded subjects that were cached previously
            drop = ismember(sub_ids_used, EXCLUDE);
            if any(drop)
                fprintf('  Dropping %d cached subject(s) per EXCLUDE: %s\n', ...
                    sum(drop), strjoin(sub_ids_used(drop),','));
                ispc_group   = ispc_group(~drop,:,:,:);
                sub_ids_used = sub_ids_used(~drop);
            end
        else
            [ispc_group, chanlocs, times_out, sub_ids_used, seed_used, cond_names] = ...
                local_compute_group_ispc(sub_ids, set_dir, seed_req, ...
                BANDS{bi}, num_frex, time_range, baseline_window, ...
                baseline_flag, use_laplacian);
            save(cache_file, 'ispc_group','chanlocs','times_out', ...
                'sub_ids_used','seed_used','cond_names', '-v7.3');
            fprintf('Cached %s (%d subs).\n', tag, size(ispc_group,1));
        end

        % --- Cluster-based permutation: Met-Lit, Met-Ano ---
        adj = local_chan_adjacency(chanlocs, neighbor_dist);

        fprintf('Cluster perm: Met-Lit (n_perm=%d)...\n', n_perm);
        [t_ML, clus_ML] = local_cluster_perm( ...
            squeeze(ispc_group(:,:,:,2)), squeeze(ispc_group(:,:,:,1)), ...
            adj, cluster_alpha, n_perm);

        fprintf('Cluster perm: Met-Ano (n_perm=%d)...\n', n_perm);
        [t_MA, clus_MA] = local_cluster_perm( ...
            squeeze(ispc_group(:,:,:,2)), squeeze(ispc_group(:,:,:,3)), ...
            adj, cluster_alpha, n_perm);

        % --- ROI test in LPC window (mean ISPC across ROI x window) ---
        % ROI: seed channel + its 75-mm neighbors (small, hypothesis-driven)
        seed_idx = find(strcmpi({chanlocs.labels}, seed_used), 1);
        roi_mask = false(numel(chanlocs),1); roi_mask(seed_idx) = true;
        roi_mask(adj(seed_idx,:)) = true;
        roi_chans = {chanlocs(roi_mask).labels};

        t_idx_lpc = times_out >= lpc_window(1) & times_out <= lpc_window(2);
        roi_lit = local_nanmean(local_nanmean(ispc_group(:,roi_mask,t_idx_lpc,1), 2), 3);
        roi_met = local_nanmean(local_nanmean(ispc_group(:,roi_mask,t_idx_lpc,2), 2), 3);
        roi_ano = local_nanmean(local_nanmean(ispc_group(:,roi_mask,t_idx_lpc,3), 2), 3);
        roi_lit = roi_lit(:); roi_met = roi_met(:); roi_ano = roi_ano(:);

        [tML_roi, pML_roi] = local_paired_t_1d(roi_met, roi_lit);
        [tMA_roi, pMA_roi] = local_paired_t_1d(roi_met, roi_ano);

        % --- Plot summary figure ---
        local_plot_combo(ispc_group, chanlocs, times_out, seed_used, ...
            band_name, BANDS{bi}, t_ML, clus_ML, t_MA, clus_MA, ...
            roi_chans, lpc_window, fig_dir, tag);

        % --- Collect summary row ---
        row.band     = band_name;
        row.band_hz  = BANDS{bi};
        row.seed     = seed_used;
        row.n_subs   = size(ispc_group,1);
        row.ML_min_p = min([clus_ML.p_cluster, NaN]);
        row.ML_n_sig = sum([clus_ML.p_cluster] < cluster_alpha);
        row.MA_min_p = min([clus_MA.p_cluster, NaN]);
        row.MA_n_sig = sum([clus_MA.p_cluster] < cluster_alpha);
        row.roi_chans      = strjoin(roi_chans, ',');
        row.ML_roi_t       = tML_roi;
        row.ML_roi_p       = pML_roi;
        row.ML_roi_meanDiff = mean(roi_met - roi_lit);
        row.MA_roi_t       = tMA_roi;
        row.MA_roi_p       = pMA_roi;
        row.MA_roi_meanDiff = mean(roi_met - roi_ano);
        results = [results, row]; %#ok<AGROW>
    end
end

%% -----------------------------------------------------------------------
%  STEP 2: Summary table
%  -----------------------------------------------------------------------
fprintf('\n\n================ SUMMARY ================\n');
T = struct2table(results);
disp(T(:, {'band','seed','n_subs', ...
          'ML_min_p','ML_n_sig','ML_roi_t','ML_roi_p','ML_roi_meanDiff', ...
          'MA_min_p','MA_n_sig','MA_roi_t','MA_roi_p','MA_roi_meanDiff'}));

writetable(T, fullfile(out_dir, 'ispc_sweep_summary.csv'));
fprintf('Summary written to %s\n', fullfile(out_dir,'ispc_sweep_summary.csv'));
fprintf('Per-combo figures saved in %s\n', fig_dir);


% =========================================================================
%  Compute per-subject ISPC for a given (band, seed)
% =========================================================================
function [ispc_group, chanlocs_out, times_out, sub_ids_used, seed_used, cond_names] = ...
    local_compute_group_ispc(sub_ids, set_dir, seed_req, band, num_frex, ...
    time_range, baseline_window, baseline_flag, use_laplacian)

cond_names  = {'Literal','Metaphor','Anomalous'};
valid_codes = [111 121 131 112 122 132 113 123 133];
n_conds     = 3;

ispc_group   = [];
chanlocs_out = [];
times_out    = [];
sub_ids_used = {};
seed_used    = '';

for s = 1:numel(sub_ids)
    sid = sub_ids{s};
    fname = sprintf('MetVis_%s_postart.set', sid);
    fpath = fullfile(set_dir, fname);
    if ~exist(fpath,'file'), continue; end

    fprintf('  [%2d/%d] sub %s ', s, numel(sub_ids), sid);
    EEG = pop_loadset('filename', fname, 'filepath', set_dir);

    % Resolve seed for this subject (Fz fallback if requested seed missing)
    sc = seed_req;
    if ~any(strcmpi({EEG.chanlocs.labels}, sc))
        if any(strcmpi({EEG.chanlocs.labels}, 'Fz')), sc = 'Fz';
        else, fprintf('-- no seed/Fz, skip\n'); continue; end
    end

    % Classify trials
    cond_of_epoch = nan(1, EEG.trials);
    for ep = 1:EEG.trials
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
        numtypes = numtypes(keep); numlats = numlats(keep);
        [~, mi] = min(abs(numlats));
        switch mod(numtypes(mi),10)
            case 1, cond_of_epoch(ep) = 1;
            case 2, cond_of_epoch(ep) = 2;
            case 3, cond_of_epoch(ep) = 3;
        end
    end
    n_per = [sum(cond_of_epoch==1), sum(cond_of_epoch==2), sum(cond_of_epoch==3)];
    if min(n_per) < 10, fprintf('-- low trial count, skip\n'); continue; end

    % Laplacian
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

    % ISPC per condition
    ispc_subj = [];
    for c = 1:n_conds
        data_c = data_for_ispc(:,:,cond_of_epoch == c);
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
    if size(ispc_subj,1) ~= size(ispc_group,2) || ...
       size(ispc_subj,2) ~= size(ispc_group,3)
        fprintf('-- size mismatch, skip\n'); continue;
    end
    ispc_group(s,:,:,:) = ispc_subj;
    sub_ids_used{end+1} = sid; %#ok<AGROW>
    fprintf('done.\n');
end

% Drop empty slots
keep = squeeze(any(any(any(~isnan(ispc_group),2),3),4));
ispc_group = ispc_group(keep,:,:,:);
end


% =========================================================================
%  Cluster-based permutation test
%   A, B: [n_subs x n_chans x n_times]
%   Returns observed t-stat map and an array of cluster structs:
%     .t_thresh  scalar
%     .clusters  per-cluster struct array
%       .channels   indices in cluster
%       .times      time-point indices
%       .mass       sum of |t| in cluster
%       .sign       +1 or -1
%       .p_cluster  cluster-level p-value
%   Output `clus` is the per-cluster struct array (flattened across signs).
% =========================================================================
function [t_obs, clus] = local_cluster_perm(A, B, adj, alpha, n_perm)

D    = A - B;                        % [n_subs x n_chans x n_times]
n    = sum(~isnan(D),1);
m    = local_nanmean(D, 1);
sd   = local_nanstd(D, 1);
t_obs = squeeze(m ./ (sd ./ sqrt(n)));   % [n_chans x n_times]
df    = squeeze(n) - 1;
df0   = max(df(:));                  % conservative for threshold lookup
t_thresh = abs(tinv(alpha/2, df0));

% Observed clusters (separate by sign)
[obs_pos, obs_neg] = local_find_clusters(t_obs,  t_thresh, adj);
obs_pos_masses = arrayfun(@(c) c.mass, obs_pos);
obs_neg_masses = arrayfun(@(c) c.mass, obs_neg);

% Permutation: sign-flip subjects
n_subs = size(D,1);
max_pos = nan(n_perm,1);
max_neg = nan(n_perm,1);
for p = 1:n_perm
    flip = (rand(n_subs,1) < 0.5) * 2 - 1;       % +1/-1
    Dp   = D .* flip;
    mp   = local_nanmean(Dp,1);
    sdp  = local_nanstd(Dp,1);
    tp   = squeeze(mp ./ (sdp ./ sqrt(n)));
    [pos_p, neg_p] = local_find_clusters(tp, t_thresh, adj);
    if isempty(pos_p), max_pos(p) = 0;
    else, max_pos(p) = max(arrayfun(@(c) c.mass, pos_p));
    end
    if isempty(neg_p), max_neg(p) = 0;
    else, max_neg(p) = max(arrayfun(@(c) c.mass, neg_p));
    end
end

% Compute cluster p-values
clus = struct([]);
for k = 1:numel(obs_pos)
    obs_pos(k).sign = +1;
    obs_pos(k).p_cluster = (1 + sum(max_pos >= obs_pos(k).mass)) / (1 + n_perm);
    clus = [clus, obs_pos(k)]; %#ok<AGROW>
end
for k = 1:numel(obs_neg)
    obs_neg(k).sign = -1;
    obs_neg(k).p_cluster = (1 + sum(max_neg >= obs_neg(k).mass)) / (1 + n_perm);
    clus = [clus, obs_neg(k)]; %#ok<AGROW>
end
end


% =========================================================================
%  Find positive and negative clusters in t-map (threshold |t|>thresh,
%  cluster by spatial adjacency + temporal contiguity)
% =========================================================================
function [pos_clus, neg_clus] = local_find_clusters(t_map, thresh, adj)
% t_map: [n_chans x n_times]
n_chans = size(t_map,1);
n_times = size(t_map,2);

mask_pos = t_map >  thresh;
mask_neg = t_map < -thresh;

pos_clus = local_extract_clusters(mask_pos, t_map, adj);
neg_clus = local_extract_clusters(mask_neg, t_map, adj);

    function clusters = local_extract_clusters(mask, tmap, adj_)
        clusters = struct('channels',{},'times',{},'mass',{});
        unvisited = mask;
        % Iterate over (chan, time) seeds, BFS-flood-fill
        while any(unvisited(:))
            [ci, ti] = find(unvisited, 1, 'first');
            queue = [ci, ti];
            chans_in = []; times_in = []; mass = 0;
            while ~isempty(queue)
                c = queue(1,1); t = queue(1,2);
                queue(1,:) = [];
                if c<1 || c>n_chans || t<1 || t>n_times, continue; end
                if ~unvisited(c,t), continue; end
                unvisited(c,t) = false;
                chans_in(end+1) = c; %#ok<AGROW>
                times_in(end+1) = t; %#ok<AGROW>
                mass = mass + abs(tmap(c,t));
                % temporal neighbors
                queue(end+1,:) = [c, t-1]; %#ok<AGROW>
                queue(end+1,:) = [c, t+1]; %#ok<AGROW>
                % spatial neighbors (channels adjacent to c)
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
end


% =========================================================================
%  Channel adjacency from chanlocs.X/Y/Z, threshold in mm
% =========================================================================
function adj = local_chan_adjacency(chanlocs, max_dist)
n = numel(chanlocs);
xyz = nan(n,3);
for k = 1:n
    if ~isempty(chanlocs(k).X) && ~isempty(chanlocs(k).Y) && ~isempty(chanlocs(k).Z)
        xyz(k,:) = [chanlocs(k).X, chanlocs(k).Y, chanlocs(k).Z];
    end
end
% chanloc XYZ in EEGLAB is unitless; rescale so head radius ~ 90 mm.
% Standard 10-05 chanlocs in dipfit are in mm already; check by max norm.
norms = sqrt(sum(xyz.^2, 2));
median_norm = median(norms(~isnan(norms)));
if median_norm < 5
    xyz = xyz * (90 / median_norm);     % rescale to ~mm if unitless
end
D = nan(n);
for i = 1:n
    for j = 1:n
        D(i,j) = norm(xyz(i,:) - xyz(j,:));
    end
end
adj = (D > 0) & (D <= max_dist);
end


% =========================================================================
%  Plotting per (band, seed) combo
% =========================================================================
function local_plot_combo(ispc_group, chanlocs, times_out, seed_used, ...
    band_name, band_hz, t_ML, clus_ML, t_MA, clus_MA, ...
    roi_chans, lpc_window, fig_dir, tag)

n_subs   = size(ispc_group,1);
n_plots  = 8;
t_idx_p  = round(linspace(1,length(times_out),n_plots));
seed_idx = find(strcmpi({chanlocs.labels}, seed_used),1);
cmap_div = local_diverging_cmap();

% Build sig-channel-per-time mask from clusters
sig_mask_ML = false(numel(chanlocs), length(times_out));
for k = 1:numel(clus_ML)
    if clus_ML(k).p_cluster < 0.05
        for j = 1:numel(clus_ML(k).channels)
            sig_mask_ML(clus_ML(k).channels(j), clus_ML(k).times(j)) = true;
        end
    end
end
sig_mask_MA = false(numel(chanlocs), length(times_out));
for k = 1:numel(clus_MA)
    if clus_MA(k).p_cluster < 0.05
        for j = 1:numel(clus_MA(k).channels)
            sig_mask_MA(clus_MA(k).channels(j), clus_MA(k).times(j)) = true;
        end
    end
end

clim_t = max([abs(t_ML(:)); abs(t_MA(:))]) * 0.9;
if ~isfinite(clim_t) || clim_t == 0, clim_t = 3; end

fig = figure('Color','w','Position',[20 20 1500 700], ...
    'Name', sprintf('Cluster perm — %s @ %s (n=%d)', band_name, seed_used, n_subs));
sgtitle(sprintf(['GROUP n=%d  |  seed=%s  |  %s (%g-%g Hz)  |  ', ...
    'cluster-perm t-maps; black dots = p_{cluster}<0.05'], ...
    n_subs, seed_used, band_name, band_hz(1), band_hz(2)), ...
    'FontSize', 12, 'FontWeight','bold');

row_titles = {'Met − Lit  (t)', 'Met − Ano  (t)'};
row_t      = {t_ML, t_MA};
row_sig    = {sig_mask_ML, sig_mask_MA};

for r = 1:2
    for p = 1:n_plots
        ax = subplot(2, n_plots, (r-1)*n_plots + p);
        ti = t_idx_p(p);
        emk = {seed_idx,'o','r',6,1.5};
        sig_here = find(row_sig{r}(:, ti));
        if ~isempty(sig_here)
            emk = {sig_here,'.','k',12,1};
        end
        topoplot(row_t{r}(:,ti), chanlocs, ...
            'maplimits', [-clim_t clim_t], 'electrodes','on', ...
            'style','both','numcontour',6, 'emarker2', emk);
        colormap(ax, cmap_div);
        if r == 1, title(sprintf('%d ms', round(times_out(ti))), 'FontSize',9); end
        if p == 1, local_row_label(ax, row_titles{r}); end
        if p == n_plots
            cb = colorbar(ax,'eastoutside'); cb.Label.String = 't'; cb.FontSize = 8;
        end
    end
end

% Save
fname = fullfile(fig_dir, sprintf('cluster_%s.png', tag));
exportgraphics(fig, fname, 'Resolution', 150);

% Also save a small ROI time-course figure
fig2 = figure('Color','w','Position',[40 40 900 400], ...
    'Name', sprintf('ROI time course — %s @ %s', band_name, seed_used));
ch_idx = ismember({chanlocs.labels}, roi_chans);
m_lit = squeeze(local_nanmean(local_nanmean(ispc_group(:,ch_idx,:,1),1),2));
m_met = squeeze(local_nanmean(local_nanmean(ispc_group(:,ch_idx,:,2),1),2));
m_ano = squeeze(local_nanmean(local_nanmean(ispc_group(:,ch_idx,:,3),1),2));
plot(times_out, m_lit,'b-', times_out, m_met,'r-', times_out, m_ano,'y-', ...
    'LineWidth', 1.6); hold on;
xline(0,'k--'); yline(0,'k:');
patch([lpc_window(1) lpc_window(2) lpc_window(2) lpc_window(1)], ...
      [min(ylim) min(ylim) max(ylim) max(ylim)], [0.85 0.85 0.85], ...
      'FaceAlpha', 0.25, 'EdgeColor','none');
xlabel('Time (ms)'); ylabel('ISPC (% change)');
title(sprintf('ROI mean (%s) — %s @ %s', strjoin(roi_chans,','), band_name, seed_used));
legend({'Lit','Met','Ano','stim','zero','LPC window'},'Location','best');
exportgraphics(fig2, fullfile(fig_dir, sprintf('roi_%s.png',tag)), 'Resolution',150);
end


% =========================================================================
%  Helpers: ID extraction, NaN-safe stats, paired t, label, colormap
% =========================================================================
function code = local_extract_code(x, valid_codes)
    code = NaN;
    if isnumeric(x) && ~isempty(x)
        if ismember(x, valid_codes), code = x; end
        return;
    end
    if ~ischar(x) || isempty(x), return; end
    s = x; if s(1) == 'S', s = s(2:end); end
    if ~isempty(s) && all(isstrprop(s,'digit'))
        v = str2double(s);
        if ismember(v, valid_codes), code = v; end
        return;
    end
    toks = regexp(x, '\d+', 'match');
    for k = 1:numel(toks)
        v = str2double(toks{k});
        if ismember(v, valid_codes), code = v; return; end
    end
end

function m = local_nanmean(X, dim)
    valid = ~isnan(X);
    Xz = X; Xz(~valid) = 0;
    n = sum(valid, dim); s = sum(Xz, dim);
    m = s ./ n; m(n==0) = NaN;
end

function sd = local_nanstd(X, dim)
    valid = ~isnan(X);
    n = sum(valid, dim);
    mu = local_nanmean(X, dim);
    Xc = X - mu; Xc(~valid) = 0;
    ss = sum(Xc.^2, dim);
    sd = sqrt(ss ./ max(n-1, 1));
    sd(n < 2) = NaN;
end

function [t, p] = local_paired_t_1d(a, b)
    d = a - b; d = d(~isnan(d));
    n = numel(d); m = mean(d); sd = std(d);
    t = m / (sd / sqrt(n));
    p = 2 * (1 - tcdf(abs(t), n-1));
end

function local_row_label(ax, txt)
    pos = get(ax,'Position');
    annotation('textbox', ...
        [0.001, pos(2)+pos(4)/2 - 0.025, pos(1)-0.005, 0.05], ...
        'String', txt, 'FontSize', 10, 'FontWeight','bold', ...
        'EdgeColor','none', 'HorizontalAlignment','right', ...
        'VerticalAlignment','middle');
end

function cmap = local_diverging_cmap()
n = 128;
low  = linspace(0,1,n+1); high = linspace(1,0,n+1);
r = [low(1:end-1), ones(1,n)];
g = [low(1:end-1), high(2:end)];
b = [ones(1,n),    high(2:end)];
cmap = [r(:), g(:), b(:)];
end
