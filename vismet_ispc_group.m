%% ========================================================================
%  VISUAL METAPHOR — Group-Level Seed-Based ISPC
%  ========================================================================
%  Computes baseline-corrected, Laplacian-transformed seed ISPC for every
%  subject in EEGSets, averages across subjects, and plots:
%     1) Per-condition group mean (Literal / Metaphor / Abstract).
%     2) Metaphor − Literal difference topography over time.
%     3) Metaphor − Abstract difference topography over time.
%     4) Time courses at probe channels (mean ± SE across subjects).
%
%  Compute caches to a .mat file so re-running just replots.
%
%  Seed:      FCz if present, else Fz (midline frontal).
%  Band:      theta [4 8] Hz, 5 log-spaced frequencies.
%  Window:    -200 to 1500 ms (matches preprocessing epoch window).
%  Baseline:  -200 to 0 ms, percent-change correction.
%  Data:      Surface Laplacian (CSD) — recommended for ISPC.
% =========================================================================

clear; close all; clc;

%% -----------------------------------------------------------------------
%  STEP 0: Paths, parameters, cache file
%  -----------------------------------------------------------------------
study_dir = 'C:\Users\Cagi\Desktop\VisualMet';
set_dir   = fullfile(study_dir, 'EEGSets');
func_dir  = 'C:\Users\Cagi\Desktop\516';
out_dir   = 'C:\Users\Cagi\Desktop\516\Day_12';
cache_file = fullfile(out_dir, 'vismet_ispc_group_cache.mat');

addpath(func_dir);                            % seed_ispc_topo lives here

if exist('eeglab','file') ~= 2
    error('EEGLAB is not on the MATLAB path. Start EEGLAB first.');
end

% Analysis parameters (edit here, then delete the cache to recompute)
seed_chan       = 'FCz';                      % falls back to Fz if missing
freq_band       = [4 8];                      % theta
num_frex        = 5;
time_range      = [-200 1500];
baseline_window = [-200 0];
baseline_flag   = 1;                          % 1 = % change (cross-subject)
use_laplacian   = true;
force_recompute = false;                      % set true to ignore cache

% Subjects to exclude from any analysis (failed QC, etc.).
EXCLUDE = {'32'};

% Subjects: glob whatever postart files exist in EEGSets
files   = dir(fullfile(set_dir, 'MetVis_*_postart.set'));
sub_ids = cell(1, numel(files));
for k = 1:numel(files)
    tok = regexp(files(k).name, 'MetVis_(\d+)_postart\.set', 'tokens', 'once');
    if ~isempty(tok), sub_ids{k} = tok{1}; end
end
sub_ids = sub_ids(~cellfun(@isempty, sub_ids));
sub_ids = sub_ids(~ismember(sub_ids, EXCLUDE));
fprintf('Found %d subjects (after excluding %s): %s\n', ...
    numel(sub_ids), strjoin(EXCLUDE,','), strjoin(sub_ids, ', '));

%% -----------------------------------------------------------------------
%  STEP 1: Compute (or load cached) per-subject ISPC
%  -----------------------------------------------------------------------
%  Result: ispc_group is [n_subs x n_chans x n_times x n_conds]
%          conds 1=Lit, 2=Met, 3=Abs

if exist(cache_file, 'file') && ~force_recompute
    fprintf('\nLoading cached results from %s\n', cache_file);
    S = load(cache_file);
    ispc_group   = S.ispc_group;
    times_out    = S.times_out;
    chanlocs     = S.chanlocs;
    seed_chan    = S.seed_chan;
    sub_ids_used = S.sub_ids_used;
    cond_names   = S.cond_names;
    n_trials_mat = S.n_trials_mat;

    % Trim cached data to drop any excluded subjects on the fly so we
    % don't have to recompute when the EXCLUDE list changes.
    drop = ismember(sub_ids_used, EXCLUDE);
    if any(drop)
        fprintf('  Dropping %d cached subjects per EXCLUDE: %s\n', ...
            sum(drop), strjoin(sub_ids_used(drop), ','));
        ispc_group   = ispc_group(~drop,:,:,:);
        sub_ids_used = sub_ids_used(~drop);
        n_trials_mat = n_trials_mat(~drop,:);
    end
    fprintf('  Using %d subjects, %d channels, %d timepoints, %d conditions.\n', ...
        size(ispc_group,1), size(ispc_group,2), size(ispc_group,3), size(ispc_group,4));
else
    cond_names   = {'Literal','Metaphor','Abstract'};
    n_conds      = 3;
    valid_codes  = [111 121 131 112 122 132 113 123 133];

    ispc_group   = [];
    n_trials_mat = nan(numel(sub_ids), n_conds);
    sub_ids_used = {};
    chanlocs     = [];
    times_out    = [];

    for s = 1:numel(sub_ids)
        sid = sub_ids{s};
        fname = sprintf('MetVis_%s_postart.set', sid);
        fpath = fullfile(set_dir, fname);
        if ~exist(fpath, 'file')
            fprintf('  [skip] %s not found.\n', fname); continue;
        end

        fprintf('\n----- Subject %s (%d/%d) -----\n', sid, s, numel(sub_ids));
        EEG = pop_loadset('filename', fname, 'filepath', set_dir);

        % --- Resolve seed channel for this subject ---
        sc = seed_chan;
        if ~any(strcmpi({EEG.chanlocs.labels}, sc))
            if any(strcmpi({EEG.chanlocs.labels}, 'Fz'))
                sc = 'Fz';
            else
                fprintf('  [skip] No FCz or Fz on this montage.\n'); continue;
            end
        end

        % --- Classify epochs by condition (last digit of the stim code) ---
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
            numtypes = numtypes(keep);
            numlats  = numlats(keep);
            [~, mi]  = min(abs(numlats));
            switch mod(numtypes(mi), 10)
                case 1, cond_of_epoch(ep) = 1;
                case 2, cond_of_epoch(ep) = 2;
                case 3, cond_of_epoch(ep) = 3;
            end
        end
        n_trials_mat(s,:) = [sum(cond_of_epoch==1), ...
                             sum(cond_of_epoch==2), sum(cond_of_epoch==3)];
        fprintf('  Trials L/M/A: %d / %d / %d\n', n_trials_mat(s,:));
        if min(n_trials_mat(s,:)) < 10
            fprintf('  [skip] one condition has <10 trials.\n'); continue;
        end

        % --- Laplacian transform (once per subject) ---
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

        % --- Run seed_ispc_topo per condition ---
        ispc_subj = [];   % [n_chans x n_times x n_conds]
        for c = 1:n_conds
            data_c = data_for_ispc(:, :, cond_of_epoch == c);
            [ispc_c, ~, t_out] = seed_ispc_topo( ...
                data_c, EEG.srate, EEG.times, EEG.chanlocs, sc, ...
                freq_band, num_frex, time_range, baseline_window, ...
                baseline_flag, 0, 0);
            if isempty(ispc_subj)
                ispc_subj = nan(size(ispc_c,1), size(ispc_c,2), n_conds);
            end
            ispc_subj(:,:,c) = ispc_c;
        end

        % --- Stash in group array ---
        if isempty(ispc_group)
            ispc_group = nan(numel(sub_ids), size(ispc_subj,1), ...
                             size(ispc_subj,2), n_conds);
            chanlocs   = EEG.chanlocs;        % use first subject's locations
            times_out  = t_out;
        end
        % Sanity: same channel count and times across subjects
        if size(ispc_subj,1) ~= size(ispc_group,2) || ...
           size(ispc_subj,2) ~= size(ispc_group,3)
            fprintf('  [skip] channel/time mismatch with first subject.\n');
            continue;
        end
        ispc_group(s,:,:,:) = ispc_subj;
        sub_ids_used{end+1} = sid;             %#ok<SAGROW>
    end

    % Drop empty subject slots
    keep = squeeze(any(any(any(~isnan(ispc_group),2),3),4));
    ispc_group   = ispc_group(keep,:,:,:);
    n_trials_mat = n_trials_mat(keep,:);

    save(cache_file, 'ispc_group', 'times_out', 'chanlocs', 'seed_chan', ...
        'sub_ids_used', 'cond_names', 'n_trials_mat', '-v7.3');
    fprintf('\nCached %d subjects to %s\n', size(ispc_group,1), cache_file);
end

n_subs = size(ispc_group,1);
fprintf('\n=== %d subjects in group analysis ===\n', n_subs);
mt = local_nanmean(n_trials_mat, 1);
fprintf('Mean trial counts L/M/A: %.1f / %.1f / %.1f\n', mt(1), mt(2), mt(3));

%% -----------------------------------------------------------------------
%  STEP 2: Group-mean topographies per condition (sanity check row)
%  -----------------------------------------------------------------------
seed_idx = find(strcmpi({chanlocs.labels}, seed_chan), 1);
n_plots  = 8;
t_idx_plot = round(linspace(1, length(times_out), n_plots));

% Group means: [n_chans x n_times x n_conds]
group_mean = squeeze(local_nanmean(ispc_group, 1));

% Differences (per subject -> mean -> SE)
diff_ML = squeeze(ispc_group(:,:,:,2) - ispc_group(:,:,:,1));   % Met - Lit
diff_MA = squeeze(ispc_group(:,:,:,2) - ispc_group(:,:,:,3));   % Met - Abs

mean_ML = squeeze(local_nanmean(diff_ML, 1));
mean_MA = squeeze(local_nanmean(diff_MA, 1));

% Paired t-stats (across subjects) for the two contrasts
[t_ML, p_ML] = local_paired_t(ispc_group(:,:,:,2), ispc_group(:,:,:,1));
[t_MA, p_MA] = local_paired_t(ispc_group(:,:,:,2), ispc_group(:,:,:,3));

%% --- Figure 1: per-condition group mean (3 rows × 8 timepoints) ---
fig1 = figure('Color','w','Position',[20 20 1500 700], ...
    'Name', sprintf('GROUP (n=%d) — per-condition mean', n_subs));
sgtitle(sprintf(['GROUP (n=%d)  |  seed=%s  |  theta %g–%g Hz  |  ' ...
    'Laplacian/CSD  |  %% change vs. [%g %g] ms baseline'], ...
    n_subs, seed_chan, freq_band(1), freq_band(2), ...
    baseline_window(1), baseline_window(2)), ...
    'FontSize', 13, 'FontWeight','bold');

cmap_div = local_diverging_cmap();
clim_cond = max(abs(group_mean(:))) * 0.9;
if clim_cond == 0, clim_cond = 1; end

for r = 1:3
    for p = 1:n_plots
        ax = subplot(3, n_plots, (r-1)*n_plots + p);
        ti = t_idx_plot(p);
        topoplot(group_mean(:,ti,r), chanlocs, ...
            'maplimits', [-clim_cond clim_cond], ...
            'electrodes','on', 'style','both', 'numcontour',6, ...
            'emarker2', {seed_idx,'o','r',6,1.5});
        colormap(ax, cmap_div);
        if r == 1, title(sprintf('%d ms', round(times_out(ti))), 'FontSize',9); end
        if p == 1, local_row_label(ax, cond_names{r}); end
        if p == n_plots
            cb = colorbar(ax,'eastoutside'); cb.Label.String = '% change'; cb.FontSize = 8;
        end
    end
end

%% --- Figure 2: Met-Lit and Met-Abs difference topographies + t-maps ---
fig2 = figure('Color','w','Position',[40 40 1500 700], ...
    'Name', sprintf('GROUP (n=%d) — Met-Lit & Met-Abs differences', n_subs));
sgtitle(sprintf(['GROUP (n=%d)  |  Metaphor − Literal  &  Metaphor − Abstract' ...
    '  |  seed=%s, theta'], n_subs, seed_chan), ...
    'FontSize', 13, 'FontWeight','bold');

clim_diff = max([abs(mean_ML(:)); abs(mean_MA(:))]) * 0.9;
if clim_diff == 0, clim_diff = 1; end
clim_t    = max([abs(t_ML(:)); abs(t_MA(:))]) * 0.9;
if clim_t == 0, clim_t = 1; end
t_crit    = abs(tinv(0.025, n_subs-1));    % two-tailed p<0.05 (uncorrected)

row_titles = {'Met − Lit  (mean % change)', ...
              'Met − Abs  (mean % change)', ...
              sprintf('Met − Lit  (t-stat, |t|>%.2f marked)', t_crit), ...
              sprintf('Met − Abs  (t-stat, |t|>%.2f marked)', t_crit)};
row_data  = {mean_ML, mean_MA, t_ML, t_MA};
row_clims = {[-clim_diff clim_diff], [-clim_diff clim_diff], ...
             [-clim_t clim_t],       [-clim_t clim_t]};
row_units = {'% change','% change','t','t'};

n_rows2 = 4;
for r = 1:n_rows2
    for p = 1:n_plots
        ax = subplot(n_rows2, n_plots, (r-1)*n_plots + p);
        ti = t_idx_plot(p);
        vals = row_data{r}(:, ti);

        % For t-stat rows, mark significant channels with a black dot
        if r >= 3
            sig_mask = abs(vals) > t_crit;
            sig_idx  = find(sig_mask);
            emk = {seed_idx,'o','r',6,1.5};
            if ~isempty(sig_idx)
                emk = {sig_idx,'.','k',10,1};
            end
            topoplot(vals, chanlocs, ...
                'maplimits', row_clims{r}, 'electrodes','on', ...
                'style','both', 'numcontour',6, ...
                'emarker2', emk);
            % overlay the seed marker on top
            hold on;
        else
            topoplot(vals, chanlocs, ...
                'maplimits', row_clims{r}, 'electrodes','on', ...
                'style','both', 'numcontour',6, ...
                'emarker2', {seed_idx,'o','r',6,1.5});
        end
        colormap(ax, cmap_div);
        if r == 1, title(sprintf('%d ms', round(times_out(ti))), 'FontSize',9); end
        if p == 1, local_row_label(ax, row_titles{r}); end
        if p == n_plots
            cb = colorbar(ax,'eastoutside');
            cb.Label.String = row_units{r};
            cb.FontSize = 8;
        end
    end
end

%% --- Figure 3: time courses at probe channels (mean +/- SE across subjects) ---
probe_chans = {'Cz','Pz','Oz','F7','F8','P3','P4','C3','C4'};
colors = lines(3);

fig3 = figure('Color','w','Position',[60 60 1400 800], ...
    'Name', sprintf('GROUP (n=%d) — ISPC time courses', n_subs));
sgtitle(sprintf('GROUP (n=%d) — ISPC time course (vs %s, theta), mean ± SE', ...
    n_subs, seed_chan), 'FontSize', 12, 'FontWeight','bold');

n_p = numel(probe_chans);
ncols = 3; nrows = ceil(n_p/ncols);
for i = 1:n_p
    subplot(nrows, ncols, i); hold on;
    ch_idx = find(strcmpi(probe_chans{i}, {chanlocs.labels}), 1);
    if isempty(ch_idx)
        title(sprintf('%s (n/a)', probe_chans{i})); continue;
    end
    for c = 1:3
        Y = squeeze(ispc_group(:, ch_idx, :, c));   % [n_subs x n_times]
        m = local_nanmean(Y, 1);
        se = local_nanstd(Y, 1) ./ sqrt(sum(~isnan(Y(:,1))));
        % shaded SE band
        x = times_out(:);
        fill([x; flipud(x)], [m(:)-se(:); flipud(m(:)+se(:))], ...
            colors(c,:), 'FaceAlpha', 0.18, 'EdgeColor','none');
        plot(x, m, '-', 'Color', colors(c,:), 'LineWidth', 1.6);
    end
    yl = ylim; xline(0,'k--','Stim','LabelVerticalAlignment','bottom');
    yline(0,'k:');
    xlabel('Time (ms)'); ylabel('ISPC (% change)');
    title(sprintf('%s (seed=%s)', probe_chans{i}, seed_chan));
    set(gca,'FontSize',10);
    if i == 1
        legend({'Lit SE','Lit mean','Met SE','Met mean','Abs SE','Abs mean'}, ...
            'Location','best','FontSize',7);
    end
end

fprintf('\n=== Group analysis complete. ===\n');
fprintf('To recompute from scratch, set force_recompute = true (or delete %s)\n', ...
    cache_file);


% =========================================================================
%  Helpers
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

function [t, p] = local_paired_t(A, B)
% A, B: [n_subs x n_chans x n_times]. Paired t across dim 1.
    D  = A - B;                              % [n_subs x n_chans x n_times]
    n  = sum(~isnan(D), 1);
    m  = local_nanmean(D, 1);
    sd = local_nanstd(D, 1);
    t  = squeeze(m ./ (sd ./ sqrt(n)));      % [n_chans x n_times]
    df = squeeze(n) - 1;
    p  = 2 * (1 - tcdf(abs(t), df));
end

function m = local_nanmean(X, dim)
% NaN-safe mean along dimension `dim`. Works on any MATLAB version
% (avoids the 'omitnan' keyword that older mean() doesn't accept).
    valid = ~isnan(X);
    Xz    = X;
    Xz(~valid) = 0;
    n  = sum(valid, dim);
    s  = sum(Xz, dim);
    m  = s ./ n;
    m(n == 0) = NaN;
end

function sd = local_nanstd(X, dim)
% NaN-safe std along dimension `dim`, using N-1 normalization.
    valid = ~isnan(X);
    n     = sum(valid, dim);
    mu    = local_nanmean(X, dim);
    Xc    = X - mu;
    Xc(~valid) = 0;
    ss    = sum(Xc.^2, dim);
    sd    = sqrt(ss ./ max(n - 1, 1));
    sd(n < 2) = NaN;
end

function local_row_label(ax, txt)
    pos = get(ax, 'Position');
    annotation('textbox', ...
        [0.001, pos(2)+pos(4)/2 - 0.025, pos(1)-0.005, 0.05], ...
        'String', txt, 'FontSize', 10, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'HorizontalAlignment','right', ...
        'VerticalAlignment','middle');
end

function cmap = local_diverging_cmap()
n = 128;
low  = linspace(0, 1, n+1); high = linspace(1, 0, n+1);
r = [low(1:end-1), ones(1,n)];
g = [low(1:end-1), high(2:end)];
b = [ones(1,n),    high(2:end)];
cmap = [r(:), g(:), b(:)];
end
