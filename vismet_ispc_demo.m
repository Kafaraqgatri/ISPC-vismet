%% ========================================================================
%  VISUAL METAPHOR — Seed-Based ISPC Topography (single-subject, exploratory)
%  ========================================================================
%  Runs seed_ispc_topo() on one subject's cleaned, epoched set from the
%  VisualMet study for each of the three conditions (Literal, Metaphor,
%  Abstract/Anomalous) plus a Metaphor-minus-Literal difference map.
%
%  Seed:      FCz  (exploratory default — midline frontal, implicated in
%                   theta-band cognitive control / semantic integration;
%                   falls back to Fz if FCz is not in the montage)
%  Data:      Surface Laplacian (CSD) — reduces volume conduction so ISPC
%             reflects genuine inter-site synchronization rather than
%             shared signal spread. Toggle with use_laplacian below.
%  Band:      theta [4 8] Hz
%  Window:    -200 to 1500 ms  (matches preprocessing epoch window)
%  Baseline:  -200 to 0 ms
%
%  Condition event codes (from Code\MetVis_bdf.txt):
%     Literal   : 111, 121, 131
%     Metaphor  : 112, 122, 132
%     Abstract  : 113, 123, 133     (Anomalous in pipeline)
% =========================================================================

clear; close all; clc;

%% -----------------------------------------------------------------------
%  STEP 0: Paths and subject selection
%  -----------------------------------------------------------------------
study_dir  = 'C:\Users\Cagi\Desktop\VisualMet';
set_dir    = fullfile(study_dir, 'EEGSets');
func_dir   = 'C:\Users\Cagi\Desktop\516';     % seed_ispc_topo.m lives here

addpath(func_dir);                            % make the function visible
if exist('eeglab','file') ~= 2
    error(['EEGLAB is not on the MATLAB path. Start EEGLAB once ', ...
           '(or addpath to eeglab2024.0) before running this demo.']);
end

sub_id        = '03';                         % <-- change to try other subjects
use_laplacian = true;                         % true = CSD/Laplacian, false = raw voltage

% Subjects to exclude from any analysis (e.g. failed QC).
EXCLUDE = {'32'};
if any(strcmp(sub_id, EXCLUDE))
    error('Subject %s is in the EXCLUDE list — pick a different sub_id.', sub_id);
end
set_file      = sprintf('MetVis_%s_postart.set', sub_id);

fprintf('Loading %s ...\n', set_file);
EEG = pop_loadset('filename', set_file, 'filepath', set_dir);
fprintf('  %d channels, %d time points, %d epochs\n', ...
    EEG.nbchan, EEG.pnts, EEG.trials);

% Show channel labels so we can pick a seed that actually exists
fprintf('\nChannel labels: %s\n', strjoin({EEG.chanlocs.labels}, ', '));

% Peek at one epoch's events so we can see how ERPLAB stored them
if EEG.trials >= 1
    fprintf('\nSample EEG.epoch(1) fields: %s\n', ...
        strjoin(fieldnames(EEG.epoch(1)), ', '));
    et = EEG.epoch(1).eventtype;
    if ~iscell(et), et = {et}; end
    try
        fprintf('EEG.epoch(1).eventtype (first few): ');
        disp(et(1:min(end,5)));
    catch
    end
end

%% -----------------------------------------------------------------------
%  STEP 1: Identify the condition of each epoch
%  -----------------------------------------------------------------------
%  For each epoch, find the event whose eventlatency is closest to 0
%  (the time-locking stimulus) and read its numeric type. Then map the
%  stimulus code to one of three conditions via its last digit.
%     ...1 = Literal, ...2 = Metaphor, ...3 = Abstract
%  (codes are 111/121/131, 112/122/132, 113/123/133 — see MetVis_bdf.txt)

cond_of_epoch = nan(1, EEG.trials);           % 1=Lit, 2=Met, 3=Abs, NaN=skip
valid_codes   = [111 121 131 112 122 132 113 123 133];

% Helper: extract the first 3-digit stim code (111..133) from any string
%         or return the number directly if numeric. Handles ERPLAB forms
%         like 'B1(111)', 'S111', '111', or numeric 111.
extract_code = @(x) local_extract_code(x, valid_codes);

for ep = 1:EEG.trials
    types = EEG.epoch(ep).eventtype;
    lats  = EEG.epoch(ep).eventlatency;
    if ~iscell(types), types = {types}; end
    if ~iscell(lats),  lats  = {lats};  end
    numlats = cell2mat(lats);

    % Pull a valid stim code from each event (NaN if none)
    numtypes = nan(1, numel(types));
    for k = 1:numel(types)
        numtypes(k) = extract_code(types{k});
    end

    keep = ~isnan(numtypes);
    if ~any(keep), continue; end
    numtypes = numtypes(keep);
    numlats  = numlats(keep);

    [~, mi] = min(abs(numlats));              % event closest to 0
    code = numtypes(mi);

    last_digit = mod(code, 10);
    switch last_digit
        case 1, cond_of_epoch(ep) = 1;        % Literal
        case 2, cond_of_epoch(ep) = 2;        % Metaphor
        case 3, cond_of_epoch(ep) = 3;        % Abstract
    end
end

n_lit = sum(cond_of_epoch == 1);
n_met = sum(cond_of_epoch == 2);
n_abs = sum(cond_of_epoch == 3);
fprintf('\nTrial counts:\n');
fprintf('  Literal  : %d\n', n_lit);
fprintf('  Metaphor : %d\n', n_met);
fprintf('  Abstract : %d\n', n_abs);
fprintf('  Unclassified / dropped: %d\n', sum(isnan(cond_of_epoch)));

if min([n_lit n_met n_abs]) < 5
    warning('At least one condition has <5 trials; ISPC will be unstable.');
end

%% -----------------------------------------------------------------------
%  STEP 1b: Surface Laplacian / CSD transform (recommended for ISPC)
%  -----------------------------------------------------------------------
%  Raw-voltage ISPC is inflated between neighboring channels by volume
%  conduction (one source smeared across the scalp). The surface Laplacian
%  acts as a spatial high-pass filter and suppresses that shared signal, so
%  the remaining ISPC better reflects genuine inter-site phase coupling.

if use_laplacian
    fprintf('\nApplying surface Laplacian (Perrin spherical spline)...\n');
    if exist('laplacian_perrinX','file') ~= 2
        error(['laplacian_perrinX not on path. Add the folder that ', ...
               'contains it (class code / 516 folder) with addpath().']);
    end
    if isempty(EEG.chanlocs(1).X) || isempty(EEG.chanlocs(1).Y)
        error('chanlocs missing XYZ — re-run pop_chanedit with lookup.');
    end
    X = [EEG.chanlocs.X];
    Y = [EEG.chanlocs.Y];
    Z = [EEG.chanlocs.Z];

    data_for_ispc = zeros(size(EEG.data));
    for tri = 1:EEG.trials
        [surf_lap, ~, ~] = laplacian_perrinX(EEG.data(:,:,tri), X, Y, Z);
        data_for_ispc(:,:,tri) = surf_lap;
    end
    fprintf('  Laplacian transform complete.\n');
else
    fprintf('\nUsing raw voltage data (no Laplacian).\n');
    data_for_ispc = EEG.data;
end

%% -----------------------------------------------------------------------
%  STEP 2: Run seed_ispc_topo on each condition
%  -----------------------------------------------------------------------
% Pick a seed. FCz is ideal but may not exist on a 32-ch montage;
% fall back to Fz (also a midline frontal site) if so.
seed_chan = 'FCz';
if ~any(strcmpi({EEG.chanlocs.labels}, seed_chan))
    seed_chan = 'Fz';
    fprintf('\nNote: FCz not in montage — using Fz as seed instead.\n');
end
freq_band       = [4 8];                      % theta
num_frex        = 5;
time_range      = [-200 1500];
baseline_window = [-200 0];
baseline_flag   = 1;                          % 1 = % change from baseline,
                                              % 0 = raw ISPC in [0,1].
                                              % Baseline correction is
                                              % strongly recommended here —
                                              % single-subject ISPC has a
                                              % high floor (~1/sqrt(N))
                                              % that hides task-locked
                                              % changes otherwise.
plot_flag       = 0;                          % we plot side-by-side ourselves
save_flag       = 0;

% Sanity check seed exists
if ~any(strcmpi({EEG.chanlocs.labels}, seed_chan))
    error('Seed channel %s not found in EEG.chanlocs.', seed_chan);
end

cond_names = {'Literal','Metaphor','Abstract'};
ispc_by_cond = cell(1,3);
pli_by_cond  = cell(1,3);

for c = 1:3
    fprintf('\n===== %s =====\n', cond_names{c});
    data_c = data_for_ispc(:, :, cond_of_epoch == c);

    [ispc_c, frex, times_out, ~, pli_c] = seed_ispc_topo( ...
        data_c, EEG.srate, EEG.times, EEG.chanlocs, seed_chan, ...
        freq_band, num_frex, time_range, baseline_window, ...
        baseline_flag, plot_flag, save_flag);

    ispc_by_cond{c} = ispc_c;
    pli_by_cond{c}  = pli_c;
end

%% -----------------------------------------------------------------------
%  STEP 3: Side-by-side topography plot (one row per condition + diff)
%  -----------------------------------------------------------------------
n_plots   = 8;
t_indices = round(linspace(1, length(times_out), n_plots));
seed_idx  = find(strcmpi({EEG.chanlocs.labels}, seed_chan), 1);

data_tag = ternary(use_laplacian, 'Laplacian/CSD', 'Raw voltage');
figure('Color','w','Position',[30 30 1500 800], ...
    'Name', sprintf('Sub %s — Seed ISPC (%s, theta, %s)', ...
    sub_id, seed_chan, data_tag));
bl_tag = ternary(baseline_flag, ...
    sprintf('%% change vs. [%g %g] ms baseline', baseline_window(1), baseline_window(2)), ...
    'raw ISPC');
sgtitle(sprintf('Subject %s  |  seed = %s  |  theta %g–%g Hz  |  %s  |  %s', ...
    sub_id, seed_chan, freq_band(1), freq_band(2), data_tag, bl_tag), ...
    'FontSize', 13, 'FontWeight', 'bold');

n_rows = 4;                                   % 3 conds + (Met - Lit)
row_labels = {cond_names{:}, 'Metaphor − Literal'};
row_data   = {ispc_by_cond{1}, ispc_by_cond{2}, ispc_by_cond{3}, ...
              ispc_by_cond{2} - ispc_by_cond{1}};

% Color-limit policy:
%   baseline_flag == 1 : symmetric diverging per row (% change)
%   baseline_flag == 0 : [0,1] parula for raw condition rows,
%                        symmetric diverging for the diff row
for r = 1:n_rows
    is_diff = (r == 4);
    if baseline_flag || is_diff
        mx = max(abs(row_data{r}(:))) * 0.9;
        if ~isfinite(mx) || mx == 0, mx = 0.01; end
        clims  = [-mx mx];
        cmap_r = make_diverging_cmap_local();
    else
        clims  = [0 1];
        cmap_r = parula(256);
    end

    for p = 1:n_plots
        ax = subplot(n_rows, n_plots, (r-1)*n_plots + p);
        ti = t_indices(p);
        topoplot(row_data{r}(:, ti), EEG.chanlocs, ...
            'maplimits',  clims, ...
            'electrodes', 'on', ...
            'style',      'both', ...
            'numcontour', 6, ...
            'emarker2',   {seed_idx, 'o', 'r', 6, 1.5});
        colormap(ax, cmap_r);

        if r == 1
            title(sprintf('%d ms', round(times_out(ti))), 'FontSize', 9);
        end

        % Row label via text() on the leftmost subplot — ylabel() is
        % swallowed by topoplot, so we drop text in normalized axes.
        if p == 1
            pos = get(ax, 'Position');
            annotation('textbox', ...
                [0.001, pos(2)+pos(4)/2 - 0.02, pos(1)-0.005, 0.04], ...
                'String', row_labels{r}, ...
                'FontSize', 11, 'FontWeight', 'bold', ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'middle');
        end

        % Per-row colorbar on the rightmost subplot (so the reader can
        % read off % change rather than just relative color).
        if p == n_plots
            cb = colorbar(ax, 'eastoutside');
            cb.FontSize = 8;
            if baseline_flag || is_diff
                cb.Label.String = '% change';
            else
                cb.Label.String = 'ISPC';
            end
        end
    end
end

%% -----------------------------------------------------------------------
%  STEP 4: Time-course at a few example channels (FCz, Cz, Pz, Oz, F7, F8)
%  -----------------------------------------------------------------------
probe_chans = {'Cz','Pz','Oz','F7','F8'};
colors = lines(3);

figure('Color','w','Position',[50 50 1200 700], ...
    'Name', sprintf('Sub %s — ISPC time course (vs %s, theta, %s)', ...
    sub_id, seed_chan, data_tag));

for i = 1:length(probe_chans)
    subplot(2, 3, i); hold on;
    ch_idx = find(strcmpi(probe_chans{i}, {EEG.chanlocs.labels}), 1);
    if isempty(ch_idx)
        title(sprintf('%s (not found)', probe_chans{i})); continue;
    end
    for c = 1:3
        plot(times_out, ispc_by_cond{c}(ch_idx, :), '-', ...
            'Color', colors(c,:), 'LineWidth', 1.6);
    end
    xline(0, 'k--', 'Stim');
    xlabel('Time (ms)');
    if baseline_flag
        ylabel('ISPC (% change)');
    else
        ylabel('ISPC');
        ylim([0 1]);
    end
    title(sprintf('%s  (seed=%s)', probe_chans{i}, seed_chan));
    set(gca, 'FontSize', 10);
    if i == 1, legend(cond_names, 'Location','northwest'); end
end

fprintf('\n=== Done. ===\n');
fprintf('Tip: change sub_id at the top to inspect other subjects.\n');


% =========================================================================
%  Local diverging colormap (same design as seed_ispc_topo's internal one,
%  duplicated here so this demo has no extra file dependencies)
% =========================================================================
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function cmap = make_diverging_cmap_local()
n = 128;
low  = linspace(0, 1, n+1);
high = linspace(1, 0, n+1);
r = [low(1:end-1),  ones(1, n)];
g = [low(1:end-1),  high(2:end)];
b = [ones(1, n),    high(2:end)];
cmap = [r(:), g(:), b(:)];
end

% -------------------------------------------------------------------------
%  extract_code: pull a valid stim code (e.g. 111, 132) out of whatever
%  ERPLAB or EEGLAB put in EEG.epoch(ep).eventtype. Handles:
%     numeric:    111                  -> 111
%     'S111'                           -> 111
%     '111'                            -> 111
%     'B1(111)'  (ERPLAB bin-tagged)   -> 111
%     'boundary' or anything else      -> NaN
%  Only returns codes that are in the valid_codes whitelist so we don't
%  accidentally match an unrelated 3-digit token.
% -------------------------------------------------------------------------
function code = local_extract_code(x, valid_codes)
    code = NaN;
    if isnumeric(x) && ~isempty(x)
        if ismember(x, valid_codes), code = x; end
        return;
    end
    if ~ischar(x) || isempty(x), return; end

    % Fast path: pure-digit string (optionally with leading 'S')
    s = x;
    if s(1) == 'S', s = s(2:end); end
    if ~isempty(s) && all(isstrprop(s, 'digit'))
        v = str2double(s);
        if ismember(v, valid_codes), code = v; end
        return;
    end

    % General path: scan for any 3-digit number in the string (e.g. 'B1(111)')
    toks = regexp(x, '\d+', 'match');
    for k = 1:numel(toks)
        v = str2double(toks{k});
        if ismember(v, valid_codes)
            code = v;
            return;
        end
    end
end
