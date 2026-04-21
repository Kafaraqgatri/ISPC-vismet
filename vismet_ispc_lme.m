%% ========================================================================
%  VISMET ISPC — Single-trial windowed PLV at FC2-ROI in LPC window,
%  fit LME with Condition + Familiarity (parallel to ERP paper §3.2 LME).
%  ========================================================================
%
%  For each trial:
%    1) wavelet-extract theta phase (5 log-spaced freqs 4-8 Hz) at FC2
%       and at every channel in the FC2 75-mm ROI, on Laplacian data
%    2) compute phase difference (FC2 minus ROI channel)
%    3) within-trial PLV = | mean_{time x freq x roi-chan} of exp(i*phase_diff) |
%       across the LPC window (800-1500 ms)
%    -> single value per trial in [0,1]
%
%  Merge with item_name / familiarity / match_status from
%  MetVis_SingleTrial_FINAL.csv, write new SingleTrial_ISPC_LPC.csv,
%  then fit:
%    lme1: ISPC_LPC ~ Condition + Familiarity + (1|Sub) + (1|Item)   (matched only)
%    lme2: same on full set
%    lme3: + Condition*Familiarity interaction
% =========================================================================

clear; close all; clc;

study_dir   = 'C:\Users\Cagi\Desktop\VisualMet';
set_dir     = fullfile(study_dir, 'EEGSets');
behav_csv   = fullfile(study_dir, 'MetVis_SingleTrial_FINAL.csv');
func_dir    = 'C:\Users\Cagi\Desktop\516';
out_dir     = 'C:\Users\Cagi\Desktop\516\Day_12';
out_csv     = fullfile(out_dir, 'SingleTrial_ISPC_LPC.csv');
addpath(func_dir);
if exist('eeglab','file')~=2, error('EEGLAB not on path.'); end

EXCLUDE   = {'32'};
seed_chan = 'FC2';
band      = [4 8];
num_frex  = 5;
range_cycles = [3 10];
lpc_window   = [800 1500];
neighbor_dist = 75;
use_laplacian = true;

% Read behavioral file (per trial: sub, orig_epoch, item, familiarity, match_status, condition)
fprintf('Reading %s\n', behav_csv);
B = readtable(behav_csv);
B.sub = strtrim(string(B.sub));
% Strip leading zeros / make numeric so it joins cleanly
B.sub_num = str2double(B.sub);

files = dir(fullfile(set_dir,'MetVis_*_postart.set'));
sub_ids = cell(1, numel(files));
for k = 1:numel(files)
    tok = regexp(files(k).name,'MetVis_(\d+)_postart\.set','tokens','once');
    if ~isempty(tok), sub_ids{k} = tok{1}; end
end
sub_ids = sub_ids(~cellfun(@isempty, sub_ids));
sub_ids = sub_ids(~ismember(sub_ids, EXCLUDE));
fprintf('Subjects: %d\n', numel(sub_ids));

% --- Build the per-trial table of within-trial PLV ---
out_rows = struct('sub',{},'sub_num',{},'orig_epoch',{},'condition',{}, ...
    'item_name',{},'familiarity',{},'match_status',{},'isPC_LPC',{});

for s = 1:numel(sub_ids)
    sid = sub_ids{s};
    sub_num = str2double(sid);
    fname = sprintf('MetVis_%s_postart.set', sid);
    if ~exist(fullfile(set_dir,fname),'file'), continue; end
    fprintf('\n[%2d/%d] sub %s ', s, numel(sub_ids), sid);
    EEG = pop_loadset('filename', fname, 'filepath', set_dir);

    % Resolve seed (FC2 should exist on all subjects)
    sc = seed_chan;
    if ~any(strcmpi({EEG.chanlocs.labels}, sc))
        fprintf('-- no FC2, skip\n'); continue;
    end
    seed_idx = find(strcmpi({EEG.chanlocs.labels}, sc), 1);

    % Get ROI = seed + neighbors within 75 mm
    adj = local_chan_adjacency(EEG.chanlocs, neighbor_dist);
    roi_idx = find(adj(seed_idx,:) | (1:numel(EEG.chanlocs))==seed_idx);
    roi_idx = setdiff(roi_idx, seed_idx);  % we'll compute phase-diff to non-seed ROI chans
    if isempty(roi_idx), fprintf('-- empty ROI, skip\n'); continue; end

    % Laplacian
    if use_laplacian
        X = [EEG.chanlocs.X]; Y = [EEG.chanlocs.Y]; Z = [EEG.chanlocs.Z];
        data_for_phase = zeros(size(EEG.data));
        for tri = 1:EEG.trials
            [sl,~,~] = laplacian_perrinX(EEG.data(:,:,tri), X, Y, Z);
            data_for_phase(:,:,tri) = sl;
        end
    else
        data_for_phase = EEG.data;
    end

    % Wavelet phases per frequency, per trial: [n_chans x n_pnts x n_trials]
    srate = EEG.srate; times_ms = EEG.times;
    n_pnts = EEG.pnts; n_trials = EEG.trials;
    frex = logspace(log10(band(1)), log10(band(2)), num_frex);
    n_cycles = logspace(log10(range_cycles(1)), log10(range_cycles(2)), num_frex);
    wavelet_time = -2 : 1/srate : 2;
    n_wav = numel(wavelet_time); half_wave = (n_wav-1)/2;
    nData = n_pnts * n_trials;
    nConv = n_wav + nData - 1;

    % We need phases at seed and ROI chans only (no need for whole scalp)
    chans_needed = unique([seed_idx roi_idx]);
    phase_all = nan(numel(chans_needed), num_frex, n_pnts, n_trials);
    for fi = 1:num_frex
        s_g = n_cycles(fi) / (2*pi*frex(fi));
        cmw = exp(2*1i*pi*frex(fi).*wavelet_time) .* exp(-wavelet_time.^2./(2*s_g^2));
        cmw = cmw ./ max(abs(cmw));
        cmw_fft = fft(cmw, nConv);
        for ci = 1:numel(chans_needed)
            ch = chans_needed(ci);
            data_cat = reshape(data_for_phase(ch,:,:), 1, nData);
            as = ifft(fft(data_cat, nConv) .* cmw_fft, nConv);
            as = as(half_wave+1 : end-half_wave);
            phase_all(ci, fi, :, :) = angle(reshape(as, n_pnts, n_trials));
        end
    end

    seed_pos = find(chans_needed == seed_idx, 1);
    roi_pos  = find(chans_needed ~= seed_idx);

    % Time mask for LPC window
    t_mask = times_ms >= lpc_window(1) & times_ms <= lpc_window(2);

    % Per-trial within-trial PLV across (roi-chan, freq, time-in-LPC)
    %   phase_diff(c,f,t,tr) = seed_phase(f,t,tr) - roi_phase(c,f,t,tr)
    %   PLV_tr = | mean over (c,f,t in LPC) of exp(i*phase_diff(c,f,t,tr)) |
    isPC_per_trial = nan(n_trials,1);
    sp_seed = squeeze(phase_all(seed_pos, :, t_mask, :));     % [freq x t_lpc x n_trials]
    sp_roi  = phase_all(roi_pos, :, t_mask, :);                % [n_roi x freq x t_lpc x n_trials]
    for tr = 1:n_trials
        seed_tr = sp_seed(:,:,tr);                             % [freq x t_lpc]
        roi_tr  = squeeze(sp_roi(:,:,:,tr));                   % [n_roi x freq x t_lpc]
        % Broadcast subtract: seed (1 x freq x t_lpc) minus roi (n_roi x freq x t_lpc)
        seed_tr_b = reshape(seed_tr, [1 size(seed_tr,1) size(seed_tr,2)]);
        pd = seed_tr_b - roi_tr;                               % [n_roi x freq x t_lpc]
        plv = abs(mean(exp(1i*pd), 'all'));
        isPC_per_trial(tr) = plv;
    end

    % Pull behavioral rows for this subject
    rows_b = B(B.sub_num == sub_num, :);
    % Map orig_epoch -> table row (the SingleTrial CSV is per-epoch already)
    % For each EEG epoch, find matching row by orig_epoch
    for ep = 1:n_trials
        % EEG.epoch(ep) loop variable ep = post-rejection sequential index,
        % which matches the CSV's epoch_idx column.
        idx = find(rows_b.epoch_idx == ep, 1);
        if isempty(idx), continue; end

        out_rows(end+1).sub = sid; %#ok<SAGROW>
        out_rows(end).sub_num = sub_num;
        out_rows(end).orig_epoch = ep;
        out_rows(end).condition  = char(rows_b.condition(idx));
        out_rows(end).item_name  = char(rows_b.item_name(idx));
        out_rows(end).familiarity= rows_b.familiarity(idx);
        out_rows(end).match_status = char(rows_b.match_status(idx));
        out_rows(end).isPC_LPC   = isPC_per_trial(ep);
    end
    fprintf('-> %d trials with merged behav.', sum([out_rows.sub_num] == sub_num));
end

T = struct2table(out_rows);
writetable(T, out_csv);
fprintf('\n\nWrote per-trial table: %s   (%d rows)\n', out_csv, height(T));

%% --- LME models ---
T.Condition = categorical(T.condition);
% Reorder: Literal as reference (matches ERP paper)
% Note: condition strings from CSV were 'Literal','Metaphor','Anomalous'
T.Condition = reordercats(T.Condition, intersect({'Literal','Metaphor','Anomalous'}, categories(T.Condition), 'stable'));
T.Sub      = categorical(T.sub);
T.Item     = categorical(T.item_name);
T.Familiarity_c = T.familiarity - mean(T.familiarity, 'omitnan');

T_matched   = T(strcmp(T.match_status,'matched') & ~isnan(T.isPC_LPC) & ~isnan(T.familiarity), :);
T_full      = T(~isnan(T.isPC_LPC) & ~isnan(T.familiarity), :);

fprintf('\n--- LME on MATCHED subset (n trials = %d) ---\n', height(T_matched));
lme_m = fitlme(T_matched, 'isPC_LPC ~ Condition + Familiarity_c + (1|Sub) + (1|Item)');
disp(lme_m)
[~, ~, FE_m] = fixedEffects(lme_m);
disp('Fixed effects:'); disp(FE_m);
anova_m = anova(lme_m);
disp('ANOVA:'); disp(anova_m);

fprintf('\n--- LME on FULL set (n trials = %d) ---\n', height(T_full));
lme_f = fitlme(T_full, 'isPC_LPC ~ Condition + Familiarity_c + (1|Sub) + (1|Item)');
disp(lme_f)
[~, ~, FE_f] = fixedEffects(lme_f);
disp('Fixed effects:'); disp(FE_f);
anova_f = anova(lme_f);
disp('ANOVA:'); disp(anova_f);

fprintf('\n--- LME with Condition*Familiarity interaction (full set) ---\n');
lme_int = fitlme(T_full, 'isPC_LPC ~ Condition*Familiarity_c + (1|Sub) + (1|Item)');
disp(lme_int)
anova_int = anova(lme_int);
disp('ANOVA:'); disp(anova_int);

% Save model summaries
diary(fullfile(out_dir,'lme_results.txt'));
fprintf('\n\n=================== LME RESULTS SUMMARY ===================\n');
fprintf('\n=== Matched subset ===\n'); disp(lme_m); disp(anova_m);
fprintf('\n=== Full set ===\n'); disp(lme_f); disp(anova_f);
fprintf('\n=== Full set + interaction ===\n'); disp(lme_int); disp(anova_int);
diary off;

fprintf('\n=== Done. CSV: %s   LME log: %s ===\n', out_csv, fullfile(out_dir,'lme_results.txt'));


% =========================================================================
function adj = local_chan_adjacency(chanlocs, max_dist)
n = numel(chanlocs); xyz = nan(n,3);
for k = 1:n
    if ~isempty(chanlocs(k).X) && ~isempty(chanlocs(k).Y) && ~isempty(chanlocs(k).Z)
        xyz(k,:) = [chanlocs(k).X, chanlocs(k).Y, chanlocs(k).Z];
    end
end
norms = sqrt(sum(xyz.^2,2)); mn = median(norms(~isnan(norms)));
if mn < 5, xyz = xyz * (90/mn); end
D = nan(n);
for i = 1:n
    for j = 1:n
        D(i,j) = norm(xyz(i,:)-xyz(j,:));
    end
end
adj = (D > 0) & (D <= max_dist);
end
