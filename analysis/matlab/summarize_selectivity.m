%% =========================================================
%  Merge neuron metadata into a single summary CSV
%  Requires:
%    - odor_selectivity.csv
%    - stimulus_responsivity.csv
%    - selectivity_v12.csv
%    - clean_units_imec0.mat
%    - ExpKeys (loaded via LoadExpKeys) for probe depth/roll
% =========================================================

%% ---------- 0. User settings ----------------------------
% Set these paths to wherever your files live
csv_dir  = 'D:\vCA1\M656\preprocessed\M656_2026_03_26_conservativeSpikeSorting\';          % e.g. '/data/session1/'  (leave '' for pwd)
mat_file = fullfile(csv_dir, 'clean_units_imec0.mat');
out_file = fullfile(csv_dir, 'neuron_summary.csv');

% Probe geometry  –– loaded from ExpKeys below; edit defaults if needed
% (these are overwritten once LoadExpKeys runs)
record_depth_mm  = NaN;   % ExpKeys.probe1_Depth  (in mm -> *1000 = microns)
insertion_roll   = 0;     % ExpKeys.probe1_insertion_roll (degrees)

%% ---------- 1. Load ExpKeys & compute recording depths --
% Assumes LoadExpKeys is on the path and the current folder is correct
LoadExpKeys;
maxdepth = 1000 * ExpKeys.probe1_Depth;          % convert to microns
roll_rad = deg2rad(ExpKeys.probe1_insertion_roll);

S = load_NPX_spike(mat_file);

% Apply the same depth transform used in the original pipeline
recording_depths = (maxdepth - S.usr.recording_depths) .* cos(roll_rad);

% Unit IDs in the .mat come as 'imec0.<n>' – build a lookup table
% S.label is assumed to be a cell array of strings like 'imec0.46'
mat_ids     = S.label(:);          % Nx1 cell of strings
mat_depths  = recording_depths(:); % Nx1 double

%% ---------- 2. Read the three CSVs ----------------------
odor_T  = readtable(fullfile(csv_dir, 'odor_selectivity.csv'),    'TextType','string');
stim_T  = readtable(fullfile(csv_dir, 'stimulus_responsivity.csv'),'TextType','string');
sel_T   = readtable(fullfile(csv_dir, 'selectivity_v12', 'selectivity_v12.csv'),      'TextType','string');

% Standardise column names (case-insensitive safety)
odor_T.Properties.VariableNames = lower(odor_T.Properties.VariableNames);
stim_T.Properties.VariableNames = lower(stim_T.Properties.VariableNames);
sel_T.Properties.VariableNames  = lower(sel_T.Properties.VariableNames);

% Master ID list comes from odor_selectivity (use as the reference)
all_ids = odor_T.global_id;   % string array, e.g. "imec0.shank0.46"
N = numel(all_ids);

%% ---------- 3. Pre-allocate output columns --------------
odor_responsivity       = strings(N,1);
odor_selectivity        = strings(N,1);
odor_selectivity_dir    = strings(N,1);
airpuff_responsive      = strings(N,1);
airpuff_direction       = strings(N,1);
reward_responsive       = strings(N,1);
reward_direction        = strings(N,1);
skip_low_baseline_fr    = strings(N,1);
region                  = strings(N,1);

%% ---------- 4. Helper: convert 'None'/'' -> 'None' ------
clean = @(s) strtrim(s);
isNone = @(s) strcmpi(s,'none') | s == "";

%% ---------- 5. Loop over neurons ------------------------
for i = 1:N
    gid = all_ids(i);   % e.g. "imec0.shank0.46"

    % ---- 5a. Odor selectivity info ----------------------
    row_o = find(odor_T.global_id == gid, 1);
    if ~isempty(row_o)
        resp = clean(odor_T.responsivity(row_o));
        sel  = clean(odor_T.selectivity(row_o));
        sdir = clean(odor_T.selectivity_type(row_o));

        odor_responsivity(i)    = ternary(isNone(resp), "None", resp);
        odor_selectivity(i)     = ternary(isNone(sel),  "None", sel);
        odor_selectivity_dir(i) = ternary(isNone(sdir), "None", sdir);
    else
        odor_responsivity(i)    = "None";
        odor_selectivity(i)     = "None";
        odor_selectivity_dir(i) = "None";
    end

    % ---- 5b. Stimulus responsivity (airpuff / reward) ---
    row_s = find(stim_T.global_id == gid, 1);
    airpuff_responsive(i) = "No";
    airpuff_direction(i)  = "None";
    reward_responsive(i)  = "No";
    reward_direction(i)   = "None";

    if ~isempty(row_s)
        stim_str = clean(stim_T.responsive_stimuli(row_s));
        type_str = clean(stim_T.responsivity_types(row_s));

        if ~isNone(stim_str)
            % Split comma-separated lists
            stims = strtrim(strsplit(stim_str, ','));
            types = strtrim(strsplit(type_str, ','));

            for k = 1:numel(stims)
                if strcmpi(stims{k}, 'airpuff')
                    airpuff_responsive(i) = "Yes";
                    if k <= numel(types)
                        airpuff_direction(i) = string(types{k});
                    end
                elseif strcmpi(stims{k}, 'reward')
                    reward_responsive(i) = "Yes";
                    if k <= numel(types)
                        reward_direction(i) = string(types{k});
                    end
                end
            end
        end
    end

    % ---- 5c. skip_low_baseline_fr flag ----------
    % row_v = find(sel_T.neuron_id == gid, 1);
    % skip_low_baseline_fr(i) = "No";
    % if ~isempty(row_v)
    %     try
    %         if ~isnan(sel_T.onset_sel_reason(row_v)) || ~ismissing(sel_T.onset_sel_reason(row_v))
    %             reason = clean(sel_T.onset_sel_reason(row_v));
    %             if strcmpi(reason, 'skip_low_baseline_fr')
    %                 skip_low_baseline_fr(i) = "Yes";
    %             end
    %         end
    %     catch
    %          if ~isnan(sel_T.onset_sel_reason{row_v}) || ~ismissing(sel_T.onset_sel_reason{row_v})
    %              reason = clean(sel_T.onset_sel_reason(row_v));
    %              if strcmpi(reason, 'skip_low_baseline_fr')
    %                  skip_low_baseline_fr(i) = "Yes";
    %              end
    %          end
    %      end
    % end

    % ---- 5c. skip_low_baseline_fr flag ----------
    row_v = find(sel_T.neuron_id == gid, 1);
    skip_low_baseline_fr(i) = "No";

    if ~isempty(row_v)
        % 1. Check if the row contains a missing value BEFORE extracting the content
        if ~ismissing(sel_T.onset_sel_reason(row_v))

            % 2. Safe to extract and clean now
            reason = clean(sel_T.onset_sel_reason{row_v});

            if strcmpi(reason, 'skip_low_baseline_fr')
                skip_low_baseline_fr(i) = "Yes";
            end
        end
    end

    % ---- 5d. Depth / region from .mat -------------------
    % Convert "imec0.shank0.46" -> "imec0.46"
    % Pattern: imec0.shankX.ID  -> imec0.ID
    short_id = regexprep(gid, 'imec0\.shank\d+\.', 'imec0_'); % I apologize for how unreadable this is

    mat_idx = find(strcmpi(mat_ids, short_id), 1);
    if ~isempty(mat_idx)
        depth = mat_depths(mat_idx);
        if depth >= 3000
            region(i) = "ventral";
        else
            region(i) = "dorsal";
        end
    else
        region(i) = "unknown";
        warning('No depth found in .mat for neuron: %s (looked for %s)', gid, short_id);
    end
end

%% ---------- 6. Build & write output table ---------------
out_T = table( ...
    all_ids,             ...
    odor_responsivity,   ...
    odor_selectivity,    ...
    odor_selectivity_dir,...
    airpuff_responsive,  ...
    airpuff_direction,   ...
    reward_responsive,   ...
    reward_direction,    ...
    skip_low_baseline_fr,    ...
    region,              ...
    'VariableNames', {   ...
        'global_id',             ...
        'odor_responsivity',     ...
        'odor_selectivity',      ...
        'odor_selectivity_direction', ...
        'airpuff_responsive',    ...
        'airpuff_direction',     ...
        'reward_responsive',     ...
        'reward_direction',      ...
        'skip_low_baseline_fr', ...
        'region'                 ...
    });

writetable(out_T, out_file);
fprintf('Saved summary to: %s\n', out_file);

clearvars
%% ---------- 7. Inline helper function -------------------
function out = ternary(cond, a, b)
% Simple ternary: returns a if cond is true, else b
    if cond
        out = a;
    else
        out = b;
    end
end