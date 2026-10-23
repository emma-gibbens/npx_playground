%% Summarize odor, reward, and airpuff encoding by region
% Reports: for each region (dorsal / ventral), what % of ALL CA1 neurons
% are odor-selective, reward-responsive, and airpuff-responsive.
%
% Denominator = all neurons in neuron_summary.csv (post-decorr-mask,
%   firing rate >= 0.3 Hz), including those flagged skip_low_baseline_fr.
%   Region (dorsal/ventral) is read directly from the CSV.
%
% Encoding numerators = neurons that additionally pass skip_low_baseline_fr.
%
% Output: per-session table + aggregate totals across all sessions.
%
% EG June 2026

%% -------------------------------------------------------------------------
%% Parameters
%% -------------------------------------------------------------------------

% folder_list = {
%     'D:\vCA1\M656\preprocessed\M656_2026_03_25'
%     'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
%     'D:\vCA1\M649\preprocessed\M649_2026_04_01'
%     'D:\vCA1\M593\preprocessed\M593_2025_12_26'
%     'D:\vCA1\M593\preprocessed\M593_2025_12_28'
%     'D:\vCA1\M653\preprocessed\M653_2026_02_04'
% };

folder_list = {
    'D:\vCA1\M560\preprocessed\M560_2025_12_19'
    'D:\vCA1\M593\preprocessed\M593_2025_12_26'
    'D:\vCA1\M593\preprocessed\M593_2025_12_28'
    'D:\vCA1\M653\preprocessed\M653_2026_02_04'
    'D:\vCA1\M656\preprocessed\M656_2026_03_25'
    'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
    'D:\vCA1\M649\preprocessed\M649_2026_04_01'
    'D:\vCA1\M656\preprocessed\M656_2026_03_26_conservativeSpikeSorting'
    };

%% -------------------------------------------------------------------------
%% Accumulators across all sessions
%% -------------------------------------------------------------------------

regions = {'Dorsal', 'Ventral'};
agg = struct();
for r = 1:length(regions)
    rn = regions{r};
    agg.(rn).total   = 0;
    agg.(rn).odor    = 0;
    agg.(rn).reward  = 0;
    agg.(rn).airpuff = 0;
end

fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('Odor / Reward / Airpuff encoding — %% of all CA1 neurons by region\n');
fprintf('%s\n\n', repmat('=', 1, 78));

%% -------------------------------------------------------------------------
%% Main loop
%% -------------------------------------------------------------------------

for iFolder = 1:length(folder_list)

    csv_path = fullfile(folder_list{iFolder}, 'neuron_summary.csv');
    if ~isfile(csv_path)
        warning('CSV not found, skipping: %s', folder_list{iFolder});
        continue;
    end
    data = readtable(csv_path, 'VariableNamingRule', 'preserve');

    sess = struct();
    for r = 1:length(regions)
        rn = regions{r};
        sess.(rn).total   = 0;
        sess.(rn).odor    = 0;
        sess.(rn).reward  = 0;
        sess.(rn).airpuff = 0;
    end

    for row = 1:height(data)

        % Region from CSV — capitalise first letter to match struct field names
        loc_raw = strtrim(data.region{row});   % e.g. 'dorsal' or 'ventral'
        loc = [upper(loc_raw(1)), lower(loc_raw(2:end))];   % → 'Dorsal' / 'Ventral'
        if ~ismember(loc, regions)
            continue;
        end

        % Every neuron in the CSV counts toward the denominator
        sess.(loc).total = sess.(loc).total + 1;

        % Encoding numerators: skip neurons flagged skip_low_baseline_fr
        skip_fr = data.skip_low_baseline_fr{row};
        if strcmpi(skip_fr, 'Yes')
            continue;
        end

        % --- read encoding fields ---
        odor_sel      = data.odor_selectivity{row};
        odor_resp_raw = data.odor_responsivity{row};
        reward_resp   = data.reward_responsive{row};
        airpuff_resp  = data.airpuff_responsive{row};

        % Odor selectivity: direct assignment OR single-odor-responsive promotion
        odor_sel_is_empty = isempty(odor_sel) || strcmpi(odor_sel, 'None') || ...
            (isscalar(odor_sel) && isnumeric(odor_sel) && isnan(odor_sel));
        odor_resp_is_empty = isempty(odor_resp_raw) || strcmpi(odor_resp_raw, 'None') || ...
            (isscalar(odor_resp_raw) && isnumeric(odor_resp_raw) && isnan(odor_resp_raw));

        is_odor_selective = ~odor_sel_is_empty;

        if ~is_odor_selective && ~odor_resp_is_empty
            responsive_odors = strtrim(strsplit(odor_resp_raw, ','));
            responsive_odors = responsive_odors(~cellfun('isempty', responsive_odors));
            if isscalar(responsive_odors)
                is_odor_selective = true;
            end
        end

        is_reward  = strcmpi(reward_resp,  'Yes');
        is_airpuff = strcmpi(airpuff_resp, 'Yes');

        if is_odor_selective,  sess.(loc).odor    = sess.(loc).odor    + 1; end
        if is_reward,          sess.(loc).reward  = sess.(loc).reward  + 1; end
        if is_airpuff,         sess.(loc).airpuff = sess.(loc).airpuff + 1; end

    end  % row loop

    % ------------------------------------------------------------------
    % Print session summary
    % ------------------------------------------------------------------
    [~, sess_name] = fileparts(folder_list{iFolder});
    fprintf('Session: %s\n', sess_name);
    fprintf('%-8s  %5s  %16s  %18s  %18s\n', ...
        'Region', 'N', 'Odor-selective', 'Reward-responsive', 'Airpuff-responsive');
    fprintf('%s\n', repmat('-', 1, 72));

    for r = 1:length(regions)
        rn = regions{r};
        N  = sess.(rn).total;
        if N == 0
            fprintf('%-8s  %5d  (no neurons)\n', rn, 0);
        else
            fprintf('%-8s  %5d  %4d (%5.1f%%)   %4d (%5.1f%%)     %4d (%5.1f%%)\n', ...
                rn, N, ...
                sess.(rn).odor,    100 * sess.(rn).odor    / N, ...
                sess.(rn).reward,  100 * sess.(rn).reward  / N, ...
                sess.(rn).airpuff, 100 * sess.(rn).airpuff / N);
        end

        % Aggregate
        agg.(rn).total   = agg.(rn).total   + sess.(rn).total;
        agg.(rn).odor    = agg.(rn).odor    + sess.(rn).odor;
        agg.(rn).reward  = agg.(rn).reward  + sess.(rn).reward;
        agg.(rn).airpuff = agg.(rn).airpuff + sess.(rn).airpuff;
    end
    fprintf('\n');

end  % folder loop

%% -------------------------------------------------------------------------
%% Aggregate summary
%% -------------------------------------------------------------------------

fprintf('%s\n', repmat('=', 1, 72));
fprintf('AGGREGATE — all sessions\n');
fprintf('%-8s  %5s  %16s  %18s  %18s\n', ...
    'Region', 'N', 'Odor-selective', 'Reward-responsive', 'Airpuff-responsive');
fprintf('%s\n', repmat('-', 1, 72));

for r = 1:length(regions)
    rn = regions{r};
    N  = agg.(rn).total;
    if N == 0
        fprintf('%-8s  %5d  (no neurons)\n', rn, 0);
    else
        fprintf('%-8s  %5d  %4d (%5.1f%%)   %4d (%5.1f%%)     %4d (%5.1f%%)\n', ...
            rn, N, ...
            agg.(rn).odor,    100 * agg.(rn).odor    / N, ...
            agg.(rn).reward,  100 * agg.(rn).reward  / N, ...
            agg.(rn).airpuff, 100 * agg.(rn).airpuff / N);
    end
end
fprintf('%s\n\n', repmat('=', 1, 72));
fprintf('Script completed.\n');