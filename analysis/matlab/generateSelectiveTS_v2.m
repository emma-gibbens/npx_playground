%% Generate spiketrain structs for all onset selective cells
% Adapted for neuron_summary.csv
% Generates separate structs for:
%   all_selective  : odor-selective OR reward-responsive OR airpuff-responsive
%   odor_only      : odor-selective only (not reward/airpuff responsive)
%   reward_only    : reward-responsive only (not odor-selective or airpuff)
%   airpuff_only   : airpuff-responsive only (not odor-selective or reward)
% KL Sept 2024
% Adapted EG May 2026 for new data structure (neuron_summary.csv)

% List paths for all sessions
folder_list = {'D:\vCA1\M656\preprocessed\M656_2026_03_25'
               'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
               'D:\vCA1\M649\preprocessed\M649_2026_04_01'
               'D:\vCA1\M593\preprocessed\M593_2025_12_26'
               'D:\vCA1\M593\preprocessed\M593_2025_12_28'
               'D:\vCA1\M653\preprocessed\M653_2026_02_04'};

% Loop through all folders
for iFolder = 1:length(folder_list)
    cd(folder_list{iFolder});
    LoadExpKeys;

    % Load the CSV file
    csv_path = 'neuron_summary.csv';
    if ~isfile(csv_path)
        warning('File %s does not exist in folder %s. Skipping.', csv_path, folder_list{iFolder});
        continue;
    end

    data = readtable(csv_path, 'VariableNamingRule', 'preserve');

    %% Parse selectivity data
    % ---- all selective (union of odor / reward / airpuff) ----
    all_selective_cells      = [];
    all_selective_odors      = [];
    all_selective_odor_dirs  = [];   % ADDED: excitatory/inhibitory/None per cell
    all_selective_reward     = [];   % 'Yes'/'No' per cell
    all_selective_airpuff    = [];   % 'Yes'/'No' per cell

    % ---- odor selective only ----
    odor_selective_cells     = [];
    odor_selective_odors     = [];
    odor_selective_dirs      = [];   % ADDED: excitatory/inhibitory per cell

    % ---- reward responsive only ----
    reward_only_cells        = [];
    reward_only_direction    = [];   % e.g. 'increase' / 'decrease'

    % ---- airpuff responsive only ----
    airpuff_only_cells       = [];
    airpuff_only_direction   = [];   % e.g. 'increase' / 'decrease'

    % Process each row (cell)
    for row = 1:height(data)

        % --- read fields from new CSV ---
        global_id_raw   = data.global_id{row};           % e.g. imec0.shank0.56
        global_id_split = split(global_id_raw, '.');
        global_id       = strcat(global_id_split{1}, '_', global_id_split{3});
        % result: 'imec0_56'  (same format as before)

        odor_sel        = data.odor_selectivity{row};          % selective odor name or 'None'
        odor_sel_dir    = data.odor_selectivity_direction{row}; % ADDED: e.g. 'excitatory', 'inhibitory', 'None'
        airpuff_resp    = data.airpuff_responsive{row};         % 'Yes' or 'No'
        airpuff_dir     = data.airpuff_direction{row};          % e.g. 'increase', 'None'
        reward_resp     = data.reward_responsive{row};          % 'Yes' or 'No'
        reward_dir      = data.reward_direction{row};           % e.g. 'decrease', 'None'
        region          = data.region{row};                     % e.g. 'ventral', 'outside'
        skip_fr         = data.skip_low_baseline_fr{row};       % 'Yes' or 'No'

        % --- exclusion criteria ---
        % Skip cells outside the target region
        if strcmpi(region, 'outside') || strcmpi(region, 'outside_CA1')
            continue;
        end
        % Skip cells with low baseline firing rate
        if strcmpi(skip_fr, 'Yes')
            continue;
        end

        % --- classification flags ---
        is_odor_selective = ~strcmpi(odor_sel,     'None') && ~isempty(odor_sel);
        is_reward         = strcmpi(reward_resp,   'Yes');
        is_airpuff        = strcmpi(airpuff_resp,  'Yes');

        % Only process cells that belong to at least one category
        if ~is_odor_selective && ~is_reward && ~is_airpuff
            continue;
        end

        % --- add to the combined "all selective" pool ---
        all_selective_cells     = [all_selective_cells;    cellstr(global_id)];
        all_selective_odors     = [all_selective_odors;    cellstr(odor_sel)];
        all_selective_odor_dirs = [all_selective_odor_dirs; cellstr(odor_sel_dir)]; % ADDED
        all_selective_reward    = [all_selective_reward;   cellstr(reward_dir)];
        all_selective_airpuff   = [all_selective_airpuff;  cellstr(airpuff_dir)];

        % --- add to outcome-specific pools ---
        if is_odor_selective && ~is_reward && ~is_airpuff
            % Odor selective only
            odor_selective_cells = [odor_selective_cells; cellstr(global_id)];
            odor_selective_odors = [odor_selective_odors; cellstr(odor_sel)];
            odor_selective_dirs  = [odor_selective_dirs;  cellstr(odor_sel_dir)]; % ADDED

        elseif is_reward && ~is_odor_selective && ~is_airpuff
            % Reward responsive only
            reward_only_cells     = [reward_only_cells;     cellstr(global_id)];
            reward_only_direction = [reward_only_direction; cellstr(reward_dir)];

        elseif is_airpuff && ~is_odor_selective && ~is_reward
            % Airpuff responsive only
            airpuff_only_cells     = [airpuff_only_cells;     cellstr(global_id)];
            airpuff_only_direction = [airpuff_only_direction; cellstr(airpuff_dir)];
        end
        % Note: cells that are selective/responsive to more than one category
        % are included in all_selective but not in any outcome-specific pool.
    end

    %% Split each category by probe and save

    probes = {'imec0', 'imec1'};

    % UPDATED: added 'odor_dir' to both odor-containing categories
    categories = {
        'all',          all_selective_cells,   {all_selective_odors, all_selective_odor_dirs, all_selective_reward, all_selective_airpuff}, {'odor','odor_dir','reward_dir','airpuff_dir'};
        'odor_only',    odor_selective_cells,  {odor_selective_odors, odor_selective_dirs},                                                 {'odor','odor_dir'};
        'reward_only',  reward_only_cells,     {reward_only_direction},                                                                     {'reward_dir'};
        'airpuff_only', airpuff_only_cells,    {airpuff_only_direction},                                                                    {'airpuff_dir'};
    };
    % columns: label | cell_list | {metadata arrays} | {usr field names}

    for iProbe = 1:length(probes)
        probe = probes{iProbe};
        input_file = sprintf('clean_units_%s.mat', probe);

        % Check whether a spike file exists for this probe
        if ~isfile(input_file)
            continue;
        end

        % Load spike data only once per probe if any category has cells
        S_probe = [];

        for iCat = 1:size(categories, 1)
            cat_label    = categories{iCat, 1};
            cat_cells    = categories{iCat, 2};
            cat_meta     = categories{iCat, 3};   % cell array of metadata arrays
            cat_fields   = categories{iCat, 4};   % cell array of usr field names

            if ~isempty(cat_cells)

                % Find cells belonging to this probe
                probe_idx = find(contains(cat_cells, probe));
                if isempty(probe_idx)
                    continue;
                end
    
                probe_cells = cat_cells(probe_idx);
                probe_meta  = cellfun(@(m) m(probe_idx), cat_meta, 'UniformOutput', false);
    
                % Load spike struct on first use for this probe
                if isempty(S_probe)
                    S_probe = load_NPX_spike(input_file);
                end
    
                % Find indices in the spike struct
                spike_idx = [];
                for iCell = 1:length(probe_cells)
                    currCell = probe_cells{iCell};
                    currIdx  = find(contains(S_probe.label, currCell));
                    if ~isempty(currIdx)
                        spike_idx = [spike_idx; currIdx(1)];
                    end
                end
    
                if isempty(spike_idx)
                    continue;
                end
    
                % Build output struct
                S_out = SelectTS([], S_probe, spike_idx);
                for iField = 1:length(cat_fields)
                    S_out.usr.(cat_fields{iField}) = probe_meta{iField};
                end
    
                % Save
                out_filename = sprintf('selective_%s_%s.mat', probe, cat_label);
                save(out_filename, 'S_out');
                fprintf('Saved %d %s cells for %s -> %s\n', ...
                    length(spike_idx), cat_label, probe, out_filename);
            end
        end
    end

    fprintf('Session %d of %d completed\n', iFolder, length(folder_list));
    clearvars -except folder_list
end

fprintf('Script completed successfully!\n');