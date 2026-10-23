% compile_neuron_summary.m
%
% Concatenates neuron_summary.csv files from multiple sessions into a
% single mega_neuron_summary.csv. A 'session' column is added to every
% row (using the folder basename, e.g. 'M656_2026_03_25'), and global_id
% values are prefixed with the session name so they remain unique across
% the combined table.
%
% USAGE
%   1. Edit folder_list below (full paths to each session folder).
%   2. Set output_path to wherever you want the merged CSV saved.
%   3. Run the script.
%
% EXPECTED INPUT
%   Each folder must contain a neuron_summary.csv with at least a
%   'global_id' column. Sessions with different column sets are handled
%   gracefully — missing columns are filled with NaN (numeric) or
%   <missing> (string).
%
% OUTPUT
%   mega_neuron_summary.csv with columns:
%       session    <-- folder basename, e.g. 'M656_2026_03_25'
%       global_id  <-- prefixed: '<session>__<original_id>'
%       <all original columns>

% -------------------------------------------------------------------------
% USER SETTINGS
% -------------------------------------------------------------------------

folder_list = {
    'D:\vCA1\M656\preprocessed\M656_2026_03_25'
    'D:\vCA1\M649\preprocessed\M649_2026_03_30_conservativeSpikeSorting'
    'D:\vCA1\M649\preprocessed\M649_2026_04_01'
    'D:\vCA1\M593\preprocessed\M593_2025_12_26'
    'D:\vCA1\M593\preprocessed\M593_2025_12_28'
    'D:\vCA1\M653\preprocessed\M653_2026_02_04'
};

% Where to save the merged output
output_path = 'D:\vCA1\mega_neuron_summary.csv';

% -------------------------------------------------------------------------
% MAIN LOOP
% -------------------------------------------------------------------------

all_tables   = cell(numel(folder_list), 1);
skipped_sess = {};

for i = 1:numel(folder_list)
    folder    = folder_list{i};
    sess_name = get_session_name(folder);   % folder basename
    csv_path  = fullfile(folder, 'neuron_summary.csv');

    if ~isfile(csv_path)
        warning('compile_neuron_summary: file not found, skipping:\n  %s', csv_path);
        skipped_sess{end+1} = sess_name; %#ok<AGROW>
        continue
    end

    fprintf('[%d/%d] Loading %s\n', i, numel(folder_list), csv_path);

    % Read table, keeping text as strings
    T = readtable(csv_path, 'TextType', 'string');

    % Prefix global_id so IDs stay unique across sessions
    if ismember('global_id', T.Properties.VariableNames)
        T.global_id = string(sess_name) + "__" + T.global_id;
    else
        warning('compile_neuron_summary: no global_id column in %s', csv_path);
    end

    % Prepend session column
    T = addvars(T, repmat(string(sess_name), height(T), 1), ...
        'NewVariableNames', 'session', ...
        'Before', 1);

    all_tables{i} = T;
end

% Drop empties (skipped sessions)
all_tables = all_tables(~cellfun(@isempty, all_tables));

if isempty(all_tables)
    error('compile_neuron_summary: no valid session files found. Check folder_list paths.');
end

% -------------------------------------------------------------------------
% CONCATENATE (handles sessions with different column sets)
% -------------------------------------------------------------------------

fprintf('\nConcatenating %d session tables...\n', numel(all_tables));
mega = vertcat_flexible(all_tables);
fprintf('Merged table: %d neurons x %d columns\n', height(mega), width(mega));

% -------------------------------------------------------------------------
% WRITE OUTPUT
% -------------------------------------------------------------------------

writetable(mega, output_path);
fprintf('Saved: %s\n', output_path);

if ~isempty(skipped_sess)
    fprintf('\nWARNING: %d session(s) skipped (neuron_summary.csv not found):\n', numel(skipped_sess));
    fprintf('  %s\n', skipped_sess{:});
end


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function name = get_session_name(folder_path)
% Return the last component of a path (works on Windows \ and Unix /)
parts = strsplit(folder_path, {'\', '/'});
parts = parts(~cellfun(@isempty, parts));   % drop any empty tokens
name  = parts{end};
end


function out = vertcat_flexible(table_list)
% Vertically concatenate tables that may have different column sets.
% Missing columns are filled with NaN (numeric) or <missing> (string).

% Build master column list in order of first appearance
all_vars = {};
for k = 1:numel(table_list)
    all_vars = union(all_vars, table_list{k}.Properties.VariableNames, 'stable');
end

% Pad each table to the master column set, then reorder
for k = 1:numel(table_list)
    T            = table_list{k};
    missing_vars = setdiff(all_vars, T.Properties.VariableNames, 'stable');

    for v = 1:numel(missing_vars)
        col  = missing_vars{v};
        fill = get_fill(table_list, col, height(T));
        T.(col) = fill;
    end

    table_list{k} = T(:, all_vars);   % enforce column order
end

out = vertcat(table_list{:});
end


function fill = get_fill(table_list, col_name, n_rows)
% Scan tables to find col_name and return an appropriately-typed fill array.
for k = 1:numel(table_list)
    if ismember(col_name, table_list{k}.Properties.VariableNames)
        example = table_list{k}.(col_name);
        if isnumeric(example) || islogical(example)
            fill = NaN(n_rows, 1);
        else
            fill = repmat(missing, n_rows, 1);   % string <missing>
        end
        return
    end
end
fill = repmat(missing, n_rows, 1);   % fallback
end