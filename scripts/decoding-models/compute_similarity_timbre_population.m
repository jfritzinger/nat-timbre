%% compute_similarity_timbre_population
clear 

%% Load in nat_data 

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

%% Get rate matrix of all data 

% Find all rows with bassoon and oboe
[sesh, num_data] = getTimbreSessions(nat_data);
T = getTimbrePopTable(nat_data, 'Rate', sesh, num_data, 'Model');
data_mat = table2array(T(:,1:end-1));

%% Compute similarity / correlation / variance explained / explainable variance? 


