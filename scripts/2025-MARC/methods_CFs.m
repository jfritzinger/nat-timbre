%% Fig0_DistributionOfData
% J. Fritzinger, updated 12/15/23
%
% This script loads in the putative neurons spreadsheet and plots the MTF
% distribution, BMF distribution, WMF distribution, hybrid BMF/WMF
% distribution, and CF distribution for all neurons
clear
save_fig = 1;

%% Load in spreadsheet

[base, datapath, ~, ppi] = getPathsNT();
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(base, spreadsheet_name), 'PreserveVariableNames',true);
num_units = size(sessions, 1);

%% Number of units for each stimulus
 
% Natural timbre datasets
NT_datasets(1,:) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
NT_datasets(2,:) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
RVF = cellfun(@(s) contains(s, 'R'), sessions.RVF);

%% Set up figure

figure('Position',[194,1045,4*ppi,2.8*ppi])
fontsize = 18;
titlesize = 20;

%% Only get sessions with synthetic timbre 

nat(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
nat(:,2) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
nat(:,3) = cellfun(@(s) contains(s, 'R'), sessions.Other);
any_synth = any(nat(:,1:3), 2);
table = sessions(any_synth, :);

%% Get CFs for each putative neuron

% WBTIN Diotic
CFs = table.CF;
edges = [0 500 1000 2000 4000 8000 13000];
names = categorical({'<0.5', '0.5-1', '1-2', '2-4', '4-8', '8+'});
names = reordercats(names,{'<0.5', '0.5-1', '1-2', '2-4', '4-8', '8+'});
CF = CFs;
CF(CF==0) = [];
[N, edges1] = histcounts(CF, edges);

% Plot
nexttile
bar(names,N,'FaceColor', 'k', 'EdgeColor','k');
grid on
ylabel('# Neurons')
xlabel('CF (kHz)')
set(gca, 'FontSize', fontsize)
title('CF Distribution', 'fontsize', titlesize)
ylim([0 115])

%% Save figure 

if save_fig == 1
	filename = 'methods_CFs';
	save_figure_MARC(filename)
end