%% NT_model_distributions 
clear 

%% Load in spreadsheet

[base, datapath, ~, ppi] = getPaths();
sheetpath = 'data-cleaning';
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(datapath, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/strf_models';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\strf_models';
end


%% Analyze

instruments = {'Oboe', 'Bassoon'};
for iinst = 2

	% Find sessions of interest
	if iinst == 1
		has_data = cellfun(@(s) contains(s, 'R'), sessions.Oboe) & cellfun(@(s) contains(s, 'R'), sessions.STRF);
	elseif  iinst == 2
		has_data = cellfun(@(s) contains(s, 'R'), sessions.Bassoon) & cellfun(@(s) contains(s, 'R'), sessions.STRF);
	end
	index = find(has_data);
	for isesh = 1:length(index) % error: 2: 30, 42

		% Load in data
		s_ind = index(isesh);
		putative_neuron = sessions.Putative_Units{s_ind};
		CF = sessions.CF(s_ind);
		load(fullfile(datapath, 'neural_data', [putative_neuron '.mat']), 'data');
		filename = sprintf('%s_STRF_%s', putative_neuron, instruments{iinst});
		try
			load(fullfile(savepath, [filename '.mat']), "STRFmodel")
			R2_all(isesh) = STRFmodel.R2;
		catch 
			R2_all(isesh) = NaN;
		end
	end
end

% Load in spreadsheet
datapath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
spreadsheet_name = 'model_r2_values_NT.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%% Plot 

figure('Position',[560,593,600,300])
tiledlayout(2, 2)
fontsize = 18;

isOboe = strcmp(sessions.Instrument, 'Oboe');
isBass = strcmp(sessions.Instrument, 'Bassoon');

nexttile
hold on
edges = linspace(0, 1, 20);
histogram(sessions.SFIE_R2(isOboe), edges, 'FaceAlpha',0.5)
histogram(sessions.SFIE_R2(isBass), edges, 'FaceAlpha',0.5)
title('SFIE')
set(gca, 'Fontsize', fontsize)
grid on
xlim([0 1])
ylabel('# Neurons')

nexttile
hold on
edges = linspace(0, 1, 20);
histogram(sessions.Energy_R2(isOboe), edges, 'FaceAlpha',0.5)
histogram(sessions.Energy_R2(isBass), edges, 'FaceAlpha',0.5)
title('Energy')
set(gca, 'Fontsize', fontsize)
grid on
xlim([0 1])

nexttile
hold on
edges = linspace(0, 1, 20);
histogram(sessions.Lat_Inh_R2(isOboe), edges, 'FaceAlpha',0.5)
histogram(sessions.Lat_Inh_R2(isBass), edges, 'FaceAlpha',0.5)
title('SFIE')
set(gca, 'Fontsize', fontsize)
legend('Oboe', 'Bassoon')
grid on
xlim([0 1])
xlabel('Variance Explained (R^2)')
ylabel('# Neurons')

nexttile
edges = linspace(0, 1, 21);
histogram(R2_all, edges)
xlabel('Variance Explained (R^2)')
title('Noise STRF')
set(gca, 'fontsize', 18)
xlim([0 1])
grid on
box off
