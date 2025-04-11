%% strf_model_predictions
clear


%% Load in spreadsheet

[base, datapath, ~, ppi] = getPaths();
spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

if ismac
	addpath('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/scripts/helper-functions', '-end')
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/strf_models';
else
	addpath('C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\scripts\helper-functions', '-end')
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

		if isesh == 30 || isesh == 42
			continue
		end

		% Load in data
		s_ind = index(isesh);
		putative_neuron = sessions.Putative_Units{s_ind};
		CF = sessions.CF(s_ind);
		load(fullfile(datapath, 'neural_data', [putative_neuron '.mat']), 'data');
		
		% STRF analysis
		params_STRF = data{4,2};
		data_STRF = analyzeSTRF(params_STRF);

		% NT analysis 
		param_NT = data{12+iinst,2};
		data_NT = analyzeNT(param_NT);

		% Calculate model response
		param_NT.Fs = 48000;
		param_NT.mnrep = param_NT.nrep;
		param_NT.dur = 0.3;
		[R2, avModel, stdModel, ratio, max_all, R] = modelNTSTRF(param_NT, data_STRF, data_NT);

		% Display progress
		fprintf('%s Done, %.2f\n', putative_neuron, isesh/length(index))

		%% Plot
		figure('Position',[3,632,623,273])
		tiledlayout(1, 2)

		% Plot STRF
		nexttile
		STRF_mat = data_STRF.H2ex_strf-data_STRF.H2in_strf;
		imagesc(data_STRF.t, data_STRF.f./1000, STRF_mat, data_STRF.clims_strf);
		set(gca,'Ydir','normal','XLim',data_STRF.tlims,'YLim',[0 10])
		hold on
		yline(CF/1000, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
		colormap(redblue)
		grid on
		xlabel('Time (s)');
		ylabel('Frequency (kHz)')

		% Plot real response
		nexttile
		hold on
		errorbar(data_NT.pitch,data_NT.rate,data_NT.rate_std/(sqrt(param_NT.nrep)), 'LineWidth',1.5);
		plot(data_NT.pitch,(avModel.*ratio), 'LineWidth',1.5);
		grid on
		xlabel('Tone Frequency (Hz)')
		ylabel('Avg Rate (sp/s)')
		ylim([0 max_all+5])
		title(sprintf('R^2 = %0.2f\n', R2))

		%% Add to struct
		temp(isesh).putative = putative_neuron;
		temp(isesh).R = R;
		temp(isesh).R2 = R2;
		temp(isesh).avModel = avModel;
		temp(isesh).stdModel = stdModel;
		temp(isesh).ratio = ratio;
		temp(isesh).max_all = max_all;
		temp(isesh).instrument = instruments{iinst};
		temp(isesh).pitch = data_NT.pitch;

		STRFmodel = temp(isesh);
		filename = sprintf('%s_STRF_%s', putative_neuron, instruments{iinst});
		save(fullfile(savepath, [filename '.mat']), "STRFmodel")
	end
end