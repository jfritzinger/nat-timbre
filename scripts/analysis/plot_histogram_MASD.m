%% plot_histogram_MASD
clear 

%% Load in MASD

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/AN_neurogram';
if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning
target = 'Bassoon';

% Get all .wav files containing the target instrument name
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};

% Extract note names and find corresponding frequencies
note_names = extractBetween(files, 'ff.', '.');
[~, ind] = ismember(note_names, tuning.Note);
F0s = tuning.Frequency(ind);

% Sort files and frequencies by pitch
[F0s, order] = sort(F0s);
files = files(order);

% Initialize variables for later use (if needed)
nfiles = numel(files);
wav_npts = zeros(1, nfiles);
wav_data = cell(1, nfiles);
for ii = 1:nfiles

	% Get single stimulus
	target = extractBefore(files{ii}, '.');
	target_F0 = extractBetween(files{ii}, 'ff.', '.wav');
	savefile = sprintf('%s_F0_%s_AN.mat', target, target_F0{1});
	load(fullfile(savepath, savefile),'AN', 'model_params', 'params')

	% Load in stimulus
	stim = params.stimulus(ii,:);

	% Plot MASD
	an_sout = squeeze(AN.an_sout);
	CFs = AN.CFs;
	t = linspace(0, params.dur, size(an_sout,2));
	spec_diff = diff(an_sout, 1);
	%spec_abs = abs(spec_diff);
	MASD(ii,:) = trapz(t, spec_diff, 2);
end


%% Load in data

[base, datapath, ~, ppi] = getPaths();
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);


% Find sessions for target synthetic timbre response
has_data = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
idata = find(has_data);
num_sessions = length(idata);

% Plot each neuron
for isesh = 1:num_sessions

	putative = sessions.Putative_Units{idata(isesh)};
	filename = sprintf('%s.mat', putative);
	load(fullfile(datapath,'neural_data', filename)), 'data';
	CF = sessions.CF(idata(isesh));
	MTF_shape = sessions.MTF{idata(isesh)};


	% Analyze natural timbre
	params_NT = data(14, 2); % 13 is oboe, 14 is bassoon
	if ~isempty(params_NT{1})
		data_NT = analyzeNT(params_NT{1});
	end

	% Read in MASD for each stimulus at the CF
	[~, CF_ind] = min(abs(CF-AN.CFs));
	if CF_ind == 100
		CF_ind = 99;
	end
	MASD_CF = MASD(:,CF_ind);

	if ~isempty(data_NT.rate)
		r = corrcoef(MASD_CF, data_NT.rate);
		R(isesh) = r(1,2);
	end
end


%% Histograms of results

figure
histogram(R)
title('Correlation between MASD and neuron rate')
ylabel('# Neurons')
xlabel('Correlation')
xlim([-1 1])
set(gca, 'fontsize', 16)






