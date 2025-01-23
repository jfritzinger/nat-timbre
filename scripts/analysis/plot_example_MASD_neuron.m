%% plot_example_MASD_neuron


%% Load in data 

[base, datapath, ~, ppi] = getPaths();
savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/AN_neurogram';
sheetpath = 'data/2025-manuscript/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);

putative = 'R29_TT3_P5_N02';
filename = sprintf('%s.mat', putative);
load(fullfile(datapath,'neural_data', filename)), 'data';
index = find(cellfun(@(s) strcmp(putative, s), sessions.Putative_Units));
CF = sessions.CF(index);
MTF_shape = sessions.MTF{index};


%% Analyze natural timbre

% Plot NT
params_NT = data(14, 2); % 13 is oboe, 14 is bassoon
data_NT = cell(2, 1);
if ~isempty(params_NT{1})
	data_NT = analyzeNT(params_NT{1});
end


%% Read in MASD for each stimulus at the CF

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
[~, index] = ismember(note_names, tuning.Note);
F0s = tuning.Frequency(index);

% Sort files and frequencies by pitch
[F0s, order] = sort(F0s);
files = files(order);

% Initialize variables for later use (if needed)
nfiles = numel(files);
wav_npts = zeros(1, nfiles);
wav_data = cell(1, nfiles);
MASD_CF = NaN(1, nfiles);
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
	MASD = trapz(t, spec_diff, 2);

	% Find 
	[~, CF_ind] = min(abs(CF-AN.CFs));
	MASD_CF(ii) = MASD(CF_ind);

end


%% Plot both 

figure

% Plot natural timbre
nexttile
hold on
errorbar(F0s,data_NT.rate, data_NT.rate_std/sqrt(params_NT{1}.nrep), 'LineWidth',2)
xlabel('Pitch (Hz)')
plot(F0s, MASD_CF);
r = corrcoef(MASD_CF, data_NT.rate);
R = r(1,2);
title([extractBefore(params_NT{1}.filename{1}, '.') ', R=' num2str(R)])
set(gca, 'fontsize', 16)


