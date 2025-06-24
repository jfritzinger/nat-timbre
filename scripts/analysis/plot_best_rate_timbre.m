[base, ~, ~, ppi] = getPathsNT();

% Load in rate models
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat');
load(filepath, "neuron_rate_timbre")

% Load in all data
load(fullfile(fullfile(base, 'model_comparisons', 'Data_NT_3.mat')), 'nat_data')

%% Set up figure

figure('Position',[50, 50, 6.7*ppi, 3*ppi])
tiledlayout(5, 5)
linewidth = 1;
fontsize = 8;
labelsize = 12;
legsize = 6;
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsTimbre = '#1b9e77';
scattersize = 12;

%% A, B. Plot examples of accurate and inaccurate


pitch = getBassoonF0s();
ylimits = [100, 40; 65, 45];

accuracy_rate = [neuron_rate_timbre.accuracy]*100;
[ac, originalpos] = sort(accuracy_rate, 'descend');
h_ind = [1, 7];
ind_high=originalpos(1:25);


putative = {neuron_rate_timbre(ind_high).putative};
ind_b = 25:40;
ind_o = [1 3:17];
for iexample = 1:25
	nexttile

	% Get rates
	pitches = pitch(ind_b);
	index = find(strcmp(putative{iexample}, {nat_data.putative}));
	oboe_rate = nat_data(index).oboe_rate(ind_o);
	oboe_std = nat_data(index).oboe_rate_std(ind_o);
	bass_rate = nat_data(index).bass_rate(ind_b);
	bass_std = nat_data(index).bass_rate_std(ind_b);
	CF = nat_data(index).CF;
	MTF = nat_data(index).MTF;

	if mean(oboe_rate)>mean(bass_rate)
		rates(iexample) = 1;
	else
		rates(iexample) = 2;
	end

	errorbar(pitches, oboe_rate, oboe_std/sqrt(20), 'color',colorsData{1}, ...
		'LineWidth', linewidth, 'CapSize',3);
	hold on
	errorbar(pitches, bass_rate, bass_std/sqrt(20), 'color',colorsData{2}, ...
		'LineWidth', linewidth, 'CapSize',3);
	box off
	grid on
	xticks([50 100 250 500 1000])
	%ylim([0 ylimits(igood, iexample)])
	yticks([0 25 50 75 100])

	ylabel('Avg Rate (sp/s)')

	xlabel('F0 (kHz)')
	xticklabels([50 100 250 500 1000]/1000)
	set(gca, 'fontsize', fontsize)

	if iexample == 1
		hleg = legend('Oboe', 'Bassoon', 'fontsize', legsize, ...
			'numcolumns', 2, 'box', 'off');
		hleg.ItemTokenSize = [8, 8];
	end

	% Annotate 
	CF_msg = ['CF=' num2str(round(CF)) ' Hz'];
	MTF_msg = MTF;
	text(0.05, 0.95, CF_msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize)
	text(0.05, 0.85, MTF_msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize)
	all_CFs(iexample) = CF;
	all_accuracy(iexample) = accuracy_rate(ind_high(iexample));
end


%% Plot CF of best and wether bassoon or oboe is higher 

figure
tiledlayout(2, 1)

nexttile
scatter(all_CFs(rates==1), all_accuracy(rates==1), 'filled')
hold on
scatter(all_CFs(rates==2), all_accuracy(rates==2), 'filled')
legend('Oboe Greater Rate', 'Bassoon Greater Rate')
set(gca, 'xscale', 'log')
xticks([50 100 200 500 1000 2000 5000 10000 12000])
xlim([500 12000])
xlabel('CF (Hz)')
ylabel('Accuracy')

load(fullfile(base, 'stimuli_avg_env_diff.mat'), 'avg_env_diff', "freq")
nexttile
plot(freq, avg_env_diff, 'k')
hold on
yline(0, 'k')
set(gca, 'xscale', 'log')
xticks([50 100 200 500 1000 2000 5000 10000 12000])
xlim([500 12000])
xlabel('Frequency (Hz)')
ylabel('Envelope Difference (dB SPL)')
scatter(all_CFs(rates==1), zeros(length(all_CFs(rates==1)), 1), 'filled', ...
	'MarkerFaceColor','#0072BD', 'MarkerEdgeColor','k')
scatter(all_CFs(rates==2), zeros(length(all_CFs(rates==2)), 1), 'filled', ...
	'MarkerFaceColor','#D95319', 'MarkerEdgeColor','k')


%% Plot of accuracy and CF region 

x1 = 899.658;
x2 = 7075;

accuracy_rate = [neuron_rate_timbre.accuracy]*100;
CFs = [neuron_rate_timbre.CF];

oboe_ind = CFs > x1 & CFs < x2;
accuracy_x1 = accuracy_rate(oboe_ind);
accuracy_x2 = accuracy_rate(~oboe_ind);

figure
swarmchart(ones(195, 1), accuracy_x1)
hold on
swarmchart(ones(51, 1)*2, accuracy_x2)
xticks([1, 2])
xticklabels({'Oboe greater', 'Bassoon greater'})
xlabel('Freq regions')
ylabel('Accuracy')


%%

function pitch = getBassoonF0s()
target = 'Bassoon';
if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning
listing = dir(fullfile(fpath, 'waveforms', '*.wav'));
target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
wav_nums =  find(cell2mat(target_WAV));
d = dir(fullfile(fpath,'waveforms', '*.wav'));
all_files = sort({d.name});
nfiles = length(wav_nums);

for i = 1:nfiles
	files{1,i} = all_files{wav_nums(i)};
end

% Sort by frequency of pitch
index = zeros(1, nfiles);
note_names = extractBetween(files, 'ff.','.');
for ii = 1:nfiles % Find index of each note in tuning spreadsheet
	index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note
[~, order] = sort(pitch_order); % Sort freqs
pitch = pitch_order(order);

end