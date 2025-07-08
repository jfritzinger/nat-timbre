%% plot_example_responses 
clear
save_fig = 0;

%% Set up figure 

% Load in rate models
[base, ~, ~, ppi] = getPathsNT();
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat');
load(filepath, "neuron_rate_timbre")

% Load in all data
load(fullfile(fullfile(base, 'model_comparisons', 'Data_NT_3.mat')), 'nat_data')

%% Set up figure

figure('Position',[50, 50, 7*ppi, 5.2*ppi])
linewidth = 1;
fontsize = 8;
labelsize = 12;
legsize = 6;
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsTimbre = '#1b9e77';
scattersize = 12;

%% A, B. Plot examples of accurate and inaccurate

pitch = getF0s('Bassoon');
pitch_oboe = getF0s('Oboe');

ylimits = [100, 40; 65, 45];
accuracy_rate = [neuron_rate_timbre.accuracy]*100;
[ac, ind_high] = sort(accuracy_rate, 'descend');
[ac_lo, ind_low] = sort(accuracy_rate, 'ascend');

putative_all = {neuron_rate_timbre(ind_high).putative};
putative_lo = {neuron_rate_timbre(ind_low).putative};

% 3, 8, 6
putative = [putative_all([1 14, 13, 2, 8, 18, 3, 6, 9]) putative_lo([3, 8, 5])];
accuracy_all = [ac([1 14, 13, 2, 8, 18, 3, 6, 9]) ac_lo([3, 8, 5])];
ind_b = 25:40;
ind_o = [1 3:17];
for iexample = 1:12
	h(iexample) = subplot(4, 4, iexample);

	% Get rates
	pitches = pitch(ind_b);
	index = find(strcmp(putative{iexample}, {nat_data.putative}));
	% oboe_rate = nat_data(index).oboe_rate(ind_o);
	% oboe_std = nat_data(index).oboe_rate_std(ind_o);
	% bass_rate = nat_data(index).bass_rate(ind_b);
	% bass_std = nat_data(index).bass_rate_std(ind_b);
	oboe_rate = nat_data(index).oboe_rate;
	oboe_std = nat_data(index).oboe_rate_std;
	bass_rate = nat_data(index).bass_rate;
	bass_std = nat_data(index).bass_rate_std;
	CF = nat_data(index).CF;
	MTF = nat_data(index).MTF;

	if mean(oboe_rate)>mean(bass_rate)
		rates(iexample) = 1;
	else
		rates(iexample) = 2;
	end

	% errorbar(pitches, oboe_rate, oboe_std/sqrt(20), 'color',colorsData{1}, ...
	% 	'LineWidth', linewidth, 'CapSize',3);
	fill([pitch(ind_b(1)) pitch(ind_b(end)) pitch(ind_b(end)) ...
		pitch(ind_b(1))], [0 0 max([oboe_rate; bass_rate]*1.3) ...
		max([oboe_rate; bass_rate]*1.3)], 'k', 'FaceAlpha',0.1, ...
		'EdgeColor','none')
	hold on
	errorbar(pitch_oboe, oboe_rate, oboe_std/sqrt(20), 'color',colorsData{1}, ...
		'LineWidth', linewidth, 'CapSize',3);
	
	% errorbar(pitches, bass_rate, bass_std/sqrt(20), 'color',colorsData{2}, ...
	% 	'LineWidth', linewidth, 'CapSize',3);
	errorbar(pitch, bass_rate, bass_std/sqrt(20), 'color',colorsData{2}, ...
		'LineWidth', linewidth, 'CapSize',3);
	box off
	grid on
	xticks([50 100 250 500 1000])
	ylim([0 max([oboe_rate; bass_rate]*1.3)])
	yticks([0 25 50 75 100])

	if ismember(iexample, [1 4 7 10]) 
		ylabel('Avg Rate (sp/s)')
	end

	if ismember(iexample, [10, 11, 12])
		xlabel('F0 (kHz)')
		xticklabels([50 100 250 500 1000]/1000)
	else
		xticklabels([])
	end
	
	set(gca, 'fontsize', fontsize, 'XScale', 'log')

	if iexample == 1
		hleg = legend('', 'Oboe', 'Bassoon', 'fontsize', legsize, ...
			'numcolumns', 2, 'box', 'off', 'location', 'southoutside');
		hleg.ItemTokenSize = [8, 8];
	end

	% Annotate 
	CF_msg = ['CF=' num2str(round(CF)) ' Hz'];
	acc_msg = sprintf('Accuracy=%0.0f%%', accuracy_all(iexample));
	text(0.05, 0.95, CF_msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize)
	text(0.05, 0.85, MTF, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize)
	text(0.5, 0.95, acc_msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize)
end

%% 

CFs = [neuron_rate_timbre.CF];
MTFs = {neuron_rate_timbre.MTF};
MTF_types = unique(MTFs);


h(13) = subplot(4, 4, 13); 
scatter(CFs, accuracy_rate, scattersize, 'filled', 'MarkerEdgeColor','k', ...
	'MarkerFaceColor',colorsTimbre, 'MarkerFaceAlpha',0.5);
hold on
set(gca, 'xscale', 'log')

mdl = fitlm(CFs, accuracy_rate);
x = linspace(300, 14000, 50);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, ':k', 'LineWidth',2)
hleg = legend('', ...
	sprintf('p=%0.04f',mdl.Coefficients{2,4}), 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hleg.ItemTokenSize = [8, 8];
ylabel('Accuracy')
xlabel('CFs')
xticks([200 500 1000 2000 5000 10000])
xticklabels([0.2 0.5 1 2 5 10])
set(gca, 'fontsize', fontsize)
grid on
ylim([38 100])

h(14) = subplot(4, 4, 14); 
[weights_ordered, order_ind] = sort(accuracy_rate);
hold on
for iMTF = 1:4
	if iMTF == 4
		ind = strcmp(MTFs, MTF_types{iMTF}) | strcmp(MTFs, MTF_types{iMTF+1});
	else
		ind = strcmp(MTFs, MTF_types{iMTF});
	end
	[weights_ordered, order_ind] = sort(accuracy_rate(ind));
	num_units = length(weights_ordered);

	RGB = hex2rgb(colorsMTF{iMTF});
	swarmchart(ones(num_units, 1)*iMTF, weights_ordered, scattersize, RGB)

	mean_vals(iMTF) = mean(weights_ordered);
	std_vals(iMTF) = std(weights_ordered)/sqrt(length(weights_ordered));
end
%errorbar(1:4, mean_vals, std_vals, 'k')
xticks(1:4)
xticklabels({'BE', 'BS', 'F', 'H'})
xlabel('MTF Groups')
%ylabel('Accuracy')
set(gca, 'fontsize', fontsize)
grid on
ylim([38 100])
yticklabels([])

% Plot envelope 
h(15) = subplot(4, 4, 15); 
load(fullfile(base, 'stimuli_avg_env_diff.mat'), 'avg_env_diff', ...
	"freq", "env_diff")

plot(freq, env_diff, 'color', [0.3 0.3 0.3])
hold on
plot(freq, avg_env_diff, 'k', 'linewidth', 2)
yline(0)
set(gca, 'xscale', 'log', 'fontsize', fontsize)
ylabel({'Oboe - Bassoon'; 'Envelope (dB SPL)'})
xlabel('Frequency (Hz)')
xticks([100 200 500 1000 2000 5000 10000])
xticklabels([0.1 0.2 0.5 1 2 5 10])
grid on
box off
cross = find(diff(avg_env_diff>0));
xline(freq(cross(1)), 'b')
xline(freq(cross(2)), 'b')
xlim([500 10000])

% Plot data 

for iexample = 1:length(ind_high)
	index = find(strcmp(putative_all{iexample}, {nat_data.putative}));
	oboe_rate = nat_data(index).oboe_rate(ind_o);
	bass_rate = nat_data(index).bass_rate(ind_b);
	CF = nat_data(index).CF;
	MTF = nat_data(index).MTF;
	all_CFs(iexample) = CF;
	all_accuracy(iexample) = accuracy_rate(ind_high(iexample));

	if mean(oboe_rate)>mean(bass_rate)
		rates(iexample) = 1;
	else
		rates(iexample) = 2;
	end

end


h(16) = subplot(4, 4, 16); 
scatter(all_CFs(rates==1), all_accuracy(rates==1), scattersize, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.9)
hold on
scatter(all_CFs(rates==2), all_accuracy(rates==2), scattersize, 'filled',...
	'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.9)
set(gca, 'xscale', 'log', 'fontsize', fontsize)
ylabel('Accuracy')
xlabel('Frequency (Hz)')
xticks([100 200 500 1000 2000 5000 10000])
xticklabels([0.1 0.2 0.5 1 2 5 10])
xlim([500 10000])
grid on
ylim([60 100])
xline(freq(cross(1)), 'b')
xline(freq(cross(2)), 'b')
legend('Oboe Rate > Bassoon Rate', 'Bassoon Rate > Oboe Rate', 'box', 'off', ...
	'Location',"northwest")

%% Arrange 

left = [linspace(0.08, 0.47, 3) 0.72 0.87];
bottom = linspace(0.08, 0.78, 4);
height = 0.19;
width = 0.16;

set(h(1), 'position', [left(1) bottom(4) width height])
set(h(2), 'position', [left(2) bottom(4) width height])
set(h(3), 'position', [left(3) bottom(4) width height])
set(h(4), 'position', [left(1) bottom(3) width height])
set(h(5), 'position', [left(2) bottom(3) width height])
set(h(6), 'position', [left(3) bottom(3) width height])
set(h(7), 'position', [left(1) bottom(2) width height])
set(h(8), 'position', [left(2) bottom(2) width height])
set(h(9), 'position', [left(3) bottom(2) width height])
set(h(10), 'position', [left(1) bottom(1) width height])
set(h(11), 'position', [left(2) bottom(1) width height])
set(h(12), 'position', [left(3) bottom(1) width height])

bottom = linspace(0.08, 0.73, 3);
height = 0.22;
width = 0.125;

set(h(13), 'position', [left(4) bottom(3) width height])
set(h(14), 'position', [left(5) bottom(3) width height])
set(h(15), 'position', [left(4) bottom(2) 0.26 height])
set(h(16), 'position', [left(4) bottom(1) 0.26 height])

%% Set labels 

labelleft= [0 0.65];
labelbottom = linspace(0.25, 0.94, 4);
annotation('textbox',[labelleft(1) labelbottom(4) 0.071 0.058],...
	'String','A','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) labelbottom(3) 0.071 0.058],...
	'String','B','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) labelbottom(2) 0.071 0.058],...
	'String','C','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(1) labelbottom(1) 0.071 0.058],...
	'String','D','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) labelbottom(4) 0.071 0.058],...
	'String','E','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) 0.62 0.071 0.058],...
	'String','F','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');
annotation('textbox',[labelleft(2) 0.28 0.071 0.058],...
	'String','G','FontWeight','bold','FontSize',labelsize,...
	'EdgeColor','none');

%% Save figure 

if save_fig == 1
	filename = 'fig3_plot_example_responses';
	save_figure(filename)
end

%%

function pitch = getF0s(target)
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