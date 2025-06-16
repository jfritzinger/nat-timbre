%% 
clear
save_fig = 1;

%% Load in data

[base, ~, ~, ppi] = getPathsNT();

% Load in rate models
filepath = fullfile(base, 'model_comparisons','Neuron_Rate_Timbre_All.mat');
load(filepath, "neuron_rate_timbre")

% Load in all data
load(fullfile(fullfile(base, 'model_comparisons', 'Data_NT_3.mat')), 'nat_data')

%% Set up figure

figure('Position',[50, 50, 5.5*ppi, 5.3*ppi])
linewidth = 2;
fontsize = 18;
titlesize = 20;
legsize = 20;
scattersize = 40;
colorsData = {"#0072BD", "#D95319"};
colorsMTF = {'#648FFF', '#DC267F', '#785EF0', '#FFB000'};
colorsTimbre = '#1b9e77';


pitch = getBassoonF0s();
ylimits = [100, 40; 65, 45];
hloc = [1, 2; 3, 4];
for igood = 1:2

	accuracy_rate = [neuron_rate_timbre.accuracy]*100;
	if igood == 1
		[ac, originalpos] = sort(accuracy_rate, 'descend');
		ind_high=originalpos([2 3]);
	else
		[ac, originalpos] = sort(accuracy_rate, 'ascend');
		ind_high=originalpos([3 8]);
	end
	
	putative = {neuron_rate_timbre(ind_high).putative};
	ind_b = 25:40;
	ind_o = [1 3:17];
	for iexample = 1:2
		h(hloc(iexample,igood)) = subplot(2, 2, hloc(iexample,igood));

		% Get rates
		pitches = pitch(ind_b);
		index = find(strcmp(putative{iexample}, {nat_data.putative}));
		oboe_rate = nat_data(index).oboe_rate(ind_o);
		oboe_std = nat_data(index).oboe_rate_std(ind_o);
		bass_rate = nat_data(index).bass_rate(ind_b);
		bass_std = nat_data(index).bass_rate_std(ind_b);
		CFs(igood, iexample) = nat_data(index).CF;
		MTFs{igood, iexample} = nat_data(index).MTF;

		errorbar(pitches, oboe_rate, oboe_std/sqrt(20), 'color',colorsData{1}, ...
			'LineWidth', linewidth, 'CapSize',3);
		hold on
		errorbar(pitches, bass_rate, bass_std/sqrt(20), 'color',colorsData{2}, ...
			'LineWidth', linewidth, 'CapSize',3);
		box off 
		grid on
		xticks([50 100 250 500 1000])
		ylim([0 ylimits(igood, iexample)])
		yticks([0 25 50 75 100])

		if iexample == 1 && igood == 1
			xticklabels([])
			ylabel('Avg Rate (sp/s)                           ')
		elseif iexample == 2
			xlabel('F0 (kHz)')
			xticklabels([50 100 250 500 1000]/1000)
		else
			xticklabels([])
		end
		set(gca, 'fontsize', fontsize)

		if iexample == 1 && igood == 2
			hleg = legend('Oboe', 'Bassoon', 'fontsize', legsize, ...
				'numcolumns', 2, 'box', 'off', 'position', [0.397905527334233,0.000030445584881,0.493044822256568,0.063131313131313]);
			hleg.ItemTokenSize = [8, 8];
		end
	end
end

% Annotate 
left = [0.16 0.61]; 
bottom = linspace(0.2, 0.62, 2);
height = 0.27;
width = 0.35;

set(h(1), 'position', [left(1) bottom(2) width height])
set(h(2), 'position', [left(2) bottom(2) width height])
set(h(3), 'position', [left(1) bottom(1) width height])
set(h(4), 'position', [left(2) bottom(1) width height])

left = [0.15 0.6];
bottom = [0.53 0.95];
for ii = 1:2 % igood
	for iii = 1:2 % iexample
		msg = sprintf('CF = %0.0f Hz', CFs(ii,iii));
		annotation('textbox',[left(ii) bottom(iii) 0.5 0.058],...
			'String',msg,'FontSize',fontsize,...
			'EdgeColor','none');
		msg = sprintf('%s MTF', MTFs{ii,iii});
		annotation('textbox',[left(ii) bottom(iii)-0.05 0.5 0.058],...
			'String',msg,'FontSize',fontsize,...
			'EdgeColor','none');
		% annotation('textbox',[0.65 0.94 0.25 0.058],...
		% 	'String','Oboe','FontSize',fontsize,...
		% 	'EdgeColor','none');
		% annotation('textbox',[0.05 0.7 0.25 0.058],...
		% 	'String','Rate','FontSize',fontsize,...
		% 	'EdgeColor','none', 'Rotation',90);
		% annotation('textbox',[0.05 0.2 0.25 0.058],...
		% 	'String','Timing','FontSize',fontsize,...
		% 	'EdgeColor','none', 'Rotation',90);
	end
end

%% Save figure 
if save_fig == 1
	filename = 'timbre_rate_examples';
	save_figure_MARC(filename)
end

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