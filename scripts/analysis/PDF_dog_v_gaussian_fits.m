%% PDF_dog_v_gaussian_fits
clear
import mlreportgen.dom.*
import mlreportgen.report.*

%% Load and initialize

if ismac
	addpath('/Users/jfritzinger/Projects/shared-models/DoG-model', '-begin')
else
	addpath('C:\Projects_JBF\shared-models\DoG-model', '-begin')
	addpath 'C:\Projects_JBF\nat-timbre\scripts\helper-functions'

end

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name),...
	'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Initialize report
filename = 'DoG_Gauss_Compare_NT'; 
images = {}; %hold all plots as images, need to delete when finished 
datetime.setDefaultFormats('default','yyyy-MM-dd_hhmmss') 
report_name = sprintf('%s/pdfs/%s_%s.pdf', savepath, datetime, filename); 
rpt = Document(report_name, 'pdf'); %initialize report document as pdf 
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.01in';
pm.PageMargins.Header = '0.01in';
pm.PageMargins.Bottom = '0.01in';
pm.PageMargins.Footer = '0.01in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';


%%

has_data = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);
linewidth = 1;
fontsize = 4.5;

% Plot each neuron
R2_dog_all = NaN(1, num_sessions);
R2_gauss_all = NaN(1, num_sessions);
for isesh = 1:num_sessions

	if ismember(isesh, [27, 117, 154, 181, 185, 187])
		continue
	end
 
	ineuron = index(order(isesh)); %indices(isesh)
	if any(has_data(ineuron))

		% Load in data
		putative = sessions.Putative_Units{ineuron};
		CF = sessions.CF(ineuron);
		MTF_shape = sessions.MTF{ineuron};
		load(fullfile(datapath, 'neural_data', [putative '.mat']))

		% Paragraph intro
		label = sprintf("%s, CF = %0.0fHz, %s\n", putative, CF, MTF_shape);
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

		% Output
		fprintf('Creating plots... %s, CF = %0.0fHz, %s\n', putative, CF, MTF_shape);

		% Calculate fits and plot

		% RM to get spont
		params_RM = data{2, 2};
		if isempty(params_RM)
			params_RM = data{2, 1};
		end
		data_RM = analyzeRM(params_RM);
		spont = data_RM.spont;

		% Natural timbre analysis
		params = data(14, 2); % 13 is oboe, 14 is bassoon
		data_NT = analyzeNT(params{1});

		% Generate stimulus
		params{1}.Fs = 100000;
		params{1}.mnrep = 1;
		params{1} = generate_NT(params{1});
		params{1}.num_stim = size(params{1}.stim, 1);

		%% Run models 

		Fs = 100000;
		observed_rate = data_NT.rate;
		r0 = spont;
		stim = params{1}.stim;
		nrep = 15;
		speed = 'slow'; % 'fast';
		dog_params = fit_dog_model(nrep, CF, Fs, stim, observed_rate, r0, speed);
		nrep = 10;
		gaussian_params = fit_gaussian_model(nrep, CF, Fs, stim, observed_rate, r0);

		%% Plot data
		fig = figure;
		tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','tight')
		nexttile
		hold on
		bar(data_NT.pitch, data_NT.rate)
		hold on
		errorbar(data_NT.pitch,data_NT.rate, data_NT.rate_std/sqrt(params{1}.nrep),...
			'LineStyle', 'none', 'Color', [0 0 0])
		xlabel('Pitch (Hz)')
		title(extractBefore(params{1}.filename{1}, '.'))
		title('Gaussian vs DoG Fits')
		ylabel('Avg. Rate (sp/s)')
		xlabel('Spectral Peak Freq. (kHz)')

		% Plot gaussian
		f = linspace(0, Fs/2, 100000);
		nstim = size(stim, 1);
		gaus_predicted = zeros(nstim, 1);
		for i = 1:nstim
			W = gaussian_model(f, gaussian_params);
			gaus_predicted(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
		end
		plot(data_NT.pitch, gaus_predicted, 'r', 'linewidth', linewidth)
		gaussian_adj_r_squared = calculate_adj_r_squared(observed_rate,...
			gaus_predicted, 3);

		% Plot DoG
		f = linspace(0, Fs/2, 100000);
		nstim = size(stim, 1);
		dog_predicted = zeros(nstim, 1);
		for i = 1:nstim
			W = dog_model(f, dog_params);
			dog_predicted(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
		end
		plot(data_NT.pitch, dog_predicted, 'g', 'linewidth', linewidth)
		dog_adj_r_squared = calculate_adj_r_squared(observed_rate,...
			dog_predicted, 5);

		% Annotations
		gaus_msg = sprintf('Gaus, \nR^2=%0.02f', gaussian_adj_r_squared);
		dog_msg = sprintf('DoG, \nR^2=%0.02f', dog_adj_r_squared);
		hleg = legend('','', gaus_msg, dog_msg, 'location', 'westoutside');
		hleg.ItemTokenSize = [4, 4];
		set(gca, 'FontSize',fontsize)

		% Get R^2 for all
		R2_dog_all(isesh) = dog_adj_r_squared;
		R2_gauss_all(isesh) = gaussian_adj_r_squared;

		%% Plot DoG Parameters
		nexttile
		hold on
		Fs = 100000;
		f = linspace(0, Fs/2, 100000);
		W = gaussian_model(f, gaussian_params);
		plot(f/1000,W, 'color', 'r')
		W2 = dog_model(f, dog_params);
		plot(f/1000,W2, 'color', 'g')

		% Plot labels
		xline(CF/1000, '--', 'linewidth', 1.5)
		title('DoG and Gaussian Kernels')
		set(gca, 'fontsize', fontsize)
		ylabel('Amplitude')
		xlabel('Frequency (kHz)')
		xlim([200 10000]/1000)
		set(gca, 'xscale', 'log')
		grid on
		xticks([0.1 0.2 0.5 1 2 5 10])

		% Add to PDF
		[plt1, images] = addtoSTPDF(images, fig, putative);
		append(rpt, plt1);

		%% Get f-test for all
		p_value = ftest(data_NT.rate, gaus_predicted, dog_predicted);

		% Struct to save out all data and fits 
		dog_gauss_analysis.putative = putative;
		dog_gauss_analysis.dog_predicted = dog_predicted;
		dog_gauss_analysis.gaus_predicted = gaus_predicted;
		dog_gauss_analysis.CF = CF;
		dog_gauss_analysis.rate = observed_rate;
		dog_gauss_analysis.R2_dog = dog_adj_r_squared;
		dog_gauss_analysis.R2_gauss = gaussian_adj_r_squared;
		dog_gauss_analysis.pitch = data_NT.pitch_num;
		dog_gauss_analysis.spont = spont;
		dog_gauss_analysis.rate_std = data_NT.rate_std;
		dog_gauss_analysis.p_value = p_value;
		dog_gauss_analysis.dog_params = dog_params;
		dog_gauss_analysis.gauss_params = gaussian_params;

		filename = [putative '.mat'];
		%savepath = '/Volumes/Synth-Timbre/data/manuscript/';
		savepath = 'C:\DataFiles_JBF\Nat-Timbre\data\manuscript';
        %savepath = '\\NSC-LCARNEY-H2\DataFiles_JBF\Nat-Timbre\data\manuscript';
		save(fullfile(savepath, 'dog_model', filename), 'dog_gauss_analysis')

	end
end

% Closes and opens PDF to view
close(rpt);
for i = 1:length(images)
	delete(images{1,i}.Path);
end
rptview(rpt)

save('R2_DOG.mat', "R2_gauss_all", "R2_dog_all")

%% FUNCTIONS

function [img, images] = addtoSTPDF(images, fig, title)
import mlreportgen.dom.*

% Set figure size, recommended
values = [4.15, 1.5];
fig.PaperSize = values;
fig.PaperPosition = [0 0 values];
fig.Units = 'inches';
fig.Position(3:4) = values;

% Add the plot to the document
name = sprintf('%s.svg', title);
print(fig, name, '-dsvg');
img = Image(name);
delete(fig) %delete plot figure window
images = [images {img}];

end
