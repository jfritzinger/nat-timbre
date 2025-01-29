%% PDF_dog_v_gaussian_fits
clear
import mlreportgen.dom.*
import mlreportgen.report.*

%% Load and initialize

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name),...
	'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Initialize report
filename = 'DoG_Gauss_Compare'; 
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
fontsize = 10;

% Plot each neuron
R2_dog_all = NaN(1, num_sessions);
R2_gauss_all = NaN(1, num_sessions);
for isesh = 1:num_sessions
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
		Fs = 100000;
		observed_rate = data_NT.rate;
		r0 = spont;

		% fmincon
		type = 2; % 1: distance, 2: MSE
		stim = params{1}.stim;
		timerVal = tic;

		% Fit gaussian
		init = [CF, 30, 100]; % Initial guess
		lb = [CF/4, 0, 0]; % Lower bounds
		ub = [CF*4, 5000, Inf]; % Upper bounds

		options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-12, ...
			'MaxFunEvals', 10^12, 'maxiterations', 1000, 'ConstraintTolerance', 1e-12, ...
			'StepTolerance', 1e-16, 'display', 'off');
		gaussian_params = fmincon(@(p) ...
			dog_objective_function(p, 'gaussian', Fs, stim, observed_rate, r0, type), ...
			init, [], [], [], [], lb, ub, [], options);


		% Fit DoG model
		dog_init = [20000, 10000, 100, 500,  CF, CF]; % Initial guess
		dog_lb = [100,   100,     10,  10,   CF/4, CF/4]; % Lower bounds
		dog_ub = [30000, 30000,   1000,1000, CF*4, CF*4]; % Upper bounds

		options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-12, ...
			'MaxFunEvals', 10^12, 'maxiterations', 1000, 'ConstraintTolerance', 1e-12, ...
			'StepTolerance', 1e-16, 'display', 'off');
		dog_params = fmincon(@(p) dog_objective_function(p, 'dog', Fs, stim, observed_rate, r0, type), ...
			dog_init, [], [], [], [], dog_lb, dog_ub, [], options);
		disp(['Model took ' num2str(toc(timerVal)) ' seconds'])

		% Plot data
		fig = figure;
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
			fc = gaussian_params(1);
			sigma = gaussian_params(2);
			g = gaussian_params(3);
			W = gaussian_model(f, fc, sigma, g);
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
		%legend('Data', 'CF', 'Spont', 'Gaussian', 'DoG', 'location', 'westoutside')
		set(gca, 'FontSize',fontsize)


		% Annotations
		gaus_msg = sprintf('Gaussian adjusted R^2=%0.02f', gaussian_adj_r_squared);
		text(0.05, 0.95, gaus_msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',fontsize)
		dog_msg = sprintf('DoG adjusted R^2=%0.02f', dog_adj_r_squared);
		text(0.05, 0.85, dog_msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',fontsize)

		% Get R^2 for all
		R2_dog_all(isesh) = dog_adj_r_squared;
		R2_gauss_all(isesh) = gaussian_adj_r_squared;

		% Add to PDF
		[plt1, images] = addtoSTPDF(images, fig, putative);
		append(rpt, plt1);

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
values = [5.5, 3];
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
