%% plot_models_pdf.m
%
% Plots each synthetic timbre response from each neuron with the model
% responses for SFIE, energy, and population SFIE overlayed. 
%
%
% Author: J. Fritzinger
% Created: 2022-09-13; Last revision: 2024-09-13
%
% -------------------------------------------------------------------------
clear

import mlreportgen.dom.*
import mlreportgen.report.*

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
modelpath = '/Volumes/Nat-Timbre/data/manuscript';
sheetpath = 'scripts/data-cleaning';
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(base, sheetpath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Initialize report
report_path = 'figures/pdfs/';
filename = 'Model_NT_BE';
images = {}; %hold all plots as images, need to delete when finished
datetime.setDefaultFormats('default','yyyy-MM-dd_hhmmss')
report_name = sprintf('%s%s_%s.pdf', fullfile(base, report_path), datetime, filename);
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

%% Plot each dataset 

% Find sessions with MTF 
MTF_target = 'BE';
isMTF = strcmp(sessions.MTF, MTF_target);
%isMTF = contains(sessions.MTF, MTF_target);

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);

bin200_MTF = bin200 & isMTF;
has_data = bin200_MTF(:,1) | bin200_MTF(:,2);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);

% Plot each neuron
dips_p = NaN(num_sessions, 4);
peaks_p = NaN(num_sessions, 4);
for isesh = 1:num_sessions
	ineuron = index(order(isesh)); %indices(isesh)
	if any(has_data(ineuron))

		% Load in data 
		putative = sessions.Putative_Units{ineuron};
		CF = sessions.CF(ineuron);
		MTF_shape = sessions.MTF{ineuron};
		load(fullfile(datapath, [putative '.mat']))

		% Load in model data 
		load(fullfile(modelpath,'SFIE_model', [putative '_SFIE.mat']), 'SFIE')
		load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')

		% Paragraph intro
		label = '';
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

		label = sprintf("%s, CF = %0.0fHz, %s\n", putative, CF, MTF_shape);
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

		% Output 
		 fprintf('Creating plots... %s, CF = %0.0fHz, %s\n', putative, CF, MTF_shape);

		% Set up figure
		fig = figure('Position',[1,571,8*ppi,5*ppi]);
		tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','tight')
		x_label = [1000 2000 3000 4000 6000 8000]./1000;
		fontsize = 7;

		% Plot RM
		params_RM = data{2, 2};
		if ~isempty(params_RM)
			data_RM = analyzeRM(params_RM);
			nexttile(1)
			hold on
			spont_color = [0.4 0.4 0.4];
			CF_color = [0.3 0.3 0.3];
			plot(data_RM.freqs./1000,data_RM.rates(:,5),'color', '#20116B','LineWidth',2) % 73 dB
			plot(data_RM.freqs./1000,data_RM.rates(:,4),'color', '#5E50A9','LineWidth',2) % 53 dB
			plot(data_RM.freqs./1000,data_RM.rates(:,3),'color', '#A49BD0','LineWidth',2) % 33 dB
			plot(data_RM.freqs([1 end])./1000,[1 1]*data_RM.spont,'-','LineWidth',2, 'Color',spont_color)
			xline(CF/1000, '--', 'Color', CF_color,'LineWidth',2);
			box on
			grid on
			hold off
			ylim([0 max(data_RM.rates, [], 'all')+10])
			set(gca,'XTick',[])
			xlim([250 14000]./1000)
			xticks(x_label)
			set(gca, 'XScale', 'log');
			set(gcf, 'color', 'w')
			set(gca,'fontsize',fontsize)
			ylabel('Avg. Rate (sp/s)')
			title('Response Map')
			xlabel('Frequency (kHz)')
		end

		% Plot MTF
		params_MTF = data{3, 2};
		if isempty(params_MTF)
			params_MTF = data{3, 1};
		end
		spont_color = [0.4 0.4 0.4];
		data_MTF = analyzeMTF(params_MTF);
		nexttile(2)
		hold on
		line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color, 'LineWidth',2);
		errorbar(data_MTF.fms,data_MTF.rate, data_MTF.rate_std/sqrt(params_MTF.nrep),'.', 'LineWidth',2, 'Color','k');
		line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', 2);
		hold off
		set(gca, 'XScale', 'log');
		xlim([data_MTF.fms(1) data_MTF.fms(end)])
		xticks([1 2 5 10 20 50 100 200 500])
		xlabel('Modulation Frequency')
		set(gca,'fontsize',fontsize)
		ylimit = ylim;
		ylim([0 ylimit(2)])
		grid on
		box on
		title(['MTF: ' MTF_target ', 200Hz=' data_MTF.at_200], 'FontSize', fontsize)
		ylabel('Avg. Rate (sp/s)')

		% Plot natural timbre (raw)
		instruments = {'Bassoon', 'Oboe'};
		data_colors = {'#82BB95', '#3F985C', '#03882F', '#034E1C'};
		for ispl = 1:2
			nexttile(2+ispl)
			if bin200_MTF(ineuron, ispl)==1
				param_ST = data(12+ispl, 2);
				data_NT = analyzeNT(param_ST{1});

				if ~isempty(data_NT.rate)

				% Z-score 
				rate = zscore(data_NT.rate);
				hold on
				plot(data_NT.pitch, rate,'-*', 'LineWidth',2);
				ylim([-4 4])

				% Normalize and plot models
				plot(data_NT.pitch, zscore(energy{ispl}.rate), 'LineWidth',1.5, 'Color','#D55E00')
				plot(data_NT.pitch, zscore(SFIE{ispl}.rate), 'LineWidth',1.5, 'Color','#009E73')				

				% Annotate SFIE model R^2
				lefts = linspace(0.04, 0.6, 2);
				message = sprintf('R^2 SFIE = %.02f', SFIE{ispl}.R2);
				annotation('textbox',[lefts(ispl) 0.38 0.2 0.0869],...
					'Color','k',...
					'String',message,...
					'FontSize',8,...
					'EdgeColor','none');

				% Annotate energy model R^2
				message = sprintf('R^2 energy = %.02f', energy{ispl}.R2);
				annotation('textbox',[lefts(ispl) 0.35 0.2 0.0869],...
					'Color','k',...
					'String',message, ...
					'FontSize',8,...
					'EdgeColor','none');

				plottitle = [instruments{ispl} ' dB SPL'];
				else
					plottitle = [instruments{ispl} ' dB SPL'];
				end
			else
				plottitle = [instruments{ispl} ' dB SPL'];
			end

			xlabel('F0 (Hz)')
			if ispl == 1
				ylabel('Z-score')
			end
			box on
			title(plottitle)
			grid on
			set(gca, 'fontsize', fontsize)
		end

		label = '';
		p = Paragraph(label);
		p.FontSize = "14pt";
		p.WhiteSpace = "preserve";
		append(rpt,p);

		
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


%% FUNCTIONS 

function [img, images] = addtoSTPDF(images, fig, title)
import mlreportgen.dom.*

% Set figure size, recommended
values = [8, 4.65];
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
