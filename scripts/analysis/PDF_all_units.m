%% PDF_all_units 

clear
import mlreportgen.dom.*
import mlreportgen.report.*

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Initialize report
filename = 'NT_Oboe';
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

target = 13; % 13 - oboe, 14 - bassoon

%% Plot each dataset 

% Find sessions for target synthetic timbre response
%bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
has_data = bin200(:,1); % | bin200(:,2);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);

% Plot each neuron
for isesh = 1:num_sessions
	% if isesh == 27 || isesh == 181 || isesh == 185
	% 	continue
	% end
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

		% Set up figure
		fig = figure('Position',[292,274,1264,420]);
		tiledlayout(4, 3, 'TileSpacing','compact', 'Padding','tight', ...
			'TileIndexing','columnmajor')
		x_label = [1000 2000 3000 4000 6000 8000]./1000;
		fontsize = 10;

		% Plot RM
		params_RM = data{2, 2};
		nexttile()
		if ~isempty(params_RM)
			data_RM = analyzeRM(params_RM);
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
			%legend('73dB SPL', '53dB SPL', '33dB SPL', 'Spont. Rate', 'fontsize',fontsize-2, 'location', 'best')
			title('Response Map')
			xlabel('Frequency (kHz)')
		end

		% Plot MTF
		params_MTF = data{3, 2};
		nexttile()
		if ~isempty(params_RM)
			data_MTF = analyzeMTF(params_MTF);
			hold on
			line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color, 'LineWidth',2);
			errorbar(data_MTF.fms,data_MTF.rate, data_MTF.rate_std/sqrt(params_MTF.nrep),'.', 'LineWidth',2, 'Color','k');
			line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', 'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', 2);
			hold off
			set(gca, 'XScale', 'log');
			xlim([data_MTF.fms(1) data_MTF.fms(end)])
			xticks([1 2 5 10 20 50 100 200 500])
			set(gca,'fontsize',fontsize)
			ylimit = ylim;
			ylim([0 ylimit(2)])
			grid on
			box on
			ylabel('Avg. Rate (sp/s)')
		end

		% Plot RVF
		params_RVF = data{5,1};
		nexttile()
		if ~isempty(params_RVF)
			plotPhysRVF([], params_RVF, []);
		end

		% Plot NT
		nexttile()
		params_NT = data(target, 2);
		if ~isempty(params_NT{1})
			yyaxis left
			data_NT = analyzeNT(params_NT{1});
			temporal = analyzeNT_Temporal(data_NT);
			hold on
			errorbar(data_NT.pitch_num,data_NT.rate, data_NT.rate_std/...
				sqrt(params_NT{1}.nrep), 'LineWidth',1.5)
			xlabel('F0 (Hz)')
			ylabel('Avg. Rate (sp/s)')
			xlim([55 600])
			xticks([55 110 220 440 600])
			title(extractBefore(params_NT{1}.filename{1}, '.'))
			hold off
			set(gca, 'xscale', 'log')

			% Plot vector strength
			yyaxis right
			hold on
			plot(data_NT.pitch_num, temporal.VS)
			xlabel('F0 (Hz)')
			ylabel('Vector Strength')
			set(gca, 'fontsize', fontsize)
			set(gca, 'xscale', 'log')
			xlim([55 600])
			xticks([55 110 220 440 600])
			ylim([0 1])
			grid on
		end

		% Plot PSTH
		max_rate = max(temporal.PSTH, [], 'all');
		note_values = round(data_NT.pitch_num);
		nexttile([4, 1])
		num_stim = length(note_values);
		hold on
		for j = 1:num_stim

			% Plot PSTHs
			counts = temporal.PSTH(j,:);
			edges = temporal.t;
			t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
			x_patch = repelem(edges, 2);
			y_patch = repelem([0; counts(:); 0]', 2);
			y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
			offset = (j-1)*max_rate; % Adjust offset amount
			patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');

			% Plot smoothed PSTH
			plot(temporal.t(1:end-1), temporal.PSTH_smooth(j,:)+offset,'k', 'LineWidth',1.5);

		end
		ylim([0 max_rate*num_stim])
		xlabel('Time (ms)')
		box on
		yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
		yticklabels(note_values)
		title('PSTH')
		set(gca, 'fontsize', fontsize)

		% Plot period histogram
		nexttile([4, 1])
		max_rate = max(temporal.p_hist, [], 'all');
		for j = 1:num_stim

			% Plot PSTHs
			counts = temporal.p_hist(j,:);
			edges = temporal.t_hist(j,:);
			t_bin = edges(1:end-1) + diff(edges)/2; % Bin centers
			x_patch = repelem(edges, 2);
			y_patch = repelem([0; counts(:); 0]', 2);
			y_patch = y_patch(2:end-1); % Creates [0 y1 y1 0 0 y2 y2 0...]
			offset = (j-1)*max_rate; % Adjust offset amount
			patch(x_patch, y_patch + offset, 'b', 'FaceAlpha',0.5, 'EdgeColor','k');
			period = 1/data_NT.pitch_num(j)*1000;
		end
		ylim([0 max_rate*num_stim])
		xlabel('Time (ms)')
		box on
		yticks(linspace(max_rate/2, max_rate*num_stim-max_rate/2, num_stim))
		yticklabels(note_values)
		xlim([0 17.5])
		%xticks(1:5)
		grid on
		title('Period Histogram')
		set(gca, 'fontsize', fontsize)


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
values = [7.75, 10.5];
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