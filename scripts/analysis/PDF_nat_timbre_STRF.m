%% PDF_nat_timbre_STRF

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
filename = 'NT_STRF_Comparison';
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

%% Plot each dataset 

% Find sessions for target synthetic timbre response
bin200(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);
bin200(:,2) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
has_data = bin200(:,1) | bin200(:,2);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);

% Plot each neuron
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

		% Set up figure
		fig = figure('Position',[292,274,1264,420]);
		tiledlayout(4, 2, 'TileSpacing','compact', 'Padding','tight')
		x_label = [1000 2000 3000 4000 6000 8000]./1000;
		fontsize = 10;

		% Plot RM
		params_RM = data{2, 2};
		nexttile
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
		nexttile
		if ~isempty(params_MTF)
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

		% Plot NT
		params_NT = data(13:14, 2);
		data_NT = cell(2, 1);
		for ii = 1:2
			nexttile
			if ~isempty(params_NT{ii})
				data_NT = analyzeNT(params_NT{ii});
				bar(data_NT.pitch, data_NT.rate)
				hold on
				errorbar(data_NT.pitch,data_NT.rate, data_NT.rate_std/sqrt(params_NT{ii}.nrep), 'LineStyle', 'none', 'Color', [0 0 0])
				xlabel('Pitch (Hz)')
				title(extractBefore(params_NT{ii}.filename{1}, '.'))
				hold off
			end
		end

		% Plot Natural timbre STRF
		for ii = 1:2
			if ~isempty(params_NT{ii})
				[~, ~, ~] = plotNatSTRF([], params_NT{ii}, []);
			else
				nexttile
			end
		end

		% Plot STRF 
		params_STRF = data{4, 2};
		if ~isempty(params_STRF)
			[params_STRF, ~, ~] = plotPhysSTRF([], params_STRF, []);
		else
			nexttile
			title('Noise STRF')
		end

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
values = [8, 10];
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

function [params, fig, data] = plotNatSTRF(cluster, params, stims)

% Process cluster
if iscell(params)
	params = params{1}; % If multiple, plots only the first instance
else
	params = params;
end
ds = params.dsid;
dsi = 1;
if isempty(stims)
	this_ds = params.stims.dsid == ds;
	cluster = params.cluster;
	stim = params.stims;
else
	if length(stims)==1
		this_ds = stims.dsid == ds(dsi);
		stim = stims;
	else
		this_ds = stims(ds).dsid == ds(dsi);
		stim = stims(ds);
	end
	params.plot_type = 'STRF';

end

% Regenerate stimuli
fs = stim.fs;
params.Fs = fs;
params.mnrep = 1;
[params] = generate_NT(params);
noises = 1:size(params.stim,1);
noise_mat = params.stim';
num_noises = length(noises);
xpts = floor(params.dur*fs);
dur = params.dur;
sc = 20e-6 * power(10,params.spl/20);

% This method uses smoothed spike time
win = 0.02;
T_pts = win*fs;
h0 = zeros(1,num_noises);
h1 = zeros(num_noises,T_pts);
h2 = zeros(num_noises,T_pts,T_pts);
stim_mx = zeros(T_pts,xpts - T_pts);
window = gausswin(ceil(1e-3*fs),1.5); % 3 std

for istim = 1:num_noises
	stimi = noise_mat(:,istim).';

	bin_width = 1/fs*1e6; % convert to microsec
	stim_time_ind = find([params.list.iwav_file] == istim);
	stim_time_onset = stim.times(stim_time_ind);
	stim_time_offset = stim_time_onset+dur*1e6;
	psth_edges = 0:bin_width:dur*1e6;
	spktrains = zeros(length(psth_edges),length(stim_time_ind));
	for itime = 1:length(stim_time_ind)
		spk_ind_this = find(cluster.t_spike >= stim_time_onset(itime) & ...
			cluster.t_spike <= stim_time_offset(itime));
		spk_time = cluster.t_spike(spk_ind_this) - stim_time_onset(itime);
		spktrains(:,itime) = histc(spk_time,psth_edges);
	end
	response = sum(spktrains,2).';
	response = conv(response(T_pts + 1:end-1),window,'same');

	for T = 0:T_pts - 1
		stim_mx(T + 1,:) = stimi((T_pts + 1:end) - T);
	end
	h0(istim) = mean(response);
	h1(istim,:) = response*stim_mx';
	h2(istim,:,:) = bsxfun(@times,response - mean(response),stim_mx)*stim_mx';
end

H0 = mean(h0);
H1 = mean(h1)./sc./length(response);
H2 = reshape(mean(h2),T_pts,T_pts)./(2*sc.^2)./length(response);

if sum(H1)==0
	fig = figure;

	% Create struct that contains processed data
	data.t = [];
	data.f = [];
	data.H2ex_strf = [];
	data.H2in_strf = [];
	data.clims_strf = [];
	data.tlims = [];
	data.flims = [];
	data.strf = [];
else

	%quickplot_wk(fs, T_pts, H1, H2)
	fn = fs/2;

	tlims = [0 T_pts/fs];
	t = (0:T_pts - 1)/fs;
	f = fn*linspace(0,1,T_pts/2 + 1);

	% Process kernels and plot
	[U,S,V] = svd(H2);
	k = sign(diag(U).*diag(V)).*diag(S);	% weights for the singular vectors (see Lewis et al. 2002)
	negInd = find(k<0);						% column indices of the negatively-weighted singular vectors
	posInd = find(k>=0);					% column indices of the positively weighted singular vectors
	U = repmat(abs(k'),T_pts,1).*U;			% The weighted vectors
	U_fft = 2*abs(fft(U)/T_pts);
	U_fft = U_fft(1:T_pts/2 + 1,:);
	H2ex = U(:,posInd)*V(:,posInd)';
	H2in = U(:,negInd)*V(:,negInd)';
	H2 = H2ex + H2in;
	env = abs(hilbert(U));

	data.U = U;
	data.S = S;
	data.V = V;
	data.H0 = H0;
	data.H1 = H1;
	data.H2 = H2;
	data.T_pts = T_pts;
	data.fs = fs;

	negInd2 = negInd(1:4);
	posInd2 = posInd(1:4);

	H1_fft = 2*abs(fft(H1)/T_pts);
	H1_fft = H1_fft(1:T_pts/2 + 1);

	H2ex_strf = 2*abs(fft(H2ex)/T_pts);
	H2ex_strf = H2ex_strf(1:T_pts/2 + 1,:);
	H2ex_strf_null = median(reshape(H2ex_strf(:,t<0.0025),1,[]));
	H2ex_strf = H2ex_strf - H2ex_strf_null;

	H2in_strf = 2*abs(fft(H2in)/T_pts);
	H2in_strf = H2in_strf(1:T_pts/2 + 1,:);
	H2in_strf_null = median(reshape(H2in_strf(:,t<0.0025),1,[]));
	H2in_strf = H2in_strf - H2in_strf_null;

	clims_strf = max(max(abs([H2ex_strf,H2in_strf])))*[-1 1];

	% Get x limits
	pos_max = [db(U_fft(:,posInd2)) db(U_fft(:,negInd2))];
	yaxes_lim = -40 + max(max(db(U_fft)));
	[ind, ~] = find(pos_max > yaxes_lim);
	max_ind = max(ind);
	f_max = f(max_ind);
	flims = [0 f_max+50];

	% Plots
	fig = 0;
	nexttile
	axx = imagesc(t,f,H2ex_strf-H2in_strf,clims_strf);
	set(gca,'Ydir','normal','XLim',tlims,'YLim',flims)
	colormap(redblue)
	grid on
	title([params.target ' STRF']);
	xlabel('Time (s)');
	ylabel('Frequency (Hz)')
	axx.HitTest = 'off';
	set(gca,'FontSize',7)

	% Create struct that contains processed data
	data.t = t;
	data.f = f;
	data.H2ex_strf = H2ex_strf;
	data.H2in_strf = H2in_strf;
	data.clims_strf = clims_strf;
	data.tlims = tlims;
	data.flims = flims;
	data.strf = H2ex_strf-H2in_strf;
end
end