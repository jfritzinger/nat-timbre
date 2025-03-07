%% plot_MTFT_heatmap
clear 


%% MTFT SESSIONS 
% 4 sessions
rabbit = 27;

neuron = 8; % Example we were looking at

switch neuron
	% SESSION 1
	case 1 % R027S531
		session = 'R027S531';
		TT = 1;
		N = 1;
		CF = 1520;
	case 2 % R027S531
		session = 'R027S531';
		TT = 1;
		N = 2;
		CF = 1720;
	case 3 % R027S531
		session = 'R027S531';
		TT = 2;
		N = 1;
		CF = 3550;
	case 4 % R027S531
		session = 'R027S531';
		TT = 4;
		N = 2;
		CF = 2000;
	
	% SESSION 2
	case 5 % R027S534
		session = 'R027S534';
		TT = 2;
		N = 1;
		CF = 9190;
	case 6 % R027S534
		session = 'R027S534';
		TT = 3;
		N = 1;
		CF = 5300;
	case 7 % R027S534
		session = 'R027S534';
		TT = 3;
		N = 3;
		CF = 5280;

	% SESSION 3
	case 8 % R027S542
		session = 'R027S542';
		TT = 1;
		N = 1;
		CF = 3500;
	case 9
		session = 'R027S542';
		TT = 2;
		N = 1;
		CF = 9190;
	case 10
		session = 'R027S542';
		TT = 3;
		N = 1;
		CF = 6030;

	% SESSION 4
	case 11 % R027S545
		session = 'R027S545';
		TT = 1;
		N = 1;
		CF = 4040;
	case 12 % R027S545
		session = 'R027S545';
		TT = 2;
		N = 1;
		CF = 9240;
	case 13
		session = 'R027S545';
		TT = 3;
		N = 1;
		CF = 6960;
end


%% Load in data 

% Identify current machine, which changes paths
[userid, base_dir, ~, report_path, data_path] = findPaths();
if ispc
    rab_str =  ['R0' num2str(rabbit)];
    base_dir = base_dir{contains(base_dir,rab_str)};
    session_dir = fullfile(base_dir, session);
else
    session_dir = fullfile(base_dir, ['R0' num2str(rabbit)], session);
end

% Load in physio session
[clusters, params, stims] = loadPhysiologySession(session_dir, session, userid);
cluster = clusters([clusters.tetrode] == TT); % Select tetrode
cluster = cluster([cluster.neuron] == N); % Select neuron

%% Get MTFT rates in matrix, sorted by carrier freq

% Find datasets with MTF and sort based on carrier frequency 
has_MTFs = cellfun(@(p) strcmp(p.type,'typMTFT'),params);
carrier_freqs = cellfun(@(s) s.carrier_freq, params(has_MTFs));
[sort_freq, order] = sort(carrier_freqs);
num_MTFs = length(carrier_freqs);

MTF_ind = find(has_MTFs);
MTF_ind = MTF_ind(order);
data = cell(num_MTFs, 1);
for ii = 1:num_MTFs

	param = params{MTF_ind(ii)}; % If multiple, plots only the first instance
	ds = param.dsid;
	stim = stims(ds);
	this_ds = stim.dsid == ds;

	if param.dur > 10
		dur = param.dur/1000; % stimulus duration in seconds.
	else
		dur = param.dur;
	end

	all_mod_depths = double([param.list.mdepth]).';
	all_mod_depths(all_mod_depths < -80) = max(all_mod_depths);
	[~,~,mdi] = unique(all_mod_depths);
	[fms,~,fmi] = unique(double([param.list.fm]).');
	if fms(1) == 0
		fms(1) = 1.2;
	end

	num_mod_freqs = length(fms);
	num_depths = length(param.all_mdepths);
	map_size = [num_mod_freqs, num_depths];

	spike_rates = cluster.num_spikes_delayed(this_ds)/...
		(dur - param.onsetWin/1000);

	[rate,rate_std,~,rlb, rub] = accumstats({fmi,mdi},spike_rates, map_size);
	rate_sm = smooth_rates(rate,rlb,rub, []);

	[BMF,WMF,MTF_shape, at_100, at_200] = MTFclassification(spike_rates,fms, fmi);

	% Create struct that contains processed data
	data{ii}.rate = rate;
	data{ii}.fms = fms;
	data{ii}.rate_sm = rate_sm;
	data{ii}.rate_std = rate_std;
	data{ii}.MTF_shape = MTF_shape;
	data{ii}.BMF = BMF;
	data{ii}.WMF = WMF;
end

%% Add rates from each dataset into a matrix

fms = data{1}.fms;
num_freq = length(fms);
rate_mat = zeros(num_MTFs, num_freq);
rate_sub_mat = zeros(num_MTFs, num_freq);
rate_sm_mat = zeros(num_MTFs, num_freq);
rate_sm_sub_mat = zeros(num_MTFs, num_freq);
for ii = 1:num_MTFs

	rate_mat(ii,:) = data{ii}.rate;
	rate_sub_mat(ii,:) = data{ii}.rate-data{ii}.rate(1);

	rate_sm_mat(ii,:) = data{ii}.rate_sm;
	rate_sm_sub_mat(ii,:) = data{ii}.rate_sm-data{ii}.rate_sm(1);
end

%% Plot 2D heatmaps 

figure('Position',[276,473,995,402])
tiledlayout(1, 2)

for iplot = 1:2
	nexttile
	fms = data{1}.fms;
	if iplot == 1
		h = pcolor(fms, sort_freq, rate_sm_mat);
		max_rate = max(rate_mat, [], "all");
		title('MTFT raw rates')
		clim([0 1]*max_rate)
	else
		h = pcolor(fms, sort_freq, rate_sm_sub_mat);
		max_rate = max(rate_sub_mat, [], "all");
		title('MTFT with noise-alone subtracted')
		clim([-1 1]*max_rate)
	end
	hold on
	yline(CF, 'k', 'linewidth', 3)
	set(h, 'EdgeColor', 'none');
	set(gca, 'xscale', 'log')
	set(gca, 'yscale', 'log')
	yticks([500 1000 2000 5000 10000])
	%colormap(redblue)
	colormap("parula")
	%shading interp
	
	colorbar
	xticks([2, 5, 10, 20, 50, 100, 200, 500 1000 2000 5000 10000])
	ylabel('Tone Frequency (Hz)')
	xlabel('Mod. Frequency (Hz)')
	set(gca, 'fontsize', 16)
end

%% Plot 3D "neural fluctuation" plots 
% 
% rate_plot = rate_sm_mat; % rate_sm_mat
% figure
% hold on
% for ii = 1:length(sort_freq)
% 
% 	% Plot noise-alone rate 
% 	plot3([sort_freq(ii) sort_freq(ii)], [fms(1) fms(end)], ...
% 		[rate_plot(ii,1) rate_plot(ii,1)], 'k', 'LineWidth',2)
% 
% 	% Plot rates 
% 	plot3(repmat(sort_freq(ii), length(fms), 1), fms, ...
% 		rate_plot(ii,:), 'k', 'LineWidth',2)
% 
% 	% Plot verical line
% 	plot3([sort_freq(ii) sort_freq(ii)], [fms(1) fms(1)], ...
% 		[0 rate_plot(ii,1)], 'k', 'LineWidth',2)
% 	plot3([sort_freq(ii) sort_freq(ii)], [fms(end) fms(end)], ...
% 		[0 rate_plot(ii,1)], 'k', 'LineWidth',2)
% 
% 	% Plot CF
% 	plot3([CF CF], [fms(1) fms(end)], [0 0], 'k', 'LineWidth',2, 'Color','r')
% 
% 	% Interpolate to get better zero crossings for excited/suppressed
% 	new_fms = linspace(1.2, 406, 1000);
% 	new_mat = interp1(fms, rate_plot(ii,:), new_fms, 'linear');
% 	x = repmat(sort_freq(ii), length(new_fms), 1)';
% 	y = new_fms;
% 	z = new_mat;
% 	pos_mat = z;
% 	pos_mat(pos_mat<pos_mat(1)) = pos_mat(1);
% 	neg_mat = z;
% 	neg_mat(neg_mat>=neg_mat(1)) = neg_mat(1);
% 
% 	% Plot enhanced and suppressed regions in red and blue
% 	fill3([x x(1)], [y 408], [pos_mat pos_mat(1)], 'r', 'FaceAlpha',0.5)
% 	fill3([x x(1)], [y 408], [neg_mat neg_mat(1)], 'b', 'FaceAlpha',0.5)
% end
% ylabel('Mod. Freq. (Hz)');
% xlabel('Carrier Freq. (Hz)');
% zlabel('Rate (sp/s)');
% title('Tone MTFs');
% view([-153, 52])
% set(gca, 'fontsize', 22)
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca, 'ydir', 'reverse')
% yticks([2 5 10 20 50 100 200 500])
% ylim([1.2 408])
% xticks([1000 2000 4000 8000])
% grid on 

