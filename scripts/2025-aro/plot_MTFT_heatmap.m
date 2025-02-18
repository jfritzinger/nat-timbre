%% plot_MTFT_heatmap
clear 

% Load in session 
rabbit = 27;
sesh = 534; %531;

% TT = 1;
% N = 1;
% CF = 1520;

% TT = 4;
% N = 2;
% CF = 2000;

TT = 3;
N = 1;

% Identify current machine, which changes paths
[userid, base_dir, ~, report_path, data_path] = findPaths();


session = sprintf('R%03dS%03d', rabbit, sesh);
if ispc
    rab_str =  ['R0' num2str(rabbit)];
    base_dir = base_dir{contains(base_dir,rab_str)};
    session_dir = fullfile(base_dir, session);
else
    session_dir = fullfile(base_dir, ['R0' num2str(rabbit)], session);
end

[clusters, params, stims] = loadPhysiologySession(session_dir, session, userid);
cluster = clusters([clusters.tetrode] == TT); % Select tetrode
cluster = cluster([cluster.neuron] == N); % Select neuron

%% Plot RVF and RM 

has_RVF = cellfun(@(p) strcmp(p.type,'RVF'),params);
%plotPhysRVF(cluster, params(has_RVF), stims)

has_RM = cellfun(@(p) strcmp(p.type,'type=RM'),params);
%plotPhysRM(cluster, params(has_RM), stims, CF)

%% Get MTFT rates in matrix, sorted by carrier freq

has_MTFs = cellfun(@(p) strcmp(p.type,'typMTFT'),params);
carrier_freqs = cellfun(@(s) s.carrier_freq, params(has_MTFs));
[sort_freq, order] = sort(carrier_freqs);
num_MTFs = length(carrier_freqs);

MTF_ind = find(has_MTFs);
MTF_ind = MTF_ind(order);

for ii = 1:num_MTFs 
	[BMF, WMF, MTF_shape, fig, data] = plotPhysMTF(cluster, params(MTF_ind(ii)), stims);
	rate_mat(ii,:) = data.rate_sm-data.rate_sm(1);
	rate_mat2(ii,:) = data.rate_sm;
end
close all

%% 

fms = data.fms;
figure
hold on
for ii = 1:length(sort_freq)

	plot3([sort_freq(ii) sort_freq(ii)], [fms(1) fms(end)], ...
		[rate_mat2(ii,1) rate_mat2(ii,1)], 'k', 'LineWidth',2)

	plot3(repmat(sort_freq(ii), length(fms), 1), fms, ...
		rate_mat2(ii,:), 'k', 'LineWidth',2)

	plot3([sort_freq(ii) sort_freq(ii)], [fms(1) fms(1)], ...
		[0 rate_mat2(ii,1)], 'k', 'LineWidth',2)
	plot3([sort_freq(ii) sort_freq(ii)], [fms(end) fms(end)], ...
		[0 rate_mat2(ii,1)], 'k', 'LineWidth',2)

	new_fms = linspace(1.2, 406, 200);
	new_mat = interp1(fms, rate_mat2(ii,:), new_fms, 'linear');

	x = repmat(sort_freq(ii), length(new_fms), 1)';
	y = new_fms;
	z = new_mat;
	pos_mat = z;
	pos_mat(pos_mat<pos_mat(1)) = pos_mat(1);
	neg_mat = z;
	neg_mat(neg_mat>=neg_mat(1)) = neg_mat(1);

	fill3([x x(1)], [y 408], [pos_mat pos_mat(1)], 'r', 'FaceAlpha',0.5)
	fill3([x x(1)], [y 408], [neg_mat neg_mat(1)], 'b', 'FaceAlpha',0.5)
end
ylabel('Mod. Freq. (Hz)');
xlabel('Carrier Freq. (Hz)');
zlabel('Rate (sp/s)');
title('Tone MTFs');
%view([-50 65]) %([-55, 60])
view([-153, 52])
set(gca, 'fontsize', 22)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ydir', 'reverse')
yticks([2 5 10 20 50 100 200 500])
ylim([1.2 408])
xticks([1000 2000 4000 8000])

%% Export 

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/2025-aro';
set(gcf, 'Renderer', 'painters')
print('-dsvg', '-vector', fullfile(savepath,'plot_MTFT_heatmap.svg'))


%% Plot heatmaps 
% 
% figure
% fms = data.fms;
% h = pcolor(fms, sort_freq, rate_mat);
% hold on
% %yline(CF, 'k')
% set(h, 'EdgeColor', 'none');
% set(gca, 'xscale', 'log')
% set(gca, 'yscale', 'log')
% yticks([500 1000 2000 5000 10000])
% colormap(redblue)
% colorbar
% max_rate = max(rate_mat, [], "all");
% clim([-1 1]*max_rate)
% xticks([2, 5, 10, 20, 50, 100, 200, 500 1000 2000 5000 10000])
% ylabel('Tone Frequency (Hz)')
% title('MTFT with noise-alone subtracted')
% xlabel('Mod. Frequency (Hz)')
% set(gca, 'fontsize', 16)


%% Plot 3D
% 
% % carrier freq
% % mod freq
% % rate
% 
% figure
% hold on
% 
% for ii = 1:length(sort_freq)
% 
% 	plot3([sort_freq(ii) sort_freq(ii)], [fms(1) fms(end)], ...
% 		[rate_mat2(ii,1) rate_mat2(ii,1)], 'k', 'LineWidth',2)
% 
% 	plot3(repmat(sort_freq(ii), length(fms), 1), fms, ...
% 		rate_mat2(ii,:), 'k', 'LineWidth',2)
% 
% end
% surf(sort_freq, [2 600], [rate_mat2(:,1) rate_mat2(:,1)]', 'EdgeColor','none', ...
% 	'FaceAlpha',0.5)
% % surf(sort_freq, fms, rate_mat2', 'EdgeColor','none', ...
% % 	'FaceAlpha',0.5)
% set(gca, 'xscale', 'log')
% set(gca, 'yscale', 'log')
% ylim([1.2 600])
% colormap(redblue)

%% 
% 
% new_fms = linspace(1.2, 406, 200);
% new_mat = interp1(fms, rate_mat2(ii,:), new_fms, 'linear');
% 
% x = new_fms;
% y = new_mat;
% pos_mat = y;
% pos_mat(pos_mat<pos_mat(1)) = pos_mat(1);
% neg_mat = y;
% neg_mat(neg_mat>=neg_mat(1)) = neg_mat(1);
% 
% area(x, pos_mat, y(1), 'FaceColor','r', 'FaceAlpha',0.5)
% area(x, neg_mat, y(1), 'FaceColor','b', 'FaceAlpha',0.5)


