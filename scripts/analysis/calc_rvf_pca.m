%% calc_rvf_pca

%% Find natural timbre bassoon datasets, CF, and MTF

% Load in spreadsheet
addpath('/Users/jfritzinger/Projects/nat-timbre/scripts/helper-functions')
[base, datapath, savepath, ppi] = getPathsNT();
modelpath = '/Volumes/Nat-Timbre/data/manuscript';
sessions = readtable(fullfile(base, 'data-cleaning', 'Data_Table.xlsx'),...
	'PreserveVariableNames',true);

%% Create matrices for bassoon and oboe separately

% Natural timbre datasets
RVF_datasets = cellfun(@(s) contains(s, 'R'), sessions.RVF);
RVF_list = find(RVF_datasets);
num_sesh = length(RVF_list);

%% Load in all data and set up matrix
% "RVFs were normalized by their individual peak rates before PCA
% analysis."

iii = 1;
for ii = 1:num_sesh

	% Load in data
	putative = sessions.Putative_Units{RVF_list(ii)};
	CF = sessions.CF(RVF_list(ii));
	MTF_shape = sessions.MTF{RVF_list(ii)};
	load(fullfile(datapath, [putative '.mat']))

	RVF_params = data{5, 1};
	data_RVF = analyzeRVF(RVF_params);

	% Order properly
	[ordered, order] = sort(data_RVF.velocities);
	rates = data_RVF.rate(order);

	% Normalize to max rate of 1
	rates = rates/max(rates);

	if size(rates, 2)==30
		RVF_mat(ii,:) = rates;
		putative_list{ii} = putative;
	else
		error_sesh(iii) = ii;
		iii = iii+1;
	end
end
velocity = ordered;

RVF_mat(error_sesh,:) = [];
putative_list(error_sesh) = [];

%% PCA RVF analysis

allTuningMats = RVF_mat';
[height, num_mats] = size(allTuningMats); %measure sizes of image stack
vectMats = reshape(allTuningMats, height, num_mats)';  % vectorize each tuning matrix
[coeff, score, latent] = pca(vectMats); %run PCA

% Plot variance explained
figure;
explained = 100*latent/sum(latent);
plot(cumsum(explained),'k-'); hold on;
plot(cumsum(explained),'k.'); hold off;
xlabel('Number of Components');
ylabel('Variance Explained (%)');
title('Cumulative Variance Explained');
grid on

% Visualize top 4 components
num_components_to_show = 4; 
figure; tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact');
for n = 1:num_components_to_show
    nexttile(n);
    pcMat = coeff(:,n);
    plot(velocity, pcMat); %clim([-0.5 0.8]);
    title(['PC ' num2str(n)]);
end

% Find images that are most aligned with each principal component
figure; tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact');
for n = 1:num_components_to_show
    [~, max_idx] = max(abs(score(:,n)));
    nexttile(n)
    orig_image = allTuningMats(:,max_idx);
    plot(velocity, orig_image); %clim([-4 8]);
    title(['Neuron most aligned with PC ' num2str(n)],'FontSize',8);
end

%% Save score for PC2 somehow to then add to decoder analysis 

% Get scores for PC2
PC2_score = score(:,2);
[~,best] = sort(abs(PC2_score), 'descend' );
best_ind = best(1:16);

% Plot 
figure; tiledlayout(4,4, 'TileSpacing','compact', 'Padding','compact');
for n = 1:16
    nexttile(n)
    orig_image = allTuningMats(:,best_ind(n));
    plot(velocity, orig_image); %clim([-4 8]);
    title(['Neuron most aligned with PC ' num2str(n)],'FontSize',8);
end

%% 

% Create table
RVF_sessions = table(putative_list', PC2_score);
save(fullfile(base, 'RVF_PC2.mat'), "RVF_sessions")

