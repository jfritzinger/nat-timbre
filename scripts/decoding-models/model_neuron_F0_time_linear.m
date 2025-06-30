%% model_neuron_time_F0_linear
clear

%% Load in data

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')


%% Get correct output of model 
target = 'Bassoon';
%target = 'Oboe';
%target = 'Invariant';

% Get bassoon stimulus
F0s = getF0s(target);
F0s = log10(F0s);

%% Get all rates for each repetition for bassoon (one example neuron)

[sesh, num_data] = getF0Sessions(nat_data, target);
for ind = 1:num_data
	index = sesh(ind);

	T = getF0NeuronTable(nat_data, target, index, F0s, 'Timing');

	for imodelrep = 1 

		% Fit using fitrlinear
		Mdl = fitrlinear(T, 'Response','BetaTolerance',0.0001, ...
			'Learner','leastsquares', 'Lambda','auto', 'Solver','lbfgs', ...
			'KFold',5, 'CrossVal','on', 'Regularization','ridge');
		pred_F0 = kfoldPredict(Mdl);
	end
	
	r = corrcoef(pred_F0, T.Response);
	r2 = r(1, 2)^2;
	actual = T.Response; % Actual response values
	mse = mean((actual - pred_F0).^2);   % Mean Squared Error
	rmse = sqrt(mse);                     % Root Mean Squared Error
	mae = mean(abs(actual - pred_F0));    % Mean Absolute Error

	% figure
	% scatter(10.^(T.Response), 10.^(pred_F0), 'filled', ...
	% 	'MarkerFaceAlpha',0.5, 'MarkerFaceColor','k')
	% set(gca, 'xscale', 'log', 'yscale', 'log')
	% hold on
	% plot(10.^T.Response, 10.^T.Response, 'k')
	% title(['R^2 = ' num2str(r2)])
	% xlabel('Actual F0s (Log10(F0s))')
	% ylabel('Predicted F0s (Log10(F0s))')
	% grid on
	
	% Save data for each
	neuron_time_F0_lin(ind).putative = nat_data(index).putative;
	neuron_time_F0_lin(ind).CF = nat_data(index).CF;
	neuron_time_F0_lin(ind).MTF = nat_data(index).MTF;
	neuron_time_F0_lin(ind).response = T.Response;
	neuron_time_F0_lin(ind).T = T;
	neuron_time_F0_lin(ind).pred_F0 = pred_F0;
	neuron_time_F0_lin(ind).r2 = r2;
	neuron_time_F0_lin(ind).rmse = rmse;
	
	fprintf('%d/%d, %0.2f%% done!\n', ind, num_data, ind/num_data*100)
end

%% Save struct of data 

save(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '_Linear.mat']), ...
	"neuron_time_F0_lin", '-v7.3')

%% Plot overall 

overall_r2 = [neuron_time_F0_lin.r2];
edges = linspace(0, 1, 51);
figure
histogram(overall_r2, edges)

%% Get best units 

[r2_sort,ind_high] = sort(overall_r2, 'descend');
ind = 1;

figure
scatter(10.^(neuron_time_F0_lin(ind_high(ind)).T.Response), 10.^(neuron_time_F0_lin(ind_high(ind)).pred_F0), 'filled', ...
	'MarkerFaceAlpha',0.5, 'MarkerFaceColor','k')
set(gca, 'xscale', 'log', 'yscale', 'log')
hold on
plot(10.^neuron_time_F0_lin(ind_high(ind)).T.Response, 10.^neuron_time_F0_lin(ind_high(ind)).T.Response, 'k')
title(['R^2 = ' num2str(r2_sort(ind))])
xlabel('Actual F0s (Log10(F0s))')
ylabel('Predicted F0s (Log10(F0s))')
grid on
