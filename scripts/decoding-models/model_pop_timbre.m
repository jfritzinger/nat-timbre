%% model_pop_VS_F0
clear

%% Load in data

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT.mat'), 'nat_data')

%% Get correct output of model

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
target = 'Bassoon';
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);

%% Get data into proper matrix
type = 'rate'; % VS

% Create array of correct responses
response = [ones(1, 20) repmat(2, 1, 20)]';

% Find all rows with bassoon and oboe
has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
sesh = find(has_bass & has_oboe);
num_data = numel(sesh);
ind_b = 25:40;
ind_o = [1 3:17];

% Sort data 
for target = 1:16

	data_mat = NaN(2*20, num_data);
	for ii = 1:num_data

		if strcmp(type, 'rate')
			X1 = nat_data(sesh(ii)).bass_raterep(:,ind_b(target));
			X2 = nat_data(sesh(ii)).oboe_raterep(:,ind_o(target));
			X = [X1; X2];
		else
			X1 = nat_data(sesh(ii)).bass_VSrep(ind_b(target),:)';
			X2 = nat_data(sesh(ii)).oboe_VSrep(ind_o(target),:)';
			X = [X1; X2];
			X(X>=0.99) = NaN;
			X(isnan(X)) = 0;
		end
		data_mat(:,ii) = X;
	end

	% Create table for model
	T = array2table(data_mat);
	T.Response = response;
end

%% Run model 







%% Run model

% ncond = 2;
% for target = 1:16
% 
	% for irep = 1:nrep
	% 
	% 	% Take out test data
	% 	ind_test = false(ncond*20, 1); % Preallocate a logical array for 800 elements (40 * 20)
	% 	for istim = 1:ncond
		% 	index = randperm(20, 4); % Randomly select 4 indices from 1 to 20
		% 	ind_test((istim-1)*20 + index) = true; % Set the selected indices to true
	% 	end
	% 
	% 	% Make training and test rows
	% 	test_mat = data_mat(ind_test);
	% 	train_mat = data_mat(~ind_test);
	% 
	% 	% Train model on training data
	% 	% mdl = fitrlinear(train_mat, response,'Regularization','lasso',...
	% 	% 	'Solver','sparsa');
	% 	%mdl = fitclinear(train_mat, response, 'Learner', 'logistic', 'KFold', 5);
	% 	mdl = fitclinear(train_mat, response, 'Learner', 'svm', ...
		% 	'OptimizeHyperparameters', {'Lambda', 'Regularization'}, 'Solver','sparsa', ...
		% 	'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations', 50, ...
		% 	'ShowPlots', false));
	% 	% mdl = fitclinear(X, Y, 'OptimizeHyperparameters', 'auto', ...
	% 	% 	'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations', 100, ...
	% 	% 	'MaxTime', 1200, ...
	% 	% 	'AcquisitionFunctionName', 'expected-improvement-plus', ...
	% 	% 	'ShowPlots', true, ...
	% 	% 	'Verbose', 1));
	% 
	% 	%mdl = fitlm(train_mat, response);
	% 
	% 	% Evaluate
	% 	output = predict(mdl, test_mat);
	% 	output_avg = mean(reshape(output, 4, ncond), 1);
	% 
	% 	% Subtract real - predicted to get accuracy
	% 	accuracy_one(irep, 1) = sum(response_test(1:4)==output(1:4))/4;
	% 	accuracy_one(irep, 2) = sum(response_test(5:8)==output(5:8))/4;
	% 	alloutput_avg(irep, :) = output_avg;
	% 
	% 	% Correlation coefficient
	% 	R_temp = corrcoef(response_test, output);
	% 	R(irep) = R_temp(1, 2);
	% 	R2(irep) = R_temp(1,2)^2;
	% 
	% 	closest(irep, :) = output;
	% 	actual(irep, :) = response_test;
	% end
	% 
	% %% Plot outputs
	% figure('Position',[231,839,1016,351])
	% nexttile 
	% 
	% % Confusion matrix for best model
	% [best_R2, best_ind] = max(R2);
	% C = confusionmat(actual(best_ind, :), closest(best_ind, :));
	% chart = confusionchart(actual(best_ind, :),closest(best_ind, :)); % Generate confusion chart
	% confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
	% accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
	% title(['Accuracy = ' num2str(accuracy)])
	% 
	% nexttile
	% hold on
	% swarmchart(ones(10,1), alloutput_avg(:,1), 'filled')
	% swarmchart(ones(10,1)*2, alloutput_avg(:,2), 'filled')
	% xlim([0.5 2.5])
	% ylim([0.5 2.5])
	% xlabel('Actual Timbre')
	% ylabel('Predicted Timbre')
	% 
	% % Plot accuracy
	% nexttile
	% swarmchart(ones(10,1), accuracy_one(:,1), 'filled')
	% hold on
	% swarmchart(ones(10,1)*2, accuracy_one(:,2), 'filled')
	% xlim([0.5 2.5])
	% %set(gca, 'XScale', 'log')
% 
% end