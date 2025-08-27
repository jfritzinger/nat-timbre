%% model_pop_VS_F0
clear 

%% Load in data 

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
% load(fullfile(base, 'model_comparisons',  'Model_NT2.mat'), 'nat_model')
% nat_data = nat_model;

%% Get correct output of model 
%target = 'Bassoon';
%target = 'Oboe';
target = 'Invariant';

% Get stimulus
F0s = getF0s(target);
[sesh, num_data] = getF0Sessions(nat_data, target);
T = getF0PopTable(nat_data, target, sesh, F0s, num_data, 'classification', 'Rate');


%% Run model 

nrep = 1;
best_accuracy = -inf;
for irep = 1:nrep

	% Separate out training and testing data
	[T_train, T_test] = splitData_Reps(T); 

	[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
		trainClassifierPopRateF0(T_train, F0s);
	C = confusionmat(T_train.Response, validationPredictions1);
	accuracy = sum(diag(C))/sum(C, 'all');

	if accuracy > best_accuracy
		trainedClassifier = trainedClassifier1;
		validationAccuracy = validationAccuracy1;
		validationPredictions = validationPredictions1;
		ibest = irep;
		best_accuracy = accuracy;
	end

	% Print out progress
	fprintf('%d/%d, %0.2f%% done!\n', irep, nrep, irep/nrep*100)

	% To make predictions with the returned 'trainedClassifier' on new data
	[yfit,scores] = trainedClassifier.predictFcn(T_test);
	C = confusionmat(T_test.Response, yfit);
	accuracy_test(irep) = sum(diag(C))/sum(C, 'all');
end

% Mdl = trainedClassifier.ClassificationLinear; % used to be Linear
% imp = permutationImportance(Mdl);

%% Save model output 

figure
confusionchart(C);
accuracy = sum(diag(C)) / sum(C(:)); % Calculate accuracy
title([target ' ' num2str(accuracy)])

pop_rate_F0.trainedClassifier = trainedClassifier;
pop_rate_F0.validationPredictions = validationPredictions;
pop_rate_F0.response = T.Response;
pop_rate_F0.T = T;
pop_rate_F0.C = C;
pop_rate_F0.validationAccuracy = validationAccuracy;
pop_rate_F0.putative = {nat_data(sesh).putative};
pop_rate_F0.CF = [nat_data(sesh).CF];
pop_rate_F0.MTF = {nat_data(sesh).MTF};
if strcmp(target, 'Bassoon')
	pop_rate_F0.rate = {nat_data(sesh).bass_rate};
	pop_rate_F0.rate_std = {nat_data(sesh).bass_rate_std};
elseif strcmp(target, 'Oboe')
	pop_rate_F0.rate = {nat_data(sesh).oboe_rate};
	pop_rate_F0.rate_std = {nat_data(sesh).oboe_rate_std};
else
	pop_rate_F0.oboe_rate = [nat_data(sesh).oboe_rate];
	pop_rate_F0.oboe_rate_std = [nat_data(sesh).oboe_rate_std];
	pop_rate_F0.bass_rate = [nat_data(sesh).bass_rate];
	pop_rate_F0.bass_rate_std = [nat_data(sesh).bass_rate_std];
end
% pop_rate_F0.imp = imp;

%% Get shuffled accuracy 
% 
% rng("shuffle"); % Seed based on current time
% nreps = 100;
% data_mat = table2array(T(:,1:end-1));
% for imodel = 1:nreps
% 
% 	% Shuffle data (rows)
% 	data_mat3 = zeros(length(F0s)*20, num_data);
% 	for ind = 1:length(F0s)*20
% 		row = data_mat(ind, :);
% 		shuffled_data = row(randperm(numel(row)));
% 		data_mat3(ind, :) = shuffled_data;
% 	end
% 	data_mat = data_mat3;
% 
% 	% Put data into table
% 	T_new = array2table(data_mat);
% 	T_new.Response = T.Response;
% 
% 	% Run model with kfold validation
% 	[trainedClassifier1, validationAccuracy1, validationPredictions1] = ...
% 		trainClassifierPopRateF0(T_new, target);
% 	shuffled_accuracy(imodel) = validationAccuracy1;
% 
% 	% Print out progress
% 	fprintf('%d/%d, %0.2f%% done!\n', imodel, nreps, imodel/nreps*100)
% end
% pop_rate_F0.shuffled_accuracy = shuffled_accuracy;

%% Plot outputs 

save(fullfile(base, 'model_comparisons', ...
	['Pop_Rate_F0_' target '2.mat']), "pop_rate_F0")
% save(fullfile(base, 'model_comparisons', ...
% 	['Model_Pop_Rate_F0_' target '2.mat']), "pop_rate_F0")

%% 

% target = 'Invariant';
% load(fullfile(base, 'model_comparisons', ...
% 	['Model_Pop_Rate_F0_' target '.mat']), "pop_rate_F0")
% figure
% confusionchart(pop_rate_F0.C);
% accuracy = sum(diag(pop_rate_F0.C)) / sum(pop_rate_F0.C(:)); % Calculate accuracy
% title([target ' ' num2str(accuracy)])