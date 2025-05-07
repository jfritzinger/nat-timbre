%% model_pop_VS_F0
clear 

%% Load in data 

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT2.mat'), 'nat_data')

%% Get correct output of model 
%target = 'Bassoon';
target = 'Oboe';

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning

% Get bassoon stimulus
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s1 = tuning.Frequency(index);
[F0s, order] = sort(F0s1);

%% Get data into proper matrix 

% Find all rows with bassoon in them
if strcmp(target, 'Bassoon') % Bassoon
	sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
	data_mat = NaN(length(F0s)*20, num_data);
	for ii = 1:num_data
		X1 = nat_data(sesh(ii)).bass_raterep';
		X2 = reshape(X1, [], 1);
		data_mat(:,ii) = X2;
	end
else % Oboe 
	sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
	data_mat = NaN(length(F0s)*20, num_data);
	for ii = 1:num_data
		X1 = nat_data(sesh(ii)).oboe_raterep';
		X2 = reshape(X1, [], 1);
		data_mat(:,ii) = X2;
	end
end
num_data = numel(sesh);

% Create array of correct responses
response = reshape(repmat(F0s, 1, 20)', 1, []);

% Create table for model
T = array2table(data_mat);
T.Response = response';

% % Use only first 20 F0s 
% T_new = T(1:400, :);

% % Use only odd F0d
% idx = [];
% for k = 0:19
%     idx = [idx, (1 + 40*k):(20 + 40*k)];
% end
% T_new2 = T(idx,:);


%% Run model 

[trainedClassifier, validationAccuracy, validationPredictions] = ...
	trainClassifierPopRateF0(T, F0s);
C = confusionmat(response, validationPredictions);

%% Save model output 

pop_rate_F0.trainedClassifier = trainedClassifier;
pop_rate_F0.validationPredictions = validationPredictions;
pop_rate_F0.response = response;
pop_rate_F0.T = T;
pop_rate_F0.C = C;
pop_rate_F0.validationAccuracy = validationAccuracy;
pop_rate_F0.putative = nat_data(sesh).putative;
pop_rate_F0.CF = nat_data(sesh).CF;
pop_rate_F0.MTF = nat_data(sesh).MTF;
pop_rate_F0.rate = nat_data(sesh).oboe_rate;
pop_rate_F0.rate_std = nat_data(sesh).oboe_rate_std;

save(['Pop_Rate_F0_' target '.mat'], "pop_rate_F0")


%% Run model
% 
% nrep = 5;
% for irep = 1:nrep
% 
% 	% Take out test data
% 	ind_test = false(length(F0s)*20, 1); % Preallocate a logical array for 800 elements (40 * 20)
% 	for istim = 1:length(F0s)
% 		index = randperm(20, 4); % Randomly select 4 indices from 1 to 20
% 		ind_test((istim-1)*20 + index) = true; % Set the selected indices to true
% 	end
% 
% 	% Make training and test rows
% 	test_mat = data_mat(ind_test);
% 	train_mat = data_mat(~ind_test);
% 
% 	% Train model on training data
% 	% mdl = fitrlinear(train_mat, response,'Regularization','lasso',...
% 	% 	'Solver','sparsa', 'Lambda','auto');
% 	%mdl = fitlm(train_mat, response);
% 	% mdl = fitrlinear(train_mat, response, 'Regularization', 'lasso', ...
% 	% 	'Lambda', 0.01, 'KFold', 5); % Linear discriminant analysis
% 	mdl = fitcecoc(train_mat, response, ...
% 		'OptimizeHyperparameters','auto'); % Multiclass SVM using Error-Correcting Output Codes
% 
% 	% If the loss is very high (close to random guessing), 
% 	% revisit your data preprocessing steps.	
% 	%cvMdl = fitcecoc(train_mat, response, 'KFold', 5);
% 	% cvLoss = kfoldLoss(mdl); % Cross-validation loss
% 	% disp(['Cross-validation loss: ', num2str(cvLoss)]);
% 
% 	% Evaluate
% 	% trainedModel = mdl.Trained{1}; % Access the first fold's trained model
% 	output = predict(mdl, test_mat);
% 	output_avg = mean(reshape(output, length(F0s), 4), 2);
% 
% 	% Subtract real - predicted to get accuracy
% 	accuracy2(irep,:) = sum(output == response_test') / length(response_test);
% 	accuracy(irep,:) = F0s - output_avg;
% 	alloutput_avg(irep, :) = output_avg;
% 
% 	% Correlation coefficient
% 	R_temp = corrcoef(response_test, output);
% 	R(irep) = R_temp(1, 2);
% 	R2(irep) = R_temp(1,2)^2;
% 
% end
% 
% %% Plot outputs 
% 
% % Confusion matrix for best model 
% [best_R2, best_ind] = max(R2);
% 
% output_hz = 10.^output_avg;
% figure('Position',[231,839,1016,351])
% nexttile
% scatter(F0s1, output_hz)
% hold on
% plot([1 2000], [1 2000], 'k')
% xlim([min(F0s1), max(F0s1)])
% ylim([min(F0s1), max(F0s1)])
% set(gca, 'XScale', 'log', 'YScale', 'log')
% xlabel('Actual F0 (Hz)')
% ylabel('Predicted F0 (Hz)')
% 
% % Plot accuracy 
% nexttile
% scatter(F0s, accuracy, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
% %set(gca, 'XScale', 'log')
% 
% 
% nexttile
% histogram(accuracy2*100)
% xlim([0 50])