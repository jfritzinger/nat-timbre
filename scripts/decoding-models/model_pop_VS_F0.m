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
[F0s1, order] = sort(F0s1);


F0s = log10(F0s1);
response = reshape(repmat(F0s, 1, 16)', 1, []);
response_test = reshape(repmat(F0s, 1, 4)', 1, []);
response2 = reshape(repmat(F0s, 1, 20)', 1, []);

%% Get data into proper matrix 

% Find all rows with bassoon in them
sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
num_data = numel(sesh);

data_mat = NaN(length(F0s)*20, num_data);
for ii = 1:num_data

	%try1 = nat_data(sesh(ii)).bass_VSrep;
	%try1 = nat_data(sesh(ii)).bass_raterep';
	try1 = nat_data(sesh(ii)).bass_VS_CF;
	try2 = reshape(try1, [], 1);
	data_mat(:,ii) = try2;
end

% [coeff, score] = pca(data_mat);
% XReduced = score(:, 1:50); 

T = array2table(data_mat);
T.Response = response2';

%% Run model

nrep = 5;
for irep = 1:nrep

	% Take out test data
	ind_test = false(length(F0s)*20, 1); % Preallocate a logical array for 800 elements (40 * 20)
	for istim = 1:length(F0s)
		index = randperm(20, 4); % Randomly select 4 indices from 1 to 20
		ind_test((istim-1)*20 + index) = true; % Set the selected indices to true
	end

	% Make training and test rows
	test_mat = data_mat(ind_test);
	train_mat = data_mat(~ind_test);

	% Train model on training data
	% mdl = fitrlinear(train_mat, response,'Regularization','lasso',...
	% 	'Solver','sparsa', 'Lambda','auto');
	%mdl = fitlm(train_mat, response);
	% mdl = fitrlinear(train_mat, response, 'Regularization', 'lasso', ...
	% 	'Lambda', 0.01, 'KFold', 5); % Linear discriminant analysis
	mdl = fitcecoc(train_mat, response, ...
		'OptimizeHyperparameters','auto'); % Multiclass SVM using Error-Correcting Output Codes

	% If the loss is very high (close to random guessing), 
	% revisit your data preprocessing steps.	
	%cvMdl = fitcecoc(train_mat, response, 'KFold', 5);
	% cvLoss = kfoldLoss(mdl); % Cross-validation loss
	% disp(['Cross-validation loss: ', num2str(cvLoss)]);

	% Evaluate
	% trainedModel = mdl.Trained{1}; % Access the first fold's trained model
	output = predict(mdl, test_mat);
	output_avg = mean(reshape(output, length(F0s), 4), 2);

	% Subtract real - predicted to get accuracy
	accuracy2(irep,:) = sum(output == response_test') / length(response_test);
	accuracy(irep,:) = F0s - output_avg;
	alloutput_avg(irep, :) = output_avg;

	% Correlation coefficient
	R_temp = corrcoef(response_test, output);
	R(irep) = R_temp(1, 2);
	R2(irep) = R_temp(1,2)^2;

end

%% Plot outputs 

% Confusion matrix for best model 
[best_R2, best_ind] = max(R2);

output_hz = 10.^output_avg;
figure('Position',[231,839,1016,351])
nexttile
scatter(F0s1, output_hz)
hold on
plot([1 2000], [1 2000], 'k')
xlim([min(F0s1), max(F0s1)])
ylim([min(F0s1), max(F0s1)])
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Actual F0 (Hz)')
ylabel('Predicted F0 (Hz)')

% Plot accuracy 
nexttile
scatter(F0s, accuracy, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
%set(gca, 'XScale', 'log')


nexttile
histogram(accuracy2*100)
xlim([0 50])