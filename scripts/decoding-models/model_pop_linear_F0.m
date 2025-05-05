clear 

%% Load in data 

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT2.mat'), 'nat_data')

%% Get correct output of model 
target = 'Bassoon';
%target = 'Oboe';

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
F0s = log10(F0s);

%% Get data into proper matrix 
type = 'rate'; % VS

% Find all rows with bassoon in them
if strcmp(target, 'Bassoon')
	sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
else
	sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
end
num_data = numel(sesh);

data_mat = NaN(length(F0s)*20, num_data);
for ii = 1:num_data
	if strcmp(type, 'rate')
		if strcmp(target, 'Bassoon')
			%X1 = nat_data(sesh(ii)).bass_raterep';
			%X1 = nat_data(sesh(ii)).bass_norm_rep';
			X1 = nat_data(sesh(ii)).bass_sync_rep';
		else
			X1 = nat_data(sesh(ii)).oboe_raterep';
		end
	else
		if strcmp(target, 'Bassoon')
			X1 = nat_data(sesh(ii)).bass_VSrep;
		else
			X1 = nat_data(sesh(ii)).oboe_VSrep;
		end
	end
	X2 = reshape(X1, [], 1);
	data_mat(:,ii) = X2;
end
data_mat = zscore(data_mat);

% Create array of correct responses
response = reshape(repmat(F0s, 1, 20)', 1, []);

% Create table for model
T = array2table(data_mat);
T.response = response';

%% Split into training and testing data 

for group = 0:39

    % Calculate group range (1-20, 21-40, etc.)
    startIdx = group*20 + 1;
    endIdx = (group+1)*20;
    
    % Random permutation within current group
    shuffled = randperm(20) + startIdx - 1;
    
    % Store indices
    trainStart = group*16 + 1;
    trainEnd = (group+1)*16;
    testStart = group*4 + 1;
    testEnd = (group+1)*4;
    
    trainIndices(trainStart:trainEnd) = shuffled(1:16);
    testIndices(testStart:testEnd) = shuffled(17:20);
end

% Create training and testing sets
T_train = T(trainIndices, :);
T_test = T(testIndices, :);

%% Fit using fitrlinear 


hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus', 'KFold',5);
[Mdl,FitInfo,HyperparameterOptimizationResults] = fitrlinear(T_test, 'response', ...
    'OptimizeHyperparameters','auto', ...
    'HyperparameterOptimizationOptions',hyperopts);


% Mdl = fitrlinear(T, 'response','BetaTolerance',0.0001, ...
% 	'Learner','leastsquares', 'Lambda','auto', 'Solver','lbfgs', ...
% 	'KFold',5, 'CrossVal','on', 'Regularization','ridge');

%%

%Tpred_F0 = kfoldPredict(Mdl);
pred_F0 = predict(Mdl, T_test);
r = corrcoef(pred_F0, T_test.response);
r2 = r(1, 2)^2;

figure('Position',[31,910,910,413])
nexttile
scatter(T_test.response, pred_F0, 'filled', 'MarkerFaceAlpha',0.5)
set(gca, 'xscale', 'log', 'yscale', 'log')
% ylim([57 588])
% xlim([57 588])
hold on
plot(T_test.response, T_test.response)
title(r2)

closest_cat = [];
for ii = 1:length(pred_F0)

	differences = abs(F0s - pred_F0(ii));
	[~, closest_column_index] = min(differences);
	closest_cat(ii) = F0s(closest_column_index);
end


% Plot confusion matrix
nexttile
C = confusionmat(T_test.response, closest_cat);
confusionchart(C)

% Calculate accuracy
chart = confusionchart(T_test.response,closest_cat); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
title(sprintf('Accuracy = %0.2f%%', accuracy*100))


%% Leave-one-out linear regression decoding

% Example data
R = data_mat'; % Neural data: [neurons x elements]
F0 = reshape(repmat(F0s, 1, 20)', 1, []);

% Transpose R for regression: [elements x neurons]
X = R';

% Leave-one-out cross-validation
predF0 = zeros(800, 1);

for v = 1:800
    train_idx = setdiff(1:800, v);
    X_train = X(train_idx, :);
    y_train = F0(train_idx);
    
    % Add intercept
    X_train_aug = [ones(length(train_idx),1), X_train];
    
    % Fit linear regression
    beta = X_train_aug \ y_train';
    
    % Predict for held-out element
    X_test = [1, X(v,:)];
    predF0(v) = X_test * beta;
end

% Compute R^2
F0_mean = mean(F0);
R2 = 1 - sum((F0 - predF0).^2) / sum((F0 - F0_mean).^2);

