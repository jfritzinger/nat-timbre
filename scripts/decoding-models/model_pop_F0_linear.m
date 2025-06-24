clear 

%% Load in data 

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT_3.mat'), 'nat_data')

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

for group = 0:34 %0:39

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


% hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus', 'KFold',5);
% [Mdl,FitInfo,HyperparameterOptimizationResults] = fitrlinear(T_train, 'response', ...
%     'OptimizeHyperparameters','auto', ...
%     'HyperparameterOptimizationOptions',hyperopts);
%pred_F0 = predict(Mdl, T_test);


Mdl = fitrlinear(T, 'response','BetaTolerance',0.0001, ...
	'Learner','leastsquares', 'Lambda','auto', 'Solver','lbfgs', ...
	'KFold',5, 'CrossVal','on', 'Regularization','ridge');
pred_F0 = kfoldPredict(Mdl);

% [trainedModel, validationRMSE, pred_F0] =...
% 	trainRegressionModelPopTimeF0(T);

%%

r = corrcoef(pred_F0, T.response);
r2 = r(1, 2)^2;

figure('Position',[31,910,910,413])
nexttile
scatter(T.response, pred_F0, 'filled', 'MarkerFaceAlpha',0.5)
set(gca, 'xscale', 'log', 'yscale', 'log')
% ylim([57 588])
% xlim([57 588])
hold on
plot(T_test.response, T_test.response)
title(['R^2 = ' num2str(r2)])
xlabel('Actual F0s (Log10(F0s))')
ylabel('Predicted F0s (Log10(F0s))')

closest_cat = [];
for ii = 1:length(pred_F0)

	differences = abs(F0s - pred_F0(ii));
	[~, closest_column_index] = min(differences);
	closest_cat(ii) = F0s(closest_column_index);
end

% Plot confusion matrix
nexttile
C = confusionmat(T.response, closest_cat);
confusionchart(C)

% Calculate accuracy
chart = confusionchart(T.response,closest_cat); % Generate confusion chart
confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
title(sprintf('Accuracy = %0.2f%%', accuracy*100))

%% Get highest weighted rates 



%% Save data 

% base = getPathsNT;
% save(fullfile(base, 'model_comparisons', ['pop_rate_F0_linear_' target '.mat']),...
% 	"pred_F0", "T_test", "T", "r2", "C", "accuracy", "F0s")

