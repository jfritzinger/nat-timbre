clear

%% Load in data

filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/model_comparisons';
load(fullfile(filepath, 'Data_NT2.mat'), 'nat_data')

%% Get correct output of model
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

%% Get data

h_all2 = [];
for ineuron = 1:3
	if ineuron == 1
		%putative = 'R29_TT2_P3_N03'; % Predictions 44, yesterday 62?
		putative = 'R29_TT3_P2_N13';
	elseif ineuron == 2
		%putative = 'R29_TT3_P5_N05'; % Not predicting well, not sure why
		putative = 'R29_TT1_P2_N17';
	else
		%putative = 'R29_TT3_P2_N06'; % Good
		putative = 'R29_TT1_P3_N04';
	end
	index = find(strcmp({nat_data.putative}, putative));

	h_all = [];
	for itarget = 1:length(F0s)
		if strcmp(target, 'Bassoon')
			spikes_bass = nat_data(index).bass_spikerate{itarget}/1000; % ms
			spikereps_bass = nat_data(index).bass_spikerep{itarget};
		else
			spikes_bass = nat_data(index).oboe_spikerate{itarget}/1000; % ms
			spikereps_bass = nat_data(index).oboe_spikerep{itarget};
		end

		% Arrange data for SVM
		min_dis = 1;
		edges = 0:min_dis:300;
		t = 0+min_dis/2:min_dis:300-min_dis/2;
		for irep = 1:20
			h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
		end
		h_all = [h_all; h_bass];
	end

	h_all2 = [h_all2, h_all];

end

% Put data into table
response = reshape(repmat(F0s, 1, 20)', 1, []);
T = array2table(h_all2);
T.response = response';
predictors = h_all2;


%% Call model (classification)

% [trainedClassifier, validationAccuracy, validationPredictions] ...
% 	= trainClassifierPopTimeF0(T, F0s);
% 
% % Plot 
% figure
% C = confusionmat(T.response, validationPredictions);
% confusionchart(C)
% title([target ', Accuracy = ' num2str(round(validationAccuracy*100)) '%'])
% 


%% Call model (linear)

[trainedModel, validationRMSE, validationPredictions] = trainRegressionModelPopTimeF0(T);

r = corrcoef(validationPredictions, T.response);
r2 = r(1, 2)^2;

figure('Position',[31,910,910,413])
nexttile
scatter(T.response, validationPredictions, 'filled', 'MarkerFaceAlpha',0.5)
%set(gca, 'xscale', 'log', 'yscale', 'log')
% ylim([57 588])
% xlim([57 588])
hold on
plot(T.response, T.response)
title(['R^2 = ' num2str(r2)])
ylabel('Predicted log10(F0)')
xlabel('Actual log10(F0)')

logF0 = log10(F0s1);
for ii = 1:length(validationPredictions)

	differences = abs(logF0 - validationPredictions(ii));
	[~, closest_column_index] = min(differences);
	closest_cat(ii) = logF0(closest_column_index);
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



