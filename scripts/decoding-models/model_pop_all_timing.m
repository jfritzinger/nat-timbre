clear

%% Load in data

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Neuron_Time_All.mat'),...
	"neuron_time_all")

%% Set responses 



%% Find all rows with bassoon and oboe

[sesh_all, ~] = getTimbreSessions(nat_data);

% Get indices of interest
accuracy = [neuron_time_all.accuracy];
[~,best_ind] = sort(abs(accuracy), 'descend' );
[~,worst_ind] = sort(abs(accuracy), 'ascend' );

%% Get data

num_neurons = 50; %[1:4 5:5:20 30 40];
nmodels = length(num_neurons);
timerVal = tic;
for imodel = 1:nmodels
	timerVal = tic;

	for inrep = 1

		index = best_ind;
		num_index = 1:num_neurons(imodel);
		num_data = num_neurons(imodel);
		sesh = sesh_all(index(num_index));
		putative = {nat_data(sesh_all).putative};

		T = getAllPopTable(nat_data(sesh), 'Timing');

		% Separate out training and testing data
		[T_train, T_test] = splitData(T); 

		% Call model (classification)
		[trainedClassifier, validationAccuracy, validationPredictions] ...
			= trainClassifierPopTimeF0(T_train, F0s);

		% Plot
		figure
		C = confusionmat(T.response, validationPredictions);
		confusionchart(C)
		title([target ', Accuracy = ' num2str(round(validationAccuracy*100)) '%'])

		% Calculate accuracy
		accuracy = sum(diag(C)) / sum(C(:)); % Calculate accuracy
		title(sprintf('%d neurons, Accuracy = %0.2f%%', num_data, accuracy*100))
		accur_all(imodel, inrep) = accuracy;
		C_all{imodel, inrep} = C;

	end

% 	pop_rate_timbre.trainedClassifier = trainedClassifier;
% 	pop_rate_timbre.accuracy = accuracy;
% 	pop_rate_timbre.T = T;
% 	pop_rate_timbre.CFs = nat_data(sesh).CF;
% 	pop_rate_timbre.putative = putative;
% 	pop_rate_timbre.sesh = sesh;
% 	pop_rate_timbre.MTF = {nat_data(sesh).MTF};

	timer = toc(timerVal);
	fprintf('Models took %0.2g minutes\n', timer/60)
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)
end

%%

mean_acc = mean(accur_all,2);
std_acc = std(accur_all, [], 2);

figure
errorbar(num_neurons, mean_acc, std_acc/sqrt(2));
hold on
xlabel('Number of Neurons in Model')
ylabel('Accuracy')
legend('Best', 'Worst')

%% Save data

save(fullfile(base, 'model_comparisons', 'Pop_Time_2.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")