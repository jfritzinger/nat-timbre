clear

%% Load in data

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'),...
	"neuron_time_timbre")

%% Run model 

% Find all rows with bassoon and oboe
[sesh_all, num_data] = getTimbreSessions(nat_data);

% Get indices of interest
accuracy = [neuron_time_timbre.accuracy];
[~,best_ind] = sort(abs(accuracy), 'descend' );
[~,worst_ind] = sort(abs(accuracy), 'ascend' );

num_neurons = [1:4 5:5:40 1:4 5:5:40];
nmodels = length(num_neurons);
timerVal = tic;
for imodel = 1:nmodels
	timerVal = tic;

	for inrep = 1:5

		if imodel < nmodels/2+1 % 6 good models
			index = best_ind;
		else % 6 bad models
			index = worst_ind;
		end
		num_index = 1:num_neurons(imodel);
		num_data = num_neurons(imodel);
		sesh = sesh_all(index(num_index));

		% Get table of PSTH data
		T = getTimbrePopTable(nat_data, 'Timing', sesh, num_data);

		% Call model (classification)
		[trainedClassifier, validationAccuracy, validationPredictions] = ...
			trainClassifierTimeTimbre(T);

		% Plot
		%figure
		C = confusionmat(T.Instrument, validationPredictions);
		%confusionchart(C)

		% Calculate accuracy
		accuracy = sum(diag(C)) / sum(C(:)); % Calculate accuracy
		%title(sprintf('%d neurons, Accuracy = %0.2f%%', num_data, accuracy*100))
		accur_all(imodel, inrep) = accuracy;
		C_all{imodel, inrep} = C;
	end
	timer = toc(timerVal);
	fprintf('Models took %0.2g minutes\n', timer/60)
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)
end

%% Plot results 

mean_acc = mean(accur_all,2);
std_acc = std(accur_all, [], 2);

figure
errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
hold on
errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
xlabel('Number of Neurons in Model')
ylabel('Accuracy')
legend('Best', 'Worst')

%% Save data

save(fullfile(base, 'model_comparisons', 'pop_timing_timbre.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")

