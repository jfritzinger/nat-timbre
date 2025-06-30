clear

%% Load in data
%target = 'Bassoon';
%target = 'Oboe';
target = 'Invariant';

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
if strcmp(target, 'Invariant')
	load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Bassoon.mat'),...
		"neuron_time_F0")
else
	load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']),...
		"neuron_time_F0")
end

%% Get correct output of model

% Get stimulus
F0s = getF0s(target);
[sesh_all, num_data_all] = getF0Sessions(nat_data, target);

% Get indices of interest
accuracy = [neuron_time_F0.accuracy];
[~,best_ind] = sort(abs(accuracy), 'descend');
[~,worst_ind] = sort(abs(accuracy), 'ascend');

%poolobj = parpool();
% figure
% bar(accuracy(originalpos))

%% Get data

num_neurons = [1:4 5:5:40 1:4 5:5:40];
nmodels = length(num_neurons);
timerVal = tic;
for imodel = 1:nmodels
	timerVal = tic;

	for inrep = 1

		if imodel < nmodels/2+1 % 6 good models
			index = best_ind;
		else % 6 bad models
			index = worst_ind;
		end
		num_index = 1:num_neurons(imodel);
		num_data = num_neurons(imodel);
		sesh = sesh_all(index(num_index));

		T = getF0PopTable(nat_data, target, sesh, F0s, num_data, [], 'Timing');

		% Call model (classification)
		[trainedClassifier, validationAccuracy, validationPredictions] ...
			= trainClassifierPopTimeF0(T, F0s);

		% Plot
		%figure
		C = confusionmat(T.response, validationPredictions);
		%confusionchart(C)
		%title([target ', Accuracy = ' num2str(round(validationAccuracy*100)) '%'])

		% Calculate accuracy
		accuracy = sum(diag(C)) / sum(C(:)); % Calculate accuracy
		accur_all(imodel, inrep) = accuracy;
		C_all{imodel, inrep} = C;
	end
	timer = toc(timerVal);
	fprintf('Models took %0.2g minutes\n', timer/60)
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)
end

%% Plot

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

save(fullfile(base, 'model_comparisons', ['Pop_Timing_F0_' target '.mat']), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")

