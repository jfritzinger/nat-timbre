clear

%% Load in data
target = 'Bassoon';
%target = 'Oboe';
%target = 'Invariant';

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
if strcmp(target, 'Invariant')
	load(fullfile(base, 'model_comparisons', 'Neuron_Time_F0_Bassoon.mat'),...
		"neuron_time_F0")
else
	load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']),...
		"neuron_time_F0")
end

% load(fullfile(base, 'model_comparisons',  'Model_NT2.mat'), 'nat_model')
% nat_data = nat_model;
% load(fullfile(base, 'model_comparisons', ['Model_Neuron_Time_F0_' target '.mat']),...
% 	"neuron_time_F0")
% %neuron_time_F0 = neuron_time_F0(1:183);
% neuron_time_F0 = neuron_time_F0(1:174);


%% Get correct output of model

% Get stimulus
F0s = getF0s(target);
[sesh_all, num_data_all] = getF0Sessions(nat_data, target);

% Get indices of interest
accuracy = [neuron_time_F0.accuracy];
[~,best_ind] = sort(abs(accuracy), 'descend');
[~,worst_ind] = sort(abs(accuracy), 'ascend');

%% Get data

%num_neurons = [1:4 5:5:40 40:10:100 1:4 5:5:40 40:10:100];
num_neurons = [1:4 5:5:30 1:4 5:5:30];

nmodels = length(num_neurons);
timerVal = tic;
for imodel = 1:nmodels
	timerVal = tic;

	for inrep = 1:2

		if imodel < nmodels/2+1 % 6 good models
			index = best_ind;
		else % 6 bad models
			index = worst_ind;
		end
		num_index = 1:num_neurons(imodel);
		num_data = num_neurons(imodel);
		sesh = sesh_all(index(num_index));

		% Get data into table
		T = getF0PopTable(nat_data, target, sesh, F0s, num_data, [], 'Timing', 'Data');

		% Separate out training and testing data
		[T_train, T_test] = splitData(T); 

		% Call model (classification)
		[trainedClassifier, validationAccuracy, validationPredictions] ...
			= trainClassifierPopTimeF0(T_train, F0s);

		% To make predictions with the returned 'trainedClassifier' on new data
		[yfit,scores] = trainedClassifier.predictFcn(T_test);
		C = confusionmat(T_test.response, yfit);
		accuracy(imodel, inrep) = sum(diag(C))/sum(C, 'all');

	end
	timer = toc(timerVal);
	fprintf('Models took %0.2g minutes, Accuracy = %0.1f%%\n', timer/60, accuracy*100)
	fprintf('%d/%d, %0.2f%% done!\n', imodel, nmodels, imodel/nmodels*100)
end

%% Plot

figure
mean_acc = mean(accuracy,2);
std_acc = std(accuracy, [], 2);
nexttile
errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(10));
hold on
errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10))
xlabel('Number of Neurons in Model')
ylabel('Accuracy')
legend('Best', 'Worst')

%% Save data

% save(fullfile(base, 'model_comparisons', ['Pop_Timing_F0_' target '.mat']), ...
% 	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")
% save(fullfile(base, 'model_comparisons', ['Model_Pop_Timing_F0_' target '.mat']), ...
%  	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")

