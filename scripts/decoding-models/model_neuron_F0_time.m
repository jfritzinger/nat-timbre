%% model_neuron_rate_F0
clear

%% Load in data

[base, datapath, savepath, ppi] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
% load(fullfile(base, 'model_comparisons',  'Model_NT2.mat'), 'nat_model')
% nat_data = nat_model;

%% Get correct output of model 
target = 'Bassoon';
%target = 'Oboe';
%target = 'Invariant';

% Get bassoon stimulus
F0s = getF0s(target);
F0s = round(F0s);

%% Get all rates for each repetition for bassoon (one example neuron)

[sesh, num_data] = getF0Sessions(nat_data, target);
min_dis = [0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 10 15 20 30 40 50 100 150 300];
accuracy = zeros(length(min_dis), num_data);
for ind = 1:num_data
	timerVal = tic;
	parfor idist = 1:length(min_dis)

		min_dis2 = min_dis(idist);
		index = sesh(ind);
		T = getF0NeuronTable(nat_data, target, index, F0s, 'Data', min_dis2);

		for imodelrep = 1 %:5

			% Split into training and testing
			[T_train, T_test] = splitData_Reps(T);

			% Train model
			%[validationPredictions] = trainClassifierNeuronTimeF0(T_train, F0s);
			[trainedClassifier, validationPredictions] = trainClassifierNeuronTimeF0(T_train, F0s);

			% To make predictions with the returned 'trainedClassifier' on new data
			[yfit,scores] = trainedClassifier.predictFcn(T_test);
			C = confusionmat(T_test.Response, yfit);
			accuracy(idist, ind) = sum(diag(C))/sum(C, 'all');
		end

		% 	figure
		% 	confusionchart(C)
		% 	title(num2str(validationAccuracy*100))

		% Save data for each
		neuron_time_F0(idist, ind).putative = nat_data(index).putative;
		neuron_time_F0(idist, ind).CF = nat_data(index).CF;
		neuron_time_F0(idist, ind).MTF = nat_data(index).MTF;
		neuron_time_F0(idist, ind).response = T.Response;
		neuron_time_F0(idist, ind).T = T;
		neuron_time_F0(idist, ind).validationPredictions = validationPredictions;
		neuron_time_F0(idist, ind).accuracy = accuracy(idist, ind);
		neuron_time_F0(idist, ind).C = C;

		% fprintf('%d/%d, %d/%d, %0.2f%% done! Acc = %0.2g%%\n', idist, ...
		% 	length(min_dis), ind, num_data, ...
		% 	ind/num_data*100, accuracy(idist, ind)*100)
	end
	time = toc(timerVal)/60;
	fprintf('%d/%d, %0.2f%% done! Time=%0.2f minutes\n',...
		ind, num_data, ind/num_data*100,time)
end

%% Save struct of data 

save(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '_All.mat']), ...
	"neuron_time_F0", '-v7.3')

%% Plot results 

figure
tiledlayout(7, 3, 'TileSpacing','tight', 'Padding','compact')
for i = 1:length(min_dis)
	nexttile
	histogram(accuracy(i,:)*100,51)
	mean_F0 = mean(accuracy(i,:));
	hold on
	xline(50, 'k')
	xline(mean_F0*100, 'r', 'LineWidth',1)
	ylabel('# Neurons')
	xlabel('Prediction Accuracy (%)')
	xlim([0 100])
	mean_all = mean(accuracy(i,:), 'all');
	fprintf('Mean for all = %0.4f\n', mean_all)
	title(['Min Dist = ' num2str(min_dis(i)) ', Mean = ' num2str(round(mean_all*100))])
end

%% Statistics 

figure
boxplot(accuracy')
xticks(1:21)
xticklabels(min_dis)
anova(accuracy')
[p,tbl,stats] = anova1(accuracy');
results = multcompare(stats);
yticklabels(min_dis)
ylabel('PSTH bin width (ms)')
yticklabels(fliplr(min_dis))
grid on