clear

%% Load in data

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
load(fullfile(base, 'model_comparisons', 'Neuron_Time_All.mat'),...
	"neuron_time_all")

%% Set responses 

% Get bassoon stimulus
F0s_b = getF0s('Bassoon');
F0s_o = getF0s('Oboe');

% Get into 'Response' 
response_b = cell(75,1);
for ii = 1:length(F0s_b)
	response_b{ii} = ['B_' num2str(F0s_b(ii))];
end
for ii = 1:length(F0s_o)
	response_b{ii+length(F0s_b)} = ['O_' num2str(F0s_o(ii))];
end
response = reshape(repmat(response_b, 1, 20)', 1, []);
response = response';

%% Find all rows with bassoon and oboe

[sesh_all, ~] = getTimbreSessions(nat_data);

% Get indices of interest
accuracy = [neuron_time_all.accuracy];
[~,best_ind] = sort(abs(accuracy), 'descend' );
[~,worst_ind] = sort(abs(accuracy), 'ascend' );



%% Get data

num_neurons = [50 75 100]; % [1:4 5:5:20 30 40];
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

		h_all2 = [];
		h_all = [];
		for ineuron = 1:num_data
			h_all_bass = [];
			for itarget = 1:length(F0s_b)

				% Arrange bassoon data for SVM
				spikes = nat_data(sesh(ineuron)).bass_spikerate{itarget}/1000; % ms
				spikereps = nat_data(sesh(ineuron)).bass_spikerep{itarget};
				min_dis = 0.25;
				edges = 0:min_dis:300;
				t = 0+min_dis/2:min_dis:300-min_dis/2;
				for irep = 1:20
					h_bass(irep, :) = histcounts(spikes(spikereps==irep), edges);
				end
				h_all_bass = [h_all_bass; h_bass];
			end

			h_all_oboe = [];
			for itarget = 1:length(F0s_o)


				spikes = nat_data(sesh(ineuron)).oboe_spikerate{itarget}/1000; % ms
				spikereps = nat_data(sesh(ineuron)).oboe_spikerep{itarget};
				for irep = 1:20
					h_oboe(irep, :) = histcounts(spikes(spikereps==irep), edges);
				end
				h_all_oboe = [h_all_oboe; h_oboe];
			end
			h_all = [h_all_bass; h_all_oboe];
			h_all2 = [h_all2, h_all];
		end
		T = array2table(h_all2);
		T.response = response;
		predictors = h_all2;

		% Call model (classification)
		[trainedClassifier, validationAccuracy, validationPredictions] ...
			= trainClassifierPopTimeF0(T, F0s);

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