%% plot_pop_timbre_timbre

%% Load in data 

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'pop_timing_timbre.mat'), ...
	"accur_all","C_all", "num_neurons", "mean_acc", "std_acc")
load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'),...
	"neuron_time_timbre")
load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')


%% Plot 

figure
tiledlayout(1, 2)

nmodels = length(num_neurons);
nexttile
errorbar(num_neurons(1:nmodels/2), mean_acc(1:nmodels/2), std_acc(1:nmodels/2)/sqrt(5), 'Color','#1b9e77', 'LineWidth',2);
hold on
errorbar(num_neurons(nmodels/2+1:end), mean_acc(nmodels/2+1:end), std_acc(nmodels/2+1:end)/sqrt(10), 'Color','k', 'LineWidth',2)
xlabel('Number of Neurons in Model')
ylabel('Accuracy')
legend('Best', 'Worst')
box off
grid on

% Plot vector strength 
ind_b = 25:40;
ind_o = [1 3:17];
accuracy = [neuron_time_timbre.accuracy];

sesh = [];
for ii = 1:length(nat_data)
	rate = nat_data(ii).bass_VS;
	rate2 = nat_data(ii).oboe_VS;
	if ~isempty(rate) && ~isempty(rate2)
		sesh = [sesh ii];
	end
end
num_data = length(sesh);

for ind = 1:num_data
	index = sesh(ind);
	VS_oboe(ind) = mean(nat_data(index).oboe_VS(ind_o));
	VS_bass(ind) = mean(nat_data(index).bass_VS(ind_b));
	VS_all(ind) = mean([VS_oboe(ind) VS_bass(ind)]);
end

% Plot VS vs accuracy
for ii = 1
	if ii == 1
		VS = VS_all;
	elseif ii == 2
		VS = VS_oboe;
	else
		VS = VS_bass;
	end
	nexttile
	scatter(VS, accuracy, 'filled',...
		'MarkerEdgeColor','k');
	hold on
	mdl = fitlm(VS, accuracy);
	x = linspace(0, 1, 20);
	y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
	plot(x, y, 'r')
	hleg = legend('Neuron', ...
		sprintf('p=%0.04f',mdl.Coefficients{2,4}));
	ylabel('Accuracy')
	xlabel('Vector Strength')
	xlim([0 1])
	ylim([0 1])
	grid on
	box off
end