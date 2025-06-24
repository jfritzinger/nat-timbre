%% plot_neuron_timbre_time

[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'),...
	"neuron_time_timbre")
neuron_time_F0 = neuron_time_timbre;


%% 

accuracy = [neuron_time_F0.accuracy]*100;
load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

% % Get subset of units
% VS_all = zeros(3,1);
% for ii = 1:297
% 	VS_1 = nat_data(ii).bass_VS;
% 	if ~isempty(VS_1)
% 		VS_all(ii,:) = mean(VS_1);
% 	end
% end
% VS_all(VS_all==0) = [];

% Get VS for all units
% Find all rows with bassoon in them
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
	for target = 1:16
		VS_1 = mean(nat_data(index).oboe_VS);
		VS_2 = mean(nat_data(index).bass_VS);
	end
	h_all = [h_bass; h_oboe];
	h = [h; h_all];
end

%% Plot VS vs accuracy 


figure
histogram(accuracy)
accuracy_all = accuracy;

nexttile
histogram(VS_all)

nexttile
scatter(VS_all, accuracy);
hold on

mdl = fitlm(VS_all, accuracy);
x = linspace(0, 0.6, 20);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')
hleg = legend('Neuron', ...
	sprintf('y = %0.2f*x+%0.2f, p=%0.04f', mdl.Coefficients{2, 1}, ...
	mdl.Coefficients{1, 1},mdl.Coefficients{2,4}));
ylabel('Accuracy (%)')
xlabel('Vector Strength')

%% Plot all histograms 


load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All_Shuffled.mat'), ...
	"neuron_time_timbre")
acc_shuffled = [neuron_time_timbre.accuracy];

load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All.mat'), ...
	"neuron_time_timbre")
acc_025ms = [neuron_time_timbre.accuracy];

load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All_Coarse100.mat'), ...
	"neuron_time_timbre")
acc_100ms = [neuron_time_timbre.accuracy];

load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All_Coarse150.mat'), ...
	"neuron_time_timbre")
acc_150ms = [neuron_time_timbre.accuracy];

load(fullfile(base, 'model_comparisons', 'Neuron_Time_Timbre_All_Coarse300.mat'), ...
	"neuron_time_timbre")
acc_300ms = [neuron_time_timbre.accuracy];

load(fullfile(base, 'model_comparisons', 'Neuron_Rate_Timbre_All.mat'), ...
	"neuron_rate_timbre")
accuracy_rate = [neuron_rate_timbre.accuracy];

%%

figure
accur_all = [acc_shuffled; acc_025ms; acc_100ms; accuracy_rate];
boxplot(accur_all')

hold on

swarmchart(ones(length(acc_shuffled), 1), acc_shuffled)
swarmchart(ones(length(acc_025ms), 1)*2, acc_025ms)
swarmchart(ones(length(acc_100ms), 1)*3, acc_100ms)
swarmchart(ones(length(accuracy_rate), 1)*4, accuracy_rate)

xlabel('Groups')
ylabel('Accuracy')
xticks(1:4)
xticklabels({'Shuffled', '0.25 ms bin', '100 ms bin', 'Rate'})

[p12, ~] = ranksum(acc_shuffled, acc_025ms);  
[p13, ~] = ranksum(acc_shuffled, acc_100ms);  
[p23, ~] = ranksum(acc_025ms, acc_100ms);  
adjusted_p = [p12, p13, p23] * 3; % Bonferroni adjustment

[p, tbl, stats] = anova1(accur_all');
[c, m, h, gnames] = multcompare(stats, 'CType', 'hsd');

