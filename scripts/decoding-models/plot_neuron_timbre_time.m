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
