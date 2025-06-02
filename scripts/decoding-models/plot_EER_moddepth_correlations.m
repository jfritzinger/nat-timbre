%% plot_EER_moddepth_correlations
% Goal: 
%
% 

%% Load in model data  

target = 'Bassoon';
[base, ~, ~, ~] = getPathsNT();
load(fullfile(base, 'model_comparisons', ['Neuron_Time_F0_' target '.mat']), ...
	"neuron_time_F0")
accuracy(1,:) = [neuron_time_F0.accuracy]*100;


%% Load in vector strength 
load(fullfile(base,'model_comparisons', 'Data_NT_3.mat'), 'nat_data')

% Get subset of units 
VS_all = zeros(3,1);
for ii = 1:297
	VS_1 = nat_data(ii).bass_VS;
	if ~isempty(VS_1)
		VS_all(ii,:) = mean(VS_1);
	end
end
VS_all(VS_all==0) = [];


%% Load in EER/mod depth stuff 

load(fullfile(base, 'ERR_ModDpth_Oboe.mat'), "ERR", "modDepth", "pitches")
ERR_oboe = ERR;
modDepth_oboe = modDepth;
pitches_oboe = pitches;

load(fullfile(base, 'ERR_ModDpth_Bassoon.mat'), "ERR", "modDepth", "pitches")


figure
nexttile
scatter(pitches, ERR, 15, 'filled')
hold on
scatter(pitches_oboe, ERR_oboe, 15, 'filled')
title('ERR')
set(gca, 'xscale', 'log')

nexttile
scatter(pitches, modDepth, 15, 'filled')
hold on
scatter(pitches_oboe, modDepth_oboe, 15, 'filled')

mdl = fitlm(pitches, modDepth);
x = linspace(55, 590, 10);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'b')
mdl = fitlm(pitches_oboe, modDepth_oboe);
x = linspace(230, 1620, 10);
y = mdl.Coefficients{2, 1}*x + mdl.Coefficients{1, 1};
plot(x, y, 'r')

title('Modulation Depth')
set(gca, 'xscale', 'log')
legend('Bassoon', 'Oboe')

%% Correlate VS and EER/mod depth 




%% Correlate accuracy and EER/mod depth 




