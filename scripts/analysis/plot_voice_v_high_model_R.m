% Natural Timbre 

% Load in spreadsheet
%[base, ~, ~, ~] = getPaths();
filepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
spreadsheet_name = 'model_r2_values_NT2.xlsx';
sessions = readtable(fullfile(filepath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Find sessions for target MTF type
MTF_target = 'BS';
isMTF = strcmp(sessions.MTF, MTF_target);

%% Plot scatter of R vs CF
% figure('position', [519,592,1100,240])
% tiledlayout(2, 2, 'TileSpacing','tight', 'Padding','compact', 'TileIndexing','columnmajor')
% instruments = {'Bassoon', 'Oboe'};
% for iNT = 1:2
% 	isInstrument = strcmp(sessions.Instrument, instruments{iNT});
% 	indices = find(isMTF & isInstrument);
% 
% 	nexttile
% 	hold on
% 	CFs = sessions.CF(indices);
% 	SFIE_R2 = sessions.SFIE_R(indices);
% 	energy_R2 = sessions.Energy_R(indices);
% 	%energy_R2 = sessions.Energy_R(indices);
% 	scatter(CFs, SFIE_R2, 20, 'filled', 'MarkerFaceColor','#009E73')
% 	set(gca, 'XScale', 'log')
% 
% 	nexttile
% 	hold on
% 	scatter(CFs, energy_R2, 20, 'filled', 'MarkerFaceColor','#D55E00')
% 	set(gca, 'XScale', 'log')
% end

%% Plot histogram of R values 

instruments = {'Bassoon', 'Oboe'};
figure('name',instruments{1}, 'position', [47,487,523,401])
tiledlayout(2, 2, 'TileSpacing','tight', 'Padding','compact', 'TileIndexing','columnmajor')
for ispl = 1:2

	MTF_target = 'BS';
	isMTF = strcmp(sessions.MTF, MTF_target);
	isInstrument = strcmp(sessions.Instrument, instruments{ispl});
	indices = find(isMTF & isInstrument);

	nexttile
	hold on
	SFIE_R2 = sessions.SFIE_R2(indices);
	energy_R2 = sessions.Energy_R2(indices);
	lat_inh_R2 = sessions.Lat_Inh_R2(indices);
	edges = linspace(0, 1, 20);
	histogram(SFIE_R2, edges, 'FaceAlpha',0.5)
	%title('BS, SFIE')

	%nexttile
	hold on
	histogram(energy_R2, edges,'FaceAlpha',0.5)
	histogram(lat_inh_R2, edges, 'FaceAlpha',0.5)
	%title('BS, energy')
	legend('SFIE', 'Energy', 'Lat Inh')
	title(['BS, ' instruments{ispl}])

	MTF_target = 'BE';
	isMTF = strcmp(sessions.MTF, MTF_target);
	isInstrument = strcmp(sessions.Instrument, instruments{ispl});
	indices = find(isMTF & isInstrument);

	nexttile
	hold on
	SFIE_R2 = sessions.SFIE_R2(indices);
	energy_R2 = sessions.Energy_R2(indices);
	lat_inh_R2 = sessions.Lat_Inh_R2(indices);
	edges = linspace(0, 1, 20);
	histogram(SFIE_R2, edges, 'FaceAlpha',0.5)
	%title('BE, SFIE')

	%nexttile
	hold on
	histogram(energy_R2, edges, 'FaceAlpha',0.5)
	histogram(lat_inh_R2, edges, 'FaceAlpha',0.5)
	title('BE, energy')
	legend('SFIE', 'Energy', 'Lat Inh')
	title(['BE, ' instruments{ispl}])
end

%% Split into voice and high pitch

instruments = 'Bassoon';
figure
tiledlayout(3, 1)

MTF_target = 'BE';
isMTF = strcmp(sessions.MTF, MTF_target);
isInstrument = strcmp(sessions.Instrument, instruments);
indices = find(isInstrument & isMTF);
edges = linspace(-1, 1, 21);

% SFIE
nexttile
hold on
SFIE_R2 = sessions.SFIE_voice_R(indices);
histogram(SFIE_R2, edges, 'FaceAlpha',0.5)
xline(mean(SFIE_R2), 'b', 'LineWidth',2)
SFIE_R2 = sessions.SFIE_high_R(indices);
histogram(SFIE_R2, edges, 'FaceAlpha',0.5)
xline(mean(SFIE_R2), 'r', 'LineWidth',2)
title('SFIE')
xlabel('Correlation')
ylabel('# Neurons')
legend('Voice F0', 'Mean Voice', 'High F0', 'Mean High')

% Broad inhibition
nexttile
hold on
lat_inh_R2 = sessions.Lat_Inh_voice_R(indices);
histogram(lat_inh_R2, edges, 'FaceAlpha',0.5)
xline(mean(lat_inh_R2), 'b', 'LineWidth',2)
lat_inh_R2 = sessions.Lat_Inh_high_R(indices);
histogram(lat_inh_R2, edges, 'FaceAlpha',0.5)
xline(mean(lat_inh_R2), 'r', 'LineWidth',2)

title('Broad Inhibition')
xlabel('Correlation')
ylabel('# Neurons')
legend('Voice F0', 'Mean Voice', 'High F0', 'Mean High')

% Energy
nexttile
hold on
energy_R2 = sessions.Energy_voice_R(indices);
histogram(energy_R2, edges,'FaceAlpha',0.5)
xline(mean(energy_R2), 'b', 'LineWidth',2)

energy_R2 = sessions.Energy_high_R(indices);
histogram(energy_R2, edges,'FaceAlpha',0.5)
xline(mean(energy_R2), 'r', 'LineWidth',2)

title('Energy')
xlabel('Correlation')
ylabel('# Neurons')
legend('Voice F0', 'Mean Voice', 'High F0', 'Mean High')


%% Scatter 

instruments = 'Bassoon';
figure
tiledlayout(1, 3)

isInstrument = strcmp(sessions.Instrument, instruments);
indices = find(isInstrument);
edges = linspace(0, 1, 20);

% SFIE
nexttile
hold on
SFIE_R2 = sessions.SFIE_voice_R(indices);
SFIE_R22 = sessions.SFIE_high_R(indices);
scatter(SFIE_R2, SFIE_R22, 'filled', 'MarkerEdgeColor','k')
plot([-1 1],[-1 1], 'k')
title('SFIE')
xlabel('Voice F0')
ylabel('High F0')
xline(0, 'k')
yline(0, 'k')
% xlim([0 0.7])
% ylim([0 0.7])

% Broad inhibition
nexttile
hold on
lat_inh_R2 = sessions.Lat_Inh_voice_R(indices);
lat_inh_R22 = sessions.Lat_Inh_high_R(indices);
scatter(lat_inh_R2, lat_inh_R22, 'filled', 'MarkerEdgeColor','k')
plot([-1 1],[-1 1], 'k')
xline(0, 'k')
yline(0, 'k')
title('Broad Inhibition')
xlabel('Voice F0')
ylabel('High F0')
% xlim([0 0.7])
% ylim([0 0.7])

% Energy
nexttile
hold on
energy_R2 = sessions.Energy_voice_R(indices);
energy_R22 = sessions.Energy_high_R(indices);
scatter(energy_R2, energy_R22,'filled', 'MarkerEdgeColor','k')
plot([-1 1],[-1 1], 'k')
title('Energy')
xlabel('Voice F0')
ylabel('High F0')
xline(0, 'k')
yline(0, 'k')
% xlim([0 0.7])
% ylim([0 0.7])

%% Scatter 

instruments = 'Bassoon';
figure
tiledlayout(2, 1)

isInstrument = strcmp(sessions.Instrument, instruments);
indices = find(isInstrument);
edges = linspace(-1, 1, 21);

% Voice
nexttile
hold on
SFIE_R2 = sessions.SFIE_voice_R(indices);
energy_R2 = sessions.Energy_voice_R(indices);
histogram(SFIE_R2, edges, 'FaceAlpha',0.5)
xline(mean(SFIE_R2), 'b', 'LineWidth',2)
histogram(energy_R2, edges, 'FaceAlpha',0.5)
xline(mean(energy_R2), 'r', 'LineWidth',2)

title('Voice F0')
xlabel('SFIE')
ylabel('Energy')

% High
nexttile
hold on
SFIE_R22 = sessions.SFIE_high_R(indices);
energy_R22 = sessions.Energy_high_R(indices);

histogram(SFIE_R22, edges, 'FaceAlpha',0.5)
xline(mean(SFIE_R22), 'b', 'LineWidth',2)
histogram(energy_R22, edges, 'FaceAlpha',0.5)
xline(mean(energy_R22), 'r', 'LineWidth',2)
title('High F0')
xlabel('SFIE')
ylabel('Energy')

%% Histogram

instruments = 'Bassoon';
figure
tiledlayout(2, 1)

isInstrument = strcmp(sessions.Instrument, instruments);
indices = find(isInstrument);
edges = linspace(0, 1, 20);

% Voice
nexttile
hold on
SFIE_R2 = sessions.SFIE_voice_R(indices);
energy_R2 = sessions.Energy_voice_R(indices);

scatter(SFIE_R2, energy_R2, 'filled', 'MarkerEdgeColor','k')
plot([-1 1],[-1 1], 'k')
title('Voice F0')
xlabel('SFIE')
ylabel('Energy')

% High
nexttile
hold on
SFIE_R22 = sessions.SFIE_high_R(indices);
energy_R22 = sessions.Energy_high_R(indices);

scatter(SFIE_R22, energy_R22,'filled', 'MarkerEdgeColor','k')
plot([-1 1],[-1 1], 'k')
title('High F0')
xlabel('SFIE')
ylabel('Energy')

%% ANOVA / box plots 

MTF_target = 'BS';
isMTF = strcmp(sessions.MTF, MTF_target);
isInstrument = strcmp(sessions.Instrument, instruments);
indices = find(isInstrument & isMTF);

SFIE_R2 = sessions.SFIE_voice_R(indices);
energy_R2 = sessions.Energy_voice_R(indices);
SFIE_R22 = sessions.SFIE_high_R(indices);
energy_R22 = sessions.Energy_high_R(indices);

matrix = [SFIE_R2, SFIE_R22, energy_R2, energy_R22];
[~,~,stats] = anova1(matrix);
[c,~,~,gnames] = multcompare(stats);

%% Scatter plots 
% 
% % energy vs SFIE
% figure('position', [78,337,939,616])
% tiledlayout(3, 4, 'TileIndexing','columnmajor')
% spls = [43, 63, 73, 83];
% for ispl = 1:2
% 	isInstrument = strcmp(sessions.Instrument, instruments{ispl});
% 
% 	isMTF = strcmp(sessions.MTF, 'BS');
% 	BS_ind = find(isMTF & isInstrument);
% 	SFIE_R2 = sessions.SFIE_R2(BS_ind);
% 	energy_R2 = sessions.Energy_R2(BS_ind);
% 	lat_inh_R2 = sessions.Lat_Inh_R2(BS_ind);
% 
% 	isMTF = strcmp(sessions.MTF, 'BE');
% 	BE_ind = find(isMTF & isInstrument);
% 	SFIE_R22 = sessions.SFIE_R2(BE_ind);
% 	energy_R22 = sessions.Energy_R2(BE_ind);
% 	lat_inh_R22 = sessions.Lat_Inh_R2(BE_ind);
% 
% 	nexttile
% 	hold on
% 	scatter(energy_R2, SFIE_R2, 'filled', 'MarkerEdgeColor','k')
% 	scatter(energy_R22, SFIE_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
% 	%[x, y] = findRegressionLine(energy_R2, SFIE_R2);
% 	%plot(x, y, 'b')
% 	% [x, y] = findRegressionLine(energy_R22, SFIE_R22);
% 	% plot(x, y, 'r')
% 	xlabel('Energy R^2')
% 	ylabel('SFIE R^2')
% 	title(instruments{ispl})
% 
% 	nexttile
% 	hold on
% 	scatter(energy_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k')
% 	scatter(energy_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
% 	%[x, y] = findRegressionLine(energy_R2, lat_inh_R2);
% 	%plot(x, y, 'b')
% 	% [x, y] = findRegressionLine(energy_R22, lat_inh_R22);
% 	% plot(x, y, 'r')
% 	xlabel('Energy R^2')
% 	ylabel('Lateral Inh R^2')
% 
% 	nexttile
% 	hold on
% 	scatter(SFIE_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k')
% 	scatter(SFIE_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
% 	%[x, y] = findRegressionLine(SFIE_R2, lat_inh_R2);
% 	%plot(x, y, 'b')
% 	% [x, y] = findRegressionLine(SFIE_R22, lat_inh_R22);
% 	% plot(x, y, 'r')
% 	xlabel('SFIE R^2')
% 	ylabel('Lateral Inh R^2')
% end