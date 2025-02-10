%% NT_model_distributions 
clear 

% Load in spreadsheet
datapath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
spreadsheet_name = 'model_r2_values_NT.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%% Plot histogram of R values

figure('Position',[560,593,1050,255])
tiledlayout(1, 3)
fontsize = 18;

isOboe = strcmp(sessions.Instrument, 'Oboe');
isBass = strcmp(sessions.Instrument, 'Bassoon');

nexttile
hold on
edges = linspace(-1, 1, 20);
histogram(sessions.SFIE_R(isOboe), edges, 'FaceAlpha',0.5)
histogram(sessions.SFIE_R(isBass), edges, 'FaceAlpha',0.5)
title('SFIE')
set(gca, 'Fontsize', fontsize)
grid on
xlim([-1 1])
xlabel('Correlation')
ylabel('# Neurons')

nexttile
hold on
edges = linspace(-1, 1, 20);
histogram(sessions.Energy_R(isOboe), edges, 'FaceAlpha',0.5)
histogram(sessions.Energy_R(isBass), edges, 'FaceAlpha',0.5)
title('Energy')
set(gca, 'Fontsize', fontsize)
grid on
xlim([-1 1])
xlabel('Correlation')

nexttile
hold on
edges = linspace(-1, 1, 20);
histogram(sessions.Lat_Inh_R(isOboe), edges, 'FaceAlpha',0.5)
histogram(sessions.Lat_Inh_R(isBass), edges, 'FaceAlpha',0.5)
title('SFIE')
set(gca, 'Fontsize', fontsize)
legend('Oboe', 'Bassoon')
grid on
xlim([-1 1])
xlabel('Correlation')
%%

figure('Position',[560,593,1050,255])
tiledlayout(1, 3)
fontsize = 18;

isOboe = strcmp(sessions.Instrument, 'Oboe');
isBass = strcmp(sessions.Instrument, 'Bassoon');

nexttile
hold on
edges = linspace(0, 1, 20);
histogram(sessions.SFIE_R2(isOboe), edges, 'FaceAlpha',0.5)
histogram(sessions.SFIE_R2(isBass), edges, 'FaceAlpha',0.5)
title('SFIE')
set(gca, 'Fontsize', fontsize)
grid on
xlim([0 1])
xlabel('Variance Explained (R^2)')
ylabel('# Neurons')

nexttile
hold on
edges = linspace(0, 1, 20);
histogram(sessions.Energy_R2(isOboe), edges, 'FaceAlpha',0.5)
histogram(sessions.Energy_R2(isBass), edges, 'FaceAlpha',0.5)
title('Energy')
set(gca, 'Fontsize', fontsize)
grid on
xlim([0 1])
xlabel('Variance Explained (R^2)')

nexttile
hold on
edges = linspace(0, 1, 20);
histogram(sessions.Lat_Inh_R2(isOboe), edges, 'FaceAlpha',0.5)
histogram(sessions.Lat_Inh_R2(isBass), edges, 'FaceAlpha',0.5)
title('SFIE')
set(gca, 'Fontsize', fontsize)
legend('Oboe', 'Bassoon')
grid on
xlim([0 1])
xlabel('Variance Explained (R^2)')
