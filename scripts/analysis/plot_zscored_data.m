%% 


clear

%% Load in data 

load('Data_NT_Matrix.mat', 'data_all')

%% Plot matrix in different ways 

% Split into MTF type 
isMTF = strcmp(data_all.MTF_shapes, 'BS');
rates_z = data_all.rates_z(isMTF,:);
CFs = data_all.CFs(isMTF);


% Plot
figure
scatter(CFs, rates_z(:,4))