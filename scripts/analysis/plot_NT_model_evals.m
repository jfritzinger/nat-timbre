% Natural Timbre 

% Load in spreadsheet
%[base, ~, ~, ~] = getPaths();
spreadsheet_name = 'model_r2_values_NT2.xlsx';
sessions = readtable(spreadsheet_name, 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Find sessions for target MTF type
MTF_target = 'BS';
isMTF = strcmp(sessions.MTF, MTF_target);

% Plot scatter of R vs CF
figure('position', [519,592,1100,240])
tiledlayout(2, 2, 'TileSpacing','tight', 'Padding','compact', 'TileIndexing','columnmajor')
instruments = {'Bassoon', 'Oboe'};
for iNT = 1:2
	isInstrument = strcmp(sessions.Instrument, instruments{iNT});
	indices = find(isMTF & isInstrument);

	nexttile
	hold on
	CFs = sessions.CF(indices);
	SFIE_R2 = sessions.SFIE_R(indices);
	energy_R2 = sessions.Energy_R(indices);
	%energy_R2 = sessions.Energy_R(indices);
	scatter(CFs, SFIE_R2, 20, 'filled', 'MarkerFaceColor','#009E73')
	set(gca, 'XScale', 'log')

	nexttile
	hold on
	scatter(CFs, energy_R2, 20, 'filled', 'MarkerFaceColor','#D55E00')
	set(gca, 'XScale', 'log')
end

%% Plot histogram of R values 
for ispl = 1:2
	figure('name',instruments{ispl}, 'position', [47,487,523,401])
	tiledlayout(1, 2, 'TileSpacing','tight', 'Padding','compact', 'TileIndexing','columnmajor')


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
	title('BS')

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
	%title('BE, energy')
	legend('SFIE', 'Energy', 'Lat Inh')
	title('BE')
end

%% Scatter plots 

% energy vs SFIE
figure('position', [78,337,939,616])
tiledlayout(3, 4, 'TileIndexing','columnmajor')
spls = [43, 63, 73, 83];
for ispl = 1:2
	isInstrument = strcmp(sessions.Instrument, instruments{ispl});

	isMTF = strcmp(sessions.MTF, 'BS');
	BS_ind = find(isMTF & isInstrument);
	SFIE_R2 = sessions.SFIE_R2(BS_ind);
	energy_R2 = sessions.Energy_R2(BS_ind);
	lat_inh_R2 = sessions.Lat_Inh_R2(BS_ind);

	isMTF = strcmp(sessions.MTF, 'BE');
	BE_ind = find(isMTF & isInstrument);
	SFIE_R22 = sessions.SFIE_R2(BE_ind);
	energy_R22 = sessions.Energy_R2(BE_ind);
	lat_inh_R22 = sessions.Lat_Inh_R2(BE_ind);

	nexttile
	hold on
	scatter(energy_R2, SFIE_R2, 'filled', 'MarkerEdgeColor','k')
	scatter(energy_R22, SFIE_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
	%[x, y] = findRegressionLine(energy_R2, SFIE_R2);
	%plot(x, y, 'b')
	% [x, y] = findRegressionLine(energy_R22, SFIE_R22);
	% plot(x, y, 'r')
	xlabel('Energy R^2')
	ylabel('SFIE R^2')
	title(instruments{ispl})

	nexttile
	hold on
	scatter(energy_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k')
	scatter(energy_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
	%[x, y] = findRegressionLine(energy_R2, lat_inh_R2);
	%plot(x, y, 'b')
	% [x, y] = findRegressionLine(energy_R22, lat_inh_R22);
	% plot(x, y, 'r')
	xlabel('Energy R^2')
	ylabel('Lateral Inh R^2')

	nexttile
	hold on
	scatter(SFIE_R2, lat_inh_R2, 'filled', 'MarkerEdgeColor','k')
	scatter(SFIE_R22, lat_inh_R22, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',"#D95319")
	%[x, y] = findRegressionLine(SFIE_R2, lat_inh_R2);
	%plot(x, y, 'b')
	% [x, y] = findRegressionLine(SFIE_R22, lat_inh_R22);
	% plot(x, y, 'r')
	xlabel('SFIE R^2')
	ylabel('Lateral Inh R^2')
end