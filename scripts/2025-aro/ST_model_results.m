%% ST_model_results

%% Synthetic Timbre
clear 

% Load in spreadsheet
[base, datapath, ~, ~] = getPaths();
spreadsheet_name = 'model_r2_values_ST.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

%% Simplified scatter plots 
fontsize = 20;

figure('position', [1796,680,600,275])
tiledlayout(1, 2)
spls = [43, 63, 73, 83];
ispl = 2;
isSPL = sessions.SPL == spls(ispl);

% not sig
for ii = 1:2

	isBS = strcmp(sessions.MTF, 'BS');
	isBE = strcmp(sessions.MTF, 'BE');
	if ii == 1
		isnotsig = sessions.p_s_e>0.05;
		x_R2 = sessions.Energy_R(isBS & ~isnotsig);
		x_R22 = sessions.Energy_R(isBE & ~isnotsig);
		x_non = sessions.Energy_R(isnotsig);
		
		y_R2 = sessions.SFIE_R(isBS & ~isnotsig);
		y_R22 = sessions.SFIE_R(isBE & ~isnotsig);
		y_non = sessions.SFIE_R(isnotsig);
	elseif ii == 2

		isnotsig = sessions.p_l_e>0.05;
		x_R2 = sessions.Energy_R(isBS & ~isnotsig);
		x_R22 = sessions.Energy_R(isBE & ~isnotsig);
		x_non = sessions.Energy_R(isnotsig);

		y_R2 = sessions.Lat_Inh_R(isBS & ~isnotsig);
		y_R22 = sessions.Lat_Inh_R(isBE & ~isnotsig);
		y_non = sessions.Lat_Inh_R(isnotsig);
	else

		isnotsig = sessions.p_l_s>0.05;
		x_R2 = sessions.SFIE_R(isBS & ~isnotsig);
		x_R22 = sessions.SFIE_R(isBE & ~isnotsig);
		x_non = sessions.SFIE_R(isnotsig);

		y_R2 = sessions.Lat_Inh_R(isBS & ~isnotsig);
		y_R22 = sessions.Lat_Inh_R(isBE & ~isnotsig);
		y_non = sessions.Lat_Inh_R(isnotsig);

	end

	nexttile
	hold on
	scatter(x_R2, y_R2, 'filled', 'MarkerEdgeColor','k', "MarkerFaceAlpha",0.7)
	scatter(x_R22, y_R22, 'filled', 'MarkerEdgeColor','k', ...
		'MarkerFaceColor',"#D95319", "MarkerFaceAlpha",0.7)
	%scatter(x_non, y_non, 'MarkerEdgeColor','k')
	plot([-1 1], [-1 1], 'k')
	xlabel('Energy R')
	ylabel('SFIE R')
	if ii == 1
		title('SFIE / Energy')
	elseif ii == 2
		title('Broad Inh / Energy')
	else
		title('Lat Inh / SFIE')
	end
	set(gca, 'fontsize', fontsize)
	yline(0)
	xline(0)
	if ii == 2
		legend('BS', 'BE', '', '', 'location', 'south')
	end
end

% %% Histograms 
% fontsize = 20;
% %Plot histograms 
% figure('position', [47,324,523,564])
% tiledlayout(3, 2, 'TileSpacing','tight', 'Padding','compact', 'TileIndexing','columnmajor')
% spls = [43, 63, 73, 83];
% for ispl = 2
% 
% 	MTF_target = 'BS';
% 	isMTF = strcmp(sessions.MTF, MTF_target);
% 	isSPL = sessions.SPL == spls(ispl);
% 	indices = find(isMTF & isSPL);
% 
% 	nexttile
% 	hold on
% 	CFs = sessions.CF(indices);
% 	SFIE_R2 = sessions.SFIE_R(indices);
% 	energy_R2 = sessions.Energy_R(indices);
% 	lat_inh_R2 = sessions.Lat_Inh_R(indices);
% 	edges = linspace(-1, 1, 20);
% 	histogram(SFIE_R2, edges, 'FaceColor','#009E73', 'FaceAlpha',0.5)
% 	title(['BS SFIE, ' num2str(spls(ispl)) ' dB SPL'])
% 	set(gca, 'fontsize', fontsize)
% 
% 	nexttile
% 	histogram(energy_R2, edges, 'FaceColor','#D55E00', 'FaceAlpha',0.5)
% 	title('BS, energy')
% 	set(gca, 'fontsize', fontsize)
% 
% 	nexttile
% 	histogram(lat_inh_R2, edges, 'FaceAlpha',0.5)
% 	title('BS, Lateral Inhibition')
% 	set(gca, 'fontsize', fontsize)
% 
% 	MTF_target = 'BE';
% 	isMTF = strcmp(sessions.MTF, MTF_target);
% 	isSPL = sessions.SPL == spls(ispl);
% 	indices = find(isMTF & isSPL);
% 
% 	nexttile
% 	hold on
% 	CFs = sessions.CF(indices);
% 	SFIE_R2 = sessions.SFIE_R(indices);
% 	energy_R2 = sessions.Energy_R(indices);
% 	lat_inh_R2 = sessions.Lat_Inh_R(indices);
% 	histogram(SFIE_R2, edges, 'FaceColor','#009E73', 'FaceAlpha',0.5)
% 	title(['BE SFIE, ' num2str(spls(ispl)) ' dB SPL'])
% 	set(gca, 'fontsize', fontsize)
% 
% 	nexttile
% 	histogram(energy_R2, edges, 'FaceColor','#D55E00', 'FaceAlpha',0.5)
% 	title('BE, energy')
% 	set(gca, 'fontsize', fontsize)
% 
% 	nexttile
% 	histogram(lat_inh_R2, edges, 'FaceAlpha',0.5)
% 	title('BE, Lateral Inhibition')
% 	set(gca, 'fontsize', fontsize)
% end