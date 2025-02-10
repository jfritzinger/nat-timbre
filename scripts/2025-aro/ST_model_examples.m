%% ST_model_examples

%% aim2_model_fit_examples

% BS EXAMPLES
% 1. Sometimes all predictions are bad
% putative = 'R29_TT1_P3_N05';
% CF = 4000;
% putative = 'R25_TT4_P7_N01';
% CF = 1516;


%% Load in examples and plot
fontsize1 = 20;
fontsize = 16;

figure('Position',[1720,621,1000,500])
tiledlayout(2, 3, 'TileIndexing','columnmajor', 'Padding','compact');
for ii = 1:6
	switch ii
		case 1 % 1. Lat inh > energy % BE EXAMPLES
			putative = 'R24_TT2_P12_N10';
			CF = 3983;
			ispl = 2;
		case 2 % 1. Energy > lat inh % BE EXAMPLES
			putative = 'R27_TT2_P8_N05';
			CF = 5278;
			ispl = 2;
		case 3
			% putative = 'R27_TT2_P8_N02';
			% CF = 2000;
			putative = 'R25_TT2_P8_N11';
			CF = 1741;
			ispl = 2;
		case 4
			putative = 'R24_TT2_P12_N05';
			CF = 1740;
			ispl = 2;
		case 5 % 3. Complexities off CF aren't explained by any model
			putative = 'R29_TT3_P5_N07';
			CF = 1320;
			ispl = 2;

		case 6 % 3. Complexities off CF aren't explained by any model
			putative = 'R29_TT3_P2_N04';
			CF = 6063;
			ispl = 1;
			% putative = 'R25_TT1_P8_N15';
			% CF = 2033;
			% ispl = 2;
	end

	% Load in data
	[base, datapath, savepath, ppi] = getPaths();
	load(fullfile(datapath, 'neural_data', [putative '.mat']))

	% Load in model data
	modelpath = '/Volumes/Synth-Timbre/data/manuscript';
	load(fullfile(modelpath,'SFIE_model', [putative '_SFIE.mat']), 'SFIE')
	load(fullfile(modelpath,'energy_model', [putative '_Energy.mat']), 'energy')
	load(fullfile(modelpath,'lat_inh_model', [putative '_Lat_Inh.mat']), 'lat_inh')

	% Plot synthetic timbre (raw)
	spls = [43, 63, 73, 83];
	
	nexttile
	param_ST = data(5+ispl, 2);

	if isempty(param_ST{1})
		continue
	end

	data_ST = analyzeST(param_ST, CF);
	data_ST = data_ST{1};
	linewidth = 2;

	% Z-score
	rate = zscore(data_ST.rate);
	rate_sm = zscore(data_ST.rates_sm);
	hold on
	plot(data_ST.fpeaks,rate, 'linewidth', 0.9, 'Color',"#0072BD");
	errorbar(data_ST.fpeaks,rate, 1/sqrt(30), 'linewidth', 0.9, 'Color','k');
	plot(data_ST.fpeaks,rate_sm, 'linewidth', linewidth,'Color','k');
	ylim([-4 4])

	% Normalize and plot models
	plot(data_ST.fpeaks, zscore(energy{ispl}.rate), 'LineWidth',linewidth, 'Color','#4634F1')
	plot(data_ST.fpeaks, zscore(SFIE{ispl}.rate), 'LineWidth',linewidth, 'Color','#009E73')
	plot(data_ST.fpeaks, zscore(lat_inh{ispl}.rate), 'LineWidth',linewidth, 'Color','#D55E00')

	% Annotate SFIE model R^2
	message = sprintf('R^2 SFIE = %.02f', SFIE{ispl}.R2);
	text(0.05, 0.22, message, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',fontsize)

	% Annotate energy model R^2
	message = sprintf('R^2 Energy = %.02f', energy{ispl}.R2);
	text(0.05, 0.12, message, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',fontsize)

	% Annotate lateral inhibition model R^2
	message = sprintf('R^2 Broad inh = %.02f', lat_inh{ispl}.R2);
	text(0.05, 0.32, message, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',fontsize)

	% Plot parameters 
	plot_range = [param_ST{1}.fpeaks(1) param_ST{1}.fpeaks(end)];
	xline(CF, '--', 'Color',[0.7 0.7 0.7], 'linewidth', linewidth)
	xlim(plot_range)
	if ii == 2 || ii == 4 || ii == 6
		xlabel('Spec. Peak Freq. (Hz)')
	end
	if ii == 1 || ii == 2
		ylabel('Z-score')
	end
	box on
	if ii == 1
		title('BE Examples')
	elseif ii == 3
		title('BS Examples')
	elseif ii == 5
		title('Complex Examples')
	end
	grid on
	set(gca, 'fontsize', fontsize1)
	% legend('', '', '','', '', '', '', '', '', '', ...
	% 	'', '', '','', '', '', '', '', '', '', ...
	% 	'', '', '','', '', '', '', '', '', '', ...
	% 	'', '', '','', '', '', '', '', '', '', '', '', 'Data', ...
	% 	'Energy', 'On-CF', 'Broad Inh.', 'Location','eastoutside')

end
