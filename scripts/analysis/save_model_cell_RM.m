%% save_model_predictions_NT
%
% Script that runs either SFIE single-cell, energy, or SFIE population
% models for each neuron with a response to synthetic timbre
%
% Author: J. Fritzinger
% Created: 2022-09-15; Last revision: 2024-09-25
%
% -------------------------------------------------------------------------
clear

%%

model_type = 'SFIE';
%model_type = 'Energy';
%model_type = 'Lat_Inh';

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPathsNT();
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(base, spreadsheet_name), 'PreserveVariableNames',true);
num_data = size(sessions, 1);

% Find sessions for target synthetic timbre response
timbre(:,1) = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
timbre(:,2) = cellfun(@(s) contains(s, 'R'), sessions.Oboe);

%%

has_data = timbre(:,1) | timbre(:,2);
indices = find(has_data);
num_index = length(indices);
for isesh = 4:num_index
	timerVal = tic;

	% Load in session
	putative = sessions.Putative_Units{indices(isesh)};
	CF = sessions.CF(indices(isesh));
	load(fullfile(datapath, [putative '.mat']))
	MTF_shape = sessions.MTF{indices(isesh)};
	params_RM = data{2, 2};

	if strcmp(MTF_shape, 'BS')
		BMF = sessions.WMF(indices(isesh));
	elseif strcmp (MTF_shape, 'BE')
		BMF = sessions.BMF(indices(isesh));
	else
		BMF = 100;
	end

	if strcmp(model_type, 'SFIE')
		%AN = cell(2, 1);
		%SFIE = cell(2,1);
	elseif strcmp(model_type, 'Energy')
		Fs = 100000;
		gamma_param.srate = Fs;
		energy = cell(2, 1);
		gamma_param.fc = CF;
	elseif strcmp(model_type, 'Lat_Inh')
		AN_lat_inh = cell(2,1);
		lat_inh = cell(2,1);
	end


	% Analysis
	if isempty(params_RM)
		% Did not record natural timbre response for this instrument
	else
		data_RM = analyzeRM(params_RM);

		if isempty(data_RM.rates) | sum(data_RM.rates) == 0

		else

			% Set up stimuli
			params_RM.Fs = 100000;
			params_RM.mnrep = 5;
			params_RM.dur = params_RM.dur/1000;
			[params_RM] = generate_RM(params_RM);
			params_RM.num_stim = size(params_RM.stim, 1);

			if strcmp(model_type, 'SFIE')

				% Model parameters
				model_params.type = 'SFIE';
				model_params.range = 2; % 1 = population model, 2 = single cell model
				model_params.species = 1; % 1 = cat, 2 = human
				model_params.BMF = BMF;
				model_params.CF_range = CF;
				model_params.num_CFs = 1;
				model_params.CFs = CF;
				model_params.nAN_fibers_per_CF = 5;
				model_params.cohc = 1; % (0-1 where 1 is normal)
				model_params.cihc = 1; % (0-1 where 1 is normal)
				model_params.nrep = 1; % how many times to run the AN model
				model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
				model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
				model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
				model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
				model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)

				% Run model
				AN_temp = modelAN(params_RM, model_params); % HSR for IC input
				SFIE_temp = wrapperIC(AN_temp.an_sout, params_RM, model_params); % SFIE output

				% Plot
				if strcmp(MTF_shape, 'BS')
					[~, rate, rate_std] = plotRM(params_RM, SFIE_temp.average_ic_sout_BS, 0);
					SFIE.rate = rate;
					SFIE.rate_std = rate_std;
					R = corrcoef(data_RM.rates,rate); % Correlation
					SFIE.R = R(1, 2);
					SFIE.R2 = R(1, 2).^2;
					SFIE.SFIE_temp = SFIE_temp;

				elseif strcmp(MTF_shape, 'BE')
					[~, rate, rate_std] = plotRM(params_RM, SFIE_temp.average_ic_sout_BE, 0);
					SFIE.rate = rate;
					SFIE.rate_std = rate_std;
					R = corrcoef(data_RM.rates,rate); % Correlation
					SFIE.R = R(1, 2);
					SFIE.R2 = R(1, 2).^2;
					SFIE.SFIE_temp = SFIE_temp;
				else
					SFIE.rate = [];
					SFIE.rate_std = [];
					SFIE.R = [];
					SFIE.R2 = [];
					SFIE.PSTH = [];
					SFIE.SFIE_temp = [];
				end

				[~, rate, rate_std] = plotRM(params_RM, AN_temp.average_AN_sout, 0);
				AN.rate = rate;
				AN.rate_std = rate_std;
				AN.AN_temp = AN_temp;

			elseif strcmp(model_type, 'Energy') % Energy model

				stimulus = [params_RM.stim zeros(size(params_RM.stim,1),0.1*Fs)];
				tvals = (1:length(stimulus))/Fs;
				gamma_IF_reg = zeros(1,length(tvals));
				impaired = 0; % 0 = not impaired; 1 = 'impaired'
				pin_gamma = zeros(size(stimulus, 1), Fs*params_RM.dur+0.1*Fs);
				for istim = 1:size(stimulus, 1)
					pin_gamma(istim,:) = gamma_filt(stimulus(istim,:),gamma_param,impaired);
				end
				pin_gamma = pin_gamma(:,1:params_RM.dur*Fs);
				energ_out = sqrt(mean(pin_gamma.^2,2));
				[rate, rate_std, pitch] = plotNT(params_RM, energ_out, 0);
				R_int = corrcoef(data_RM.rate,rate);
				energy.energ_out = energ_out;
				energy.rate = rate;
				energy.rate_std = rate_std;
				energy.target = params_RM.target;
				energy.pitch = pitch;
				energy.R = R_int(1,2);
				energy.R2 =  R_int(1, 2).^2;

			elseif strcmp(model_type, 'Lat_Inh')
				if strcmp(MTF_shape, 'BS')
					S = 0.25; % Strength, S =
					D = 0; % Delay, D =
					oct_range = 0.75; % CF range =
				else
					S = 0.4; % Strength, S =
					D = 0; % Delay, D =
					oct_range = 0.5; % CF range =
				end

				% Model parameters
				model_params.type = 'Lateral Model';
				model_params.range = 2; % 1 = population model, 2 = single cell model
				model_params.species = 1; % 1 = cat, 2 = human
				model_params.BMF = 100;
				model_params.num_CFs = 1;
				model_params.nAN_fibers_per_CF = 10;
				model_params.cohc = 1; % (0-1 where 1 is normal)
				model_params.cihc = 1; % (0-1 where 1 is normal)
				model_params.nrep = 1; % how many times to run the AN model
				model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
				model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
				model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
				model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
				model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
				model_params.lateral_CF = [CF*2^(-1*oct_range), CF, CF*2^oct_range];
				model_params.CFs = model_params.lateral_CF;
				model_params.CF_range = model_params.CFs(2);

				% Lateral model parameters
				model_params.config_type = 'BS inhibited by off-CF BS';
				lm_params = [S S D];

				% Run model
				AN_temp = modelLateralAN(params_RM, model_params);
				latinh_temp = modelLateralSFIE(params_RM, model_params,...
					AN_temp.an_sout, AN_temp.an_sout_lo, AN_temp.an_sout_hi,...
					'CS_params', lm_params);

				if strcmp(MTF_shape, 'BS') || strcmp(MTF_shape, 'BE')
					[rate, rate_std, pitch] = plotMTF(params_RM, latinh_temp.avIC, 0);
					lat_inh.rate = rate;
					lat_inh.rate_std = rate_std;
					lat_inh.target = params_RM.target;
					lat_inh.pitch = pitch;
					R = corrcoef(data_RM.rate,rate); % Correlation
					lat_inh.R = R(1, 2);
					lat_inh.R2 = R(1, 2).^2;
					lat_inh.PSTH = plotNT_PSTH(params_RM, latinh_temp.ic, 0);
				else
					lat_inh.rate = [];
					lat_inh.rate_std = [];
					lat_inh.R = [];
					lat_inh.R2 = [];
					lat_inh.PSTH = [];
				end
				lat_inh.MTF_shape = MTF_shape;
				lat_inh.BMF = BMF;

				[rate, rate_std, pitch] = plotMTF(params_RM, AN_temp.average_AN_sout, 0);
				AN_lat_inh.rate = rate;
				AN_lat_inh.rate_std = rate_std;
				AN_lat_inh.target = params_RM.target;
				AN_lat_inh.pitch = pitch;
				R = corrcoef(data_RM.rate,rate);
				AN_lat_inh.R = R(1, 2);
				AN_lat_inh.R2 = R(1, 2).^2;
				AN_lat_inh.PSTH = plotNT_PSTH(params_RM, AN_temp.an_sout, 0);
			end
			if isfield(params_RM, 'stim')
				params_RM = rmfield(params_RM, 'stim');
			end
			% Save model
			filename = [putative '_' model_type '_RM_actual.mat'];
			if strcmp(model_type, 'Energy')
				save(fullfile('C:\DataFiles_JBF\Nat-Timbre\data\manuscript',...
					'energy_model', filename), 'params_RM', 'energy')
			elseif strcmp(model_type, 'SFIE')
				save(fullfile('C:\DataFiles_JBF\Nat-Timbre\data\manuscript\',...
					'SFIE_model', filename), 'params_RM', 'SFIE', 'model_params')
			elseif strcmp(model_type, 'Lat_Inh')
				save(fullfile('C:\DataFiles_JBF\Nat-Timbre\data\manuscript\',...
					'lat_inh_model', filename), 'params_RM', 'lat_inh', 'model_params')
			end

		end
	end
	elapsedTime = toc(timerVal)/60;
	disp([putative ' Model took ' num2str(elapsedTime) ' minutes'])
end
