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

%model_type = 'SFIE';
%model_type = 'Energy';
model_type = 'Lat_Inh';

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
for isesh = 46:num_index
	timerVal = tic;

	% Load in session
	putative = sessions.Putative_Units{indices(isesh)};
	CF = sessions.CF(indices(isesh));
	load(fullfile(datapath, [putative '.mat']))
	MTF_shape = sessions.MTF{indices(isesh)};
	params_NT = data(6:7, 2);

	if strcmp(MTF_shape, 'BS')
		BMF = sessions.WMF(indices(isesh));
	elseif strcmp (MTF_shape, 'BE')
		BMF = sessions.BMF(indices(isesh));
	else
		BMF = 100;
	end

	if strcmp(model_type, 'SFIE')
		AN = cell(2, 1);
		SFIE = cell(2,1);
	elseif strcmp(model_type, 'Energy')
		Fs = 100000;
		gamma_param.srate = Fs;
		energy = cell(2, 1);
		gamma_param.fc = CF;
	elseif strcmp(model_type, 'Lat_Inh')
		AN_lat_inh = cell(2,1);
		lat_inh = cell(2,1);
	end

	for iNT = 1:2
	
		% Analysis
		if isempty(params_NT{iNT})
			% Did not record natural timbre response for this instrument
		else
			data_NT = analyzeNT(params_NT{iNT});

			if isempty(data_NT.rate)

			else

			% Set up stimuli 
			params_NT{iNT}.Fs = 100000;
			params_NT{iNT}.mnrep = 20;
			params_NT{iNT}.dur = 0.2;
			if ~isfield(params_NT{iNT}, 'target')
				params_NT{iNT}.target = extractBefore(params_NT{iNT}.list(1).wav_file, '.');
			end
			[params_NT{iNT}] = generate_NT(params_NT{iNT});
			params_NT{iNT}.num_stim = size(params_NT{iNT}.stim, 1);

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
				AN_temp = modelAN(params_NT{iNT}, model_params); % HSR for IC input
				SFIE_temp = wrapperIC(AN_temp.an_sout, params_NT{iNT}, model_params); % SFIE output

				% Plot
				if strcmp(MTF_shape, 'BS')
					[rate, rate_std, pitch] = plotNT(params_NT{iNT}, SFIE_temp.average_ic_sout_BS, 0);
					SFIE{iNT}.rate = rate;
					SFIE{iNT}.rate_std = rate_std;
					SFIE{iNT}.target = params_NT{iNT}.target;
					SFIE{iNT}.pitch = pitch;
					R = corrcoef(data_NT.rate,rate); % Correlation
					SFIE{iNT}.R = R(1, 2);
					SFIE{iNT}.R2 = R(1, 2).^2;
					SFIE{iNT}.PSTH = plotNT_PSTH(params_NT{iNT}, SFIE_temp.ic_BS, 0);
					SFIE{iNT}.SFIE_temp = SFIE_temp;

				elseif strcmp(MTF_shape, 'BE')
					[rate, rate_std, pitch] = plotNT(params_NT{iNT}, SFIE_temp.average_ic_sout_BE, 0);
					SFIE{iNT}.rate = rate;
					SFIE{iNT}.rate_std = rate_std;
					SFIE{iNT}.target = params_NT{iNT}.target;
					SFIE{iNT}.pitch = pitch;
					R = corrcoef(data_NT.rate,rate); % Correlation
					SFIE{iNT}.R = R(1, 2);
					SFIE{iNT}.R2 = R(1, 2).^2;
					SFIE{iNT}.PSTH = plotNT_PSTH(params_NT{iNT}, SFIE_temp.ic_BE, 0);
					SFIE{iNT}.SFIE_temp = SFIE_temp;
				else
					SFIE{iNT}.rate = [];
					SFIE{iNT}.rate_std = [];
					SFIE{iNT}.R = [];
					SFIE{iNT}.R2 = [];
					SFIE{iNT}.PSTH = [];
					SFIE{iNT}.SFIE_temp = [];
				end
				SFIE{iNT}.MTF_shape = MTF_shape;
				SFIE{iNT}.BMF = BMF;

				[rate, rate_std, pitch] = plotNT(params_NT{iNT}, AN_temp.average_AN_sout, 0);
				AN{iNT}.rate = rate;
				AN{iNT}.rate_std = rate_std;
				AN{iNT}.target = params_NT{iNT}.target;
				AN{iNT}.pitch = pitch;
				R = corrcoef(data_NT.rate,rate);
				AN{iNT}.R = R(1, 2);
				AN{iNT}.R2 = R(1, 2).^2;
				AN{iNT}.PSTH = plotNT_PSTH(params_NT{iNT}, AN_temp.an_sout, 0);
				AN{iNT}.AN_temp = AN_temp;

			elseif strcmp(model_type, 'Energy') % Energy model

				stimulus = [params_NT{iNT}.stim zeros(size(params_NT{iNT}.stim,1),0.1*Fs)];
				tvals = (1:length(stimulus))/Fs;
				gamma_IF_reg = zeros(1,length(tvals));
				impaired = 0; % 0 = not impaired; 1 = 'impaired'
				pin_gamma = zeros(size(stimulus, 1), Fs*params_NT{iNT}.dur+0.1*Fs);
				for istim = 1:size(stimulus, 1)
					pin_gamma(istim,:) = gamma_filt(stimulus(istim,:),gamma_param,impaired);
				end
				pin_gamma = pin_gamma(:,1:params_NT{iNT}.dur*Fs);
				energ_out = sqrt(mean(pin_gamma.^2,2));
				[rate, rate_std, pitch] = plotNT(params_NT{iNT}, energ_out, 0);
				R_int = corrcoef(data_NT.rate,rate);
				energy{iNT}.PSTH = plotNT_PSTH(params_NT{iNT}, pin_gamma, 0);

				energy{iNT}.energ_out = energ_out;
				energy{iNT}.rate = rate;
				energy{iNT}.rate_std = rate_std;
				energy{iNT}.target = params_NT{iNT}.target;
				energy{iNT}.pitch = pitch;
				energy{iNT}.R = R_int(1,2);
				energy{iNT}.R2 =  R_int(1, 2).^2;

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
				AN_temp = modelLateralAN(params_NT{iNT}, model_params);
				latinh_temp = modelLateralSFIE(params_NT{iNT}, model_params,...
					AN_temp.an_sout, AN_temp.an_sout_lo, AN_temp.an_sout_hi,...
					'CS_params', lm_params);

				if strcmp(MTF_shape, 'BS') || strcmp(MTF_shape, 'BE')
					[rate, rate_std, pitch] = plotNT(params_NT{iNT}, latinh_temp.avIC, 0);
					lat_inh{iNT}.rate = rate;
					lat_inh{iNT}.rate_std = rate_std;
					lat_inh{iNT}.target = params_NT{iNT}.target;
					lat_inh{iNT}.pitch = pitch;
					R = corrcoef(data_NT.rate,rate); % Correlation
					lat_inh{iNT}.R = R(1, 2);
					lat_inh{iNT}.R2 = R(1, 2).^2;
					lat_inh{iNT}.PSTH = plotNT_PSTH(params_NT{iNT}, latinh_temp.ic, 0);
					lat_inh{iNT}.latinh_temp = latinh_temp;
				else
					lat_inh{iNT}.rate = [];
					lat_inh{iNT}.rate_std = [];
					lat_inh{iNT}.R = [];
					lat_inh{iNT}.R2 = [];
					lat_inh{iNT}.PSTH = [];
					lat_inh{iNT}.latinh_temp = [];
				end
				lat_inh{iNT}.MTF_shape = MTF_shape;
				lat_inh{iNT}.BMF = BMF;

				[rate, rate_std, pitch] = plotNT(params_NT{iNT}, AN_temp.average_AN_sout, 0);
				AN_lat_inh{iNT}.rate = rate;
				AN_lat_inh{iNT}.rate_std = rate_std;
				AN_lat_inh{iNT}.target = params_NT{iNT}.target;
				AN_lat_inh{iNT}.pitch = pitch;
				R = corrcoef(data_NT.rate,rate);
				AN_lat_inh{iNT}.R = R(1, 2);
				AN_lat_inh{iNT}.R2 = R(1, 2).^2;
				AN_lat_inh{iNT}.PSTH = plotNT_PSTH(params_NT{iNT}, AN_temp.an_sout, 0);
			end	
			end
			if isfield(params_NT{iNT}, 'stim')
				params_NT{iNT} = rmfield(params_NT{iNT}, 'stim');
			end
		end
	end

	% Save model
	filename = [putative '_' model_type '20reps.mat'];
	if strcmp(model_type, 'Energy')
		save(fullfile('C:\DataFiles_JBF\Nat-Timbre\data\manuscript',...
			'energy_model', filename), 'params_NT', 'energy')
	elseif strcmp(model_type, 'SFIE')
		save(fullfile('C:\DataFiles_JBF\Nat-Timbre\data\manuscript\',...
			'SFIE_model', filename), 'params_NT', 'SFIE', 'model_params')
	elseif strcmp(model_type, 'Lat_Inh')
		save(fullfile('C:\DataFiles_JBF\Nat-Timbre\data\manuscript\',...
			'lat_inh_model', filename), 'params_NT', 'lat_inh', 'model_params')
	end
	elapsedTime = toc(timerVal)/60;
	disp([putative ' Model took ' num2str(elapsedTime) ' minutes'])
end
