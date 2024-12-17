%% save_model_responses
%
% Script that runs either SFIE population, energy, or SFIE population
% models for each stimulus
%
% Author: J. Fritzinger
% Created: 2022-09-15; Last revision: 2024-10-14
%
% -------------------------------------------------------------------------
clear

%% Generate stimulus

% Set up stimuli 
params.Fs = 100000;
params.mnrep = 1;
params.dur = 0.2;
params.target = 'Bassoon';
params.signal_onset_delay = 0;
params.noise_ramp_dur = 0.02;
params.signal_spls = 73;
params.nrep = 1;
params.reptim = 0.6;
params.noise_state = 0;
params.noise_shape = [];
params.noise_spls = 0;
params.signal_ramp_dur = 0.02;
[params] = generate_NT(params);
params.num_stim = size(params.stim, 1);

%% Save SFIE Response
timerVal = tic;

% Model parameters
CFs = logspace(log10(125), log10(10000), 100);
model_params.type = 'SFIE';
model_params.range = 1; % 1 = population model, 2 = single cell model
model_params.species = 2; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.CFs = CFs;
model_params.nAN_fibers_per_CF = 10;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 1; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)

%params.stim = params.stim(1, :);
%params.num_stim = 1;
AN_temp = modelAN(params, model_params); % HSR for IC input
SFIE = wrapperIC(AN_temp.an_sout, params, model_params); % SFIE output

elapsedTime = toc(timerVal)/60;
disp(['Model took ' num2str(elapsedTime) ' minutes'])

%% Lateral inhibition 
timerVal = tic;

% Model params 
CFs = logspace(log10(220), log10(10000), 100);

% Model parameters
model_params.type = 'Lateral Model';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 2; % 1 = cat, 2 = human
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
model_params.config_type = 'BS inhibited by off-CF BS';


% Run model
for iMTF = 2
	if iMTF == 1 % BE
		S = 0.4; % Strength, S =
		D = 0; % Delay, D =
		oct_range = 0.5; % CF range =
	else % BS
		S = 0.25; % Strength, S =
		D = 0; % Delay, D =
		oct_range = 0.75; % CF range =
	end
	lm_params = [S S D];
	for iCF = 1:length(CFs)
		timerVal2 = tic;

		model_params.lateral_CF = [CFs(iCF)*2^(-1*oct_range), CFs(iCF), CFs(iCF)*2^oct_range];
		model_params.CFs = model_params.lateral_CF;
		model_params.CF_range = model_params.CFs(2);
		AN_temp = modelLateralAN(params, model_params);
		latinh_temp = modelLateralSFIE(params, model_params,...
			AN_temp.an_sout, AN_temp.an_sout_lo, AN_temp.an_sout_hi,...
			'CS_params', lm_params);
		ic(:, iCF) = latinh_temp.avIC;

		elapsedTime = toc(timerVal2)/60;
		disp(['Model took ' num2str(elapsedTime) ' minutes'])	
	end
	if iMTF == 1
		latinh.BE = ic;
		latinh.BE_CFs = CFs;
	else
		latinh.BS = ic;
		latinh.BS_CFs = CFs;
	end 
end

elapsedTime = toc(timerVal)/60;
disp(['Model took ' num2str(elapsedTime) ' minutes'])		

% Interpolate data to match the same CF range as previous 
% CFs_ideal = logspace(log10(125), log10(10000), 100);
% BS_rate = interp1(latinh.BS(1,:), latinh.BS_CFs, CFs_ideal, 'linear');
% 
% figure
% plot(latinh.BS_CFs, latinh.BS(1,:))
% hold on
% plot(CFs_ideal, BS_rate)

save('Model_NT_human_latinh.mat', "latinh");


%% Energy calculations 
species = 1; % 1 for cat, 2 for human

stimulus = [params.stim zeros(size(params.stim,1),0.1*params.Fs)];
tvals = (1:length(stimulus))/params.Fs;
gamma_IF_reg = zeros(1,length(tvals));
impaired = 0; % 0 = not impaired; 1 = 'impaired'
gamma_param.srate = params.Fs;
pin_gamma = zeros(size(stimulus, 1), params.Fs*params.dur+0.1*params.Fs);
for iCF = 1:length(CFs)
	gamma_param.fc = CFs(iCF);
	for istim = 1:size(stimulus, 1)
		[pin_gamma(istim,:), erb(iCF)] = gamma_filt(stimulus(istim,:),gamma_param,impaired, species);
	end
	pin_gamma2 = pin_gamma(:,1:params.dur*params.Fs);
	RMS(:,iCF) = sqrt(mean(pin_gamma2.^2,2));
end

%% Save responses 

SFIE_BE = SFIE.average_ic_sout_BE;
SFIE_BS = SFIE.average_ic_sout_BS;
AN = AN_temp.average_AN_sout;
save('Model_NT_cat.mat', "AN", "SFIE_BS", "SFIE_BE", "CFs", "RMS", "erb");

%% Plot 
% 
% figure
% tiledlayout(2, 1)
% 
% % Stimulus spectrum
% % nexttile
% % title('Stimulus')
% % y2 = fft(params.stim);
% % m = abs(y2);
% % mdB = 20*log10(m);
% % f = (0:length(y2)-1).*params.stim/length(y2);
% % mdB(mdB<0) = 0;
% % f(f>params.Fs/2) = [];
% % mdB = mdB(1:length(f))';
% % plot(f, mdB, 'LineWidth', 1.5);
% % xlim([125 10000])
% % yticklabels([])
% % set(gca, 'XScale', 'log')
% % hold on
% % box on
% 
% % AN 
% nexttile
% plot(CFs, AN_temp.average_AN_sout);
% set(gca, 'XScale', 'log')
% xlim([125 10000])
% title('AN Response')
% xlabel('CFs (Hz)')
% ylim([0 250])
% 
% % SFIE 
% nexttile
% plot(CFs, SFIE.average_ic_sout_BS);
% set(gca, 'XScale', 'log')
% xlim([125 10000])
% title('BS SFIE Response')
% xlabel('CFs (Hz)')
% ylim([0 40])




