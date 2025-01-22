%% run_AN_model

%% Creates filename for analyzed data
if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/AN_neurogram';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\AN_neurogram';
end

%% Set up AN model to run for all stimuli 

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

%% Save AN

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
model_params.nrep = 10; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)

nstim = params.num_stim;
for ii = 1:nstim
	timerVal = tic;

	% Get single stimulus
	params.stim = params.stim(ii, :);
	params.num_stim = ii;
	target = extractBefore(params.mlist(ii).wav_file, '.');
	F0 = extractBetween(params.mlist(ii).wav_file, 'ff.', '.wav');
	savefile = sprintf('%s_F0_%s_AN.mat', target, F0{1});

	AN = modelAN(params, model_params); % HSR for IC input
	save(fullfile(savepath, savefile),'AN')

	elapsedTime = toc(timerVal)/60;
	disp(['Model took ' num2str(elapsedTime) ' minutes'])
end

