%% save_dog_v_gaussian_fits
clear

%% Load and initialize

% Load in spreadsheet
[base, datapath, savepath, ppi] = getPaths();
spreadsheet_name = 'PutativeTable.xlsx';
sessions = readtable(fullfile(datapath, 'data-cleaning', spreadsheet_name),...
	'PreserveVariableNames',true);
num_data = size(sessions, 1);

%%

has_data = cellfun(@(s) contains(s, 'R'), sessions.Bassoon);
index = find(has_data);

% Sort by CF
CF_list = sessions.CF(has_data);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);
linewidth = 1;
fontsize = 10;

% Plot each neuron
R2_dog_all = NaN(1, num_sessions);
R2_gauss_all = NaN(1, num_sessions);
for isesh = 1:num_sessions
	ineuron = index(order(isesh)); %indices(isesh)
	if any(has_data(ineuron))

		% Load in data
		putative = sessions.Putative_Units{ineuron};
		CF = sessions.CF(ineuron);
		MTF_shape = sessions.MTF{ineuron};
		load(fullfile(datapath, 'neural_data', [putative '.mat']))

		% Calculate fits and plot

		% RM to get spont
		params_RM = data{2, 2};
		data_RM = analyzeRM(params_RM);
		spont = data_RM.spont;

		% Natural timbre analysis
		params = data(14, 2); % 13 is oboe, 14 is bassoon
		data_NT = analyzeNT(params{1});

		% Generate stimulus
		params{1}.Fs = 100000;
		params{1}.mnrep = 1;
		params{1} = generate_NT(params{1});
		params{1}.num_stim = size(params{1}.stim, 1);
		Fs = 100000;
		observed_rate = data_NT.rate;
		r0 = spont;

		% fmincon
		type = 2; % 1: distance, 2: MSE
		stim = params{1}.stim;
		timerVal = tic;

		% Fit gaussian
		init = [CF, 30, 100]; % Initial guess
		lb = [CF/4, 0, 0]; % Lower bounds
		ub = [CF*4, 5000, Inf]; % Upper bounds

		options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-12, ...
			'MaxFunEvals', 10^12, 'maxiterations', 1000, 'ConstraintTolerance', 1e-12, ...
			'StepTolerance', 1e-16, 'display', 'off');
		gaussian_params = fmincon(@(p) ...
			dog_objective_function(p, 'gaussian', Fs, stim, observed_rate, r0, type), ...
			init, [], [], [], [], lb, ub, [], options);
		f = linspace(0, Fs/2, 100000);
		nstim = size(stim, 1);
		gaus_predicted = zeros(nstim, 1);
		for i = 1:nstim
			fc = gaussian_params(1);
			sigma = gaussian_params(2);
			g = gaussian_params(3);
			W = gaussian_model(f, fc, sigma, g);
			gaus_predicted(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
		end

		% Fit DoG model
		dog_init = [20000, 10000, 100, 500,  CF, CF]; % Initial guess
		dog_lb = [100,   100,     10,  10,   CF/4, CF/4]; % Lower bounds
		dog_ub = [30000, 30000,   1000,1000, CF*4, CF*4]; % Upper bounds

		options = optimoptions('fmincon', 'Algorithm','sqp','TolX', 1e-12, ...
			'MaxFunEvals', 10^12, 'maxiterations', 1000, 'ConstraintTolerance', 1e-12, ...
			'StepTolerance', 1e-16, 'display', 'off');
		dog_params = fmincon(@(p) dog_objective_function(p, 'dog', Fs, stim, observed_rate, r0, type), ...
			dog_init, [], [], [], [], dog_lb, dog_ub, [], options);
		disp(['Model took ' num2str(toc(timerVal)) ' seconds'])
		f = linspace(0, Fs/2, 100000);
		nstim = size(stim, 1);
		dog_predicted = zeros(nstim, 1);
		for i = 1:nstim
			W = dog_model(f, dog_params);
			dog_predicted(i) = compute_firing_rate(stim(i, :), Fs, W, f, r0);
		end

		% Get R^2 for all
		gaussian_adj_r_squared = calculate_adj_r_squared(observed_rate,...
			gaus_predicted, 3);
		dog_adj_r_squared = calculate_adj_r_squared(observed_rate,...
			dog_predicted, 6);
		R2_dog_all(isesh) = dog_adj_r_squared;
		R2_gauss_all(isesh) = gaussian_adj_r_squared;

		% Struct to save out all data and fits 
		dog_analysis(isesh).dog_predicted = dog_predicted;
		dog_analysis(isesh).gaus_predicted = gaus_predicted;
		dog_analysis(isesh).CF = CF;
		dog_analysis(isesh).rate = observed_rate;
		dog_analysis(isesh).R2_dog = dog_adj_r_squared;
		dog_analysis(isesh).R2_gauss = gaussian_adj_r_squared;
		dog_analysis(isesh).pitch = data_NT.pitch;
		dog_analysis(isesh).pitch_num = data_NT.pitch_num;
		dog_analysis(isesh).target = params{1}.target;
		dog_analysis(isesh).putative = putative;

	end
end

save('R2_DOG.mat', "R2_gauss_all", "R2_dog_all")

