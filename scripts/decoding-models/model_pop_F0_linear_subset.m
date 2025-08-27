%% Load in data
clear 

base = getPathsNT();
load(fullfile(base, 'model_comparisons', 'Data_NT_3.mat'), 'nat_data')
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning

%% Run linear model

targets = {'Bassoon', 'Oboe', 'Invariant'};
for itarget = 1:3
	timerVal = tic;

	% Get data and F0s
	target = targets{itarget};
	F0s = getF0s(target);
	F0s = log10(F0s);
	[sesh, num_data] = getF0Sessions(nat_data, target);
	T = getF0PopTable(nat_data, target, sesh, F0s, num_data, 'linear', 'Rate');

	% %% Test all parameter combinations
	% % Split into training/testing
	% [T_train, T_test] = splitData(T);
	% [best_model, best_params, all_results] = optimize_linear_model(T_train, T_test);

	parfor irep = 1:500

		% Split into training/testing
		% [T_train, T_test] = splitData_Reps(T);
		% F0s_test = F0s;

		% % Split data to train on only some F0s
		[T_train, T_test ,F0s_test] = splitData_F0s(T, F0s);

		% % Split data to train on F0s for one instrument
		% train_instrument = 'Oboe';
		% [T_train, T_test] = splitData_Instru(T, train_instrument);

		% Define hyperparameter search space
		lambda_values = logspace(-3, -0.5, 30);  % Range of regularization values
		best_lambda = [];
		best_cv_mse = inf;

		% Cross-validation for each lambda
		for lambda = lambda_values
			Mdl_cv = fitrlinear(T_train, 'Response', 'BetaTolerance', 1.0e-06, ...
				'Learner', 'leastsquares', 'Lambda', lambda, 'Solver', 'lbfgs', ...
				'KFold', 5, 'CrossVal', 'on', 'Regularization', 'ridge');
			cv_mse = kfoldLoss(Mdl_cv);
			if cv_mse < best_cv_mse
				best_cv_mse = cv_mse;
				best_lambda = lambda;
			end
		end

		%fprintf('Best lambda: %.6f, CV MSE: %.5f\n', best_lambda, best_cv_mse);

		% Train final model with best hyperparameters
		Mdl_final = fitrlinear(T_train, 'Response', 'BetaTolerance', 1.0e-06, ...
			'Learner', 'leastsquares', 'Lambda', best_lambda, 'Solver', 'lbfgs', ...
			'Regularization', 'ridge');

		% Test set evaluation
		pred_F0_test = predict(Mdl_final, T_test);
		test_mse = mean((T_test.Response - pred_F0_test).^2);
		%fprintf('Test MSE with best lambda: %.5f\n', test_mse);
		r = corrcoef(pred_F0_test, T_test.Response);
		r2 = r(1, 2)^2;

		% Get vals
		results(itarget, irep).target = target;
		results(itarget, irep).predicted = pred_F0_test;
		results(itarget, irep).actual = T_test.Response;
		results(itarget, irep).test_mse = test_mse;
		results(itarget, irep).r2 = r2;

		% Plot
		% figure
		% tiledlayout(1, 3)
		% nexttile
		% scatter(10.^(T_test.Response), 10.^(pred_F0_test), 20, 'filled', ...
		% 	'MarkerFaceAlpha',0.5, 'MarkerFaceColor','k')
		% set(gca, 'xscale', 'log', 'yscale', 'log')
		% ylim(10.^[F0s(1) F0s(end)])
		% xlim(10.^[F0s(1) F0s(end)])
		% xticks([100 200 500 1000 1500])
		% hold on
		% plot(10.^T_test.Response, 10.^T_test.Response, 'k')
		% plot(T_test.Response, T_test.Response)
		% title(['R^2 = ' num2str(r2)])
		% xlabel('Actual F0s (Log10(F0s))')
		% ylabel('Predicted F0s (Log10(F0s))')
		% grid on

		% closest_cat = [];
		% for ii = 1:length(pred_F0_test)
		% 	differences = abs(F0s - pred_F0_test(ii));
		% 	[~, closest_column_index] = min(differences);
		% 	closest_cat(ii) = F0s(closest_column_index);
		% end

		% % Plot confusion matrix
		% %nexttile
		% C = confusionmat(T_test.Response, closest_cat);
		% %confusionchart(C)
		%
		% % Calculate accuracy
		% chart = confusionchart(T_test.Response,closest_cat); % Generate confusion chart
		% confusionMatrix = chart.NormalizedValues; % Get the normalized confusion matrix
		% accuracy = sum(diag(confusionMatrix)) / sum(confusionMatrix(:)); % Calculate accuracy
		% title(sprintf('Accuracy = %0.2f%%', accuracy*100))
		%
		% % Residuals analysis
		% nexttile
		% residuals = T_test.Response - pred_F0_test;
		% scatter(pred_F0_test, residuals);
		% xlabel('Predicted Values');
		% ylabel('Residuals');
	end
	times = toc(timerVal);
	fprintf('%d/%d, took %0.02f minutes\n', itarget, 3, times/60)
end

save(fullfile(base, 'model_comparisons', 'Pop_Rate_F0_Linear_F0s.mat'),...
	"results")

%% Plotssss 

r2_bassoon = [results(1,:).r2];
r2_oboe = [results(2,:).r2];
r2_invariant = [results(3,:).r2];
mse_bassoon = [results(1,:).test_mse];
mse_oboe = [results(2,:).test_mse];
mse_invariant = [results(3,:).test_mse];

figure
nexttile
hold on
swarmchart(ones(500,1), r2_bassoon, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
swarmchart(ones(500,1)*2, r2_oboe, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
swarmchart(ones(500,1)*3, r2_invariant, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)

nexttile
hold on
swarmchart(ones(500,1), mse_bassoon, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
swarmchart(ones(500,1)*2, mse_oboe, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
swarmchart(ones(500,1)*3, mse_invariant, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)

%% Train on one instrument, test on the other instrument

% Choose dataset to model

% Get data and F0s
F0s = getF0s(target);
F0s = log10(F0s);
[sesh, num_data] = getF0Sessions(nat_data, target);
T = getF0PopTable(nat_data, target, sesh, F0s, num_data, 'linear', 'Rate');
CFs = [nat_data(sesh).CF];
MTFs = {nat_data(sesh).MTF};

targets = {'Bassoon', 'Oboe'};
for itarget = 1:2
	timerVal = tic;

	train_instrument = targets{itarget};
	parfor irep = 1:500

		% % Split data to train on F0s for one instrument
		[T_train, T_test] = splitData_Instru(T, train_instrument);

		% Define hyperparameter search space
		lambda_values = logspace(-3, -0.5, 30);  % Range of regularization values
		best_lambda = [];
		best_cv_mse = inf;

		% Cross-validation for each lambda
		for lambda = lambda_values
			Mdl_cv = fitrlinear(T_train, 'Response', 'BetaTolerance', 1.0e-06, ...
				'Learner', 'leastsquares', 'Lambda', lambda, 'Solver', 'lbfgs', ...
				'KFold', 5, 'CrossVal', 'on', 'Regularization', 'ridge');

			cv_mse = kfoldLoss(Mdl_cv);

			if cv_mse < best_cv_mse
				best_cv_mse = cv_mse;
				best_lambda = lambda;
			end
		end

		%fprintf('Best lambda: %.6f, CV MSE: %.5f\n', best_lambda, best_cv_mse);

		% Train final model with best hyperparameters
		Mdl_final = fitrlinear(T_train, 'Response', 'BetaTolerance', 1.0e-06, ...
			'Learner', 'leastsquares', 'Lambda', best_lambda, 'Solver', 'lbfgs', ...
			'Regularization', 'ridge');

		% Test set evaluation
		pred_F0_test = predict(Mdl_final, T_test);
		test_mse = mean((T_test.Response - pred_F0_test).^2);
		%fprintf('Test MSE with best lambda: %.5f\n', test_mse);
		r = corrcoef(pred_F0_test, T_test.Response);
		r2 = r(1, 2)^2;

		% Get vals
		results(itarget,irep).train_instrument = train_instrument;
		results(itarget,irep).predicted = pred_F0_test;
		results(itarget,irep).actual = T_test.Response;
		results(itarget,irep).test_mse = test_mse;
		results(itarget,irep).r2 = r2;
	end
	times = toc(timerVal);
	fprintf('%d/%d, took %0.02f minutes\n', itarget, 2, times/60)
end

save(fullfile(base, 'model_comparisons', 'Pop_Rate_F0_Linear_Instruments.mat'),...
	"results")

%% Plotsss 

r2_bassoon_train = [results(1,:).r2];
r2_oboe_train = [results(2,:).r2];
mse_bassoon_train = [results(1,:).test_mse];
mse_oboe_train = [results(2,:).test_mse];

figure
nexttile
hold on
swarmchart(ones(500,1), r2_bassoon_train, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
swarmchart(ones(500,1)*2, r2_oboe_train, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)

nexttile
hold on
swarmchart(ones(500,1), mse_bassoon_train, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
swarmchart(ones(500,1)*2, mse_oboe_train, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)