function [best_model, best_params, results] = optimize_linear_model(T_train, T_test)
    % Define search space
    regularizations = {'ridge', 'lasso'};
    lambda_values = [0, logspace(-6, 2, 20)];
    solvers = {'lbfgs', 'sgd'};
    beta_tolerances = [1e-6, 1e-4, 1e-3];
    
    % Initialize results storage
    results = [];
    best_cv_mse = inf;
    best_params = struct();
    best_model = [];
    
    fprintf('Testing %d parameter combinations...\n', ...
        length(regularizations) * length(lambda_values) * ...
        length(solvers) * length(beta_tolerances));
    
    combo_count = 0;
    for reg = regularizations
        for lambda = lambda_values
            for solver = solvers
                for beta_tol = beta_tolerances
                    combo_count = combo_count + 1;
                    
                    try
                        % Skip incompatible combinations
                        if strcmp(solver{1}, 'sgd') && lambda == 0
                            continue; % SGD needs regularization
                        end
                        
                        % Cross-validation
                        Mdl_cv = fitrlinear(T_train, 'Response', ...
                            'Regularization', reg{1}, ...
                            'Lambda', lambda, ...
                            'Solver', solver{1}, ...
                            'BetaTolerance', beta_tol, ...
                            'Learner', 'leastsquares', ...
                            'KFold', 5, 'CrossVal', 'on');
                        
                        cv_mse = kfoldLoss(Mdl_cv);
                        
                        % Store results
                        result = struct('regularization', reg{1}, ...
                                      'lambda', lambda, ...
                                      'solver', solver{1}, ...
                                      'beta_tolerance', beta_tol, ...
                                      'cv_mse', cv_mse);
                        results = [results; result];
                        
                        % Update best model
                        if cv_mse < best_cv_mse
                            best_cv_mse = cv_mse;
                            best_params = result;
                        end
                        
                        if mod(combo_count, 10) == 0
                            fprintf('Completed %d combinations, best CV MSE so far: %.6f\n', ...
                                combo_count, best_cv_mse);
                        end
                        
                    catch ME
                        fprintf('Failed combination: %s, lambda=%.1e, %s, tol=%.1e\n', ...
                            reg{1}, lambda, solver{1}, beta_tol);
                    end
                end
            end
        end
    end
    
    % Train best model on full training set
    best_model = fitrlinear(T_train, 'Response', ...
        'Regularization', best_params.regularization, ...
        'Lambda', best_params.lambda, ...
        'Solver', best_params.solver, ...
        'BetaTolerance', best_params.beta_tolerance, ...
        'Learner', 'leastsquares');
    
    % Test set evaluation
    pred_test = predict(best_model, T_test);
    test_mse = mean((T_test.Response - pred_test).^2);
    
    fprintf('\nBest Parameters:\n');
    fprintf('Regularization: %s\n', best_params.regularization);
    fprintf('Lambda: %.6f\n', best_params.lambda);
    fprintf('Solver: %s\n', best_params.solver);
    fprintf('Beta Tolerance: %.1e\n', best_params.beta_tolerance);
    fprintf('CV MSE: %.6f\n', best_params.cv_mse);
    fprintf('Test MSE: %.6f\n', test_mse);
end