function [T_train, T_test] = splitData_Instru(T, train_instrument)

% Create logical indices for train/test split
if strcmp(train_instrument, 'Bassoon') % Train with bassoon
	pattern = [ones(1, 20), zeros(1, 20)];  % [1,1,1...1,0,0,0...0] (20 ones, 20 zeros)
	train_indices = repmat(pattern, 1, 16);     % Repeat 16 times to get 640 elements
else % Train with oboe
	pattern = [zeros(1, 20), ones(1, 20)];
	train_indices = repmat(pattern, 1, 16);
end
test_indices = ~train_indices;

% Split your data (assuming your table is called T)
T_train = T(logical(train_indices), :);
T_test = T(test_indices, :);

end