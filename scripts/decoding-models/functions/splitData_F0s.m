function [T_train, T_test, F0s_test] = splitData_F0s(T, F0s)

% Your data parameters
num_F0s = length(F0s);
if num_F0s == 16
	reps_per_F0 = 40;
else
	reps_per_F0 = 20;
end
num_test_F0s = floor(num_F0s*0.2);  % 20% of 40 F0s

% Create F0 labels for each sample
F0_labels = repelem(1:num_F0s, reps_per_F0);  % [1,1,1...1,2,2,2...2,...,40,40,40...40]

% Method 1: Random selection of F0s for test set
test_F0s = randperm(num_F0s, num_test_F0s);  % Randomly select 8 F0s
train_F0s = setdiff(1:num_F0s, test_F0s);   % Remaining 32 F0s

% Create logical indices for train/test split
test_indices = ismember(F0_labels, test_F0s);
train_indices = ~test_indices;

% Split your data (assuming your table is called T)
T_train = T(train_indices, :);
T_test = T(test_indices, :);
F0s_test = F0s(test_F0s);

% fprintf('Training set: %d samples from %d F0s\n', sum(train_indices), length(train_F0s));
% fprintf('Test set: %d samples from %d F0s\n', sum(test_indices), length(test_F0s));
% fprintf('Test F0s: %s\n', mat2str(sort(round(10.^F0s(test_F0s)))));

end