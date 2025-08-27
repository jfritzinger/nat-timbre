function [T_train, T_test] = splitData_Reps(T)

% Determine number of rows for training and testing data
nrows = size(T, 1);
test_rows = nrows*0.2;
train_rows = nrows*0.8;

% Preallocate indices
trainIndices = zeros(train_rows, 1);
testIndices = zeros(test_rows, 1);

% Create indices using systematic selection
ngroups = nrows/20; % 20 repetitions of each stimulus
for group = 0:ngroups-1

    % Calculate group range (1-20, 21-40, etc.)
    startIdx = group*20 + 1;
    endIdx = (group+1)*20;
    
    % Random permutation within current group
    shuffled = randperm(20) + startIdx - 1;
    
    % Store indices
    trainStart = group*16 + 1;
    trainEnd = (group+1)*16;
    testStart = group*4 + 1;
    testEnd = (group+1)*4;
    
    trainIndices(trainStart:trainEnd) = shuffled(1:16);
    testIndices(testStart:testEnd) = shuffled(17:20);
end

% Create training and testing sets
T_train = T(trainIndices, :);
T_test = T(testIndices, :);

% Verify sizes
% disp(['Training size: ', num2str(size(T_train))]);
% disp(['Testing size: ', num2str(size(T_test))]);

end