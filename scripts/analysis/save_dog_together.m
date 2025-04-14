%% save_dog_together
clear 

%% Load 

%modelpath = '/Volumes/Nat-Timbre/data/manuscript/dog_model';
modelpath = 'C:\DataFiles_JBF\Nat-Timbre\data\manuscript\dog_model';
list = dir(modelpath);

filenames = {list(endsWith({list.name}, '.mat')).name};
num_files = size(filenames,2);

for ii = 1:num_files 

	% Load in 
	filename = filenames{ii};
	load(fullfile(modelpath, filename))

	% Struct to save out all data and fits
	dog_analysis(ii).putative = dog_gauss_analysis.putative;
	dog_analysis(ii).dog_predicted = dog_gauss_analysis.dog_predicted;
	dog_analysis(ii).gaus_predicted = dog_gauss_analysis.gaus_predicted;
	dog_analysis(ii).CF = dog_gauss_analysis.CF;
	dog_analysis(ii).rate = dog_gauss_analysis.rate;
	dog_analysis(ii).R2_dog = dog_gauss_analysis.R2_dog;
	dog_analysis(ii).R2_gauss = dog_gauss_analysis.R2_gauss;
	dog_analysis(ii).pitch = dog_gauss_analysis.pitch;
	dog_analysis(ii).spont = dog_gauss_analysis.spont;
	dog_analysis(ii).rate_std = dog_gauss_analysis.rate_std;
	dog_analysis(ii).p_value = dog_gauss_analysis.p_value;
	if isfield(dog_gauss_analysis, 'dog_params')
		dog_analysis(ii).dog_params = dog_gauss_analysis.dog_params;
	end 
	if isfield(dog_gauss_analysis, 'gauss_params')
		dog_analysis(ii).gauss_params = dog_gauss_analysis.gauss_params;
	end

	% Save variables 
	R2_dog_all(ii) = dog_gauss_analysis.R2_dog;
	R2_gauss_all(ii) = dog_gauss_analysis.R2_gauss;

end

%% Save file

%datapath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/';
datapath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data';
save(fullfile(datapath, 'dog_analysis2.mat'), "dog_analysis", "R2_gauss_all", "R2_dog_all")

