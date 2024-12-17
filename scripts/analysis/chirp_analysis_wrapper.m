%% chirp_analysis_wrapper.m
%
% Loads in natural timbre stimuli and runs through chirp_analysis 
%
% 
% 
clear 

%% Try vowel example
% 
% firstparts = [{'/Users/jfritzinger/Physio/DCP_naq\centerVoweltokens\F0095M03'},...
% 	{'/Users/jfritzinger/Physio/DCP_naq\centerVoweltokens\F0202W39'}];
% finalpart = {'48828center200ms.wav'};
% vow_F0s = [95, 202];
% all_twelve_ticks = [{'IY'}, {'IH'}, {'EI'}, {'EH'}, {'AE'}, {'AH'},...
% 	{'AW'}, {'OO'}, {'UW'}, {'ER'}, {'OA'}, {'UH'}];
% 
% vow_str = all_twelve_ticks{7};
% vow_F0 = vow_F0s(1);
% file = append(firstparts{1},vow_str,finalpart{1});
% 
% vowel_decomp(vow_str,vow_F0,file)


%% Get list of all timbre stimuli (bassoon)

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning 

target = 'Oboe';
listing = dir(fullfile(fpath, 'waveforms', '*.wav'));
target_WAV = arrayfun(@(n) contains(listing(n).name, target), 1:numel(listing), 'UniformOutput', false);
wav_nums =  find(cell2mat(target_WAV));

d = dir(fullfile(fpath,'waveforms', '*.wav'));
all_files = sort({d.name});
nfiles = length(wav_nums);
wav_npts = zeros(1,nfiles);
wav_data = cell(1,nfiles);

for i = 1:nfiles
   files{1,i} = all_files{wav_nums(i)}; 
end
% files = all_files;
% nfiles = size(all_files, 2);

% Sort by frequency of pitch
index = [];
note_names = extractBetween(files, 'ff.','.');

for ii = 1:nfiles % Find index of each note in tuning spreadsheet
    index(ii) = find(strcmp(note_names(ii), tuning.Note));
end
pitch_order = tuning.Frequency(index); % Get freqs of each note
F0s = pitch_order;

%% 

% Loop through all stimuli
for ind = 1:nfiles

	target = extractBefore(files{ind}, '.');

	target_F0 = F0s(ind);
	target_file = fullfile(fpath,'waveforms', files{ind});


	chirp_analysis(target,target_F0,target_file)

end







