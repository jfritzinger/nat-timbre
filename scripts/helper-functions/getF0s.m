function F0s = getF0s(target)

[base, ~, ~, ~] = getPathsNT();
tuning = readtable(fullfile(base, 'Tuning.xlsx')); % Load in tuning

if strcmp(target, 'Invariant') || strcmp(target, 'Timbre')

	listing = dir(fullfile(base, 'waveforms', ['*' 'Bassoon' '*.wav']));
	files = {listing.name};
	note_names = extractBetween(files, 'ff.', '.');
	[~, index] = ismember(note_names, tuning.Note);
	F0s1 = tuning.Frequency(index);
	[F0s, ~] = sort(F0s1);
	ind_b = 25:40;
	F0s = F0s(ind_b);
else
	listing = dir(fullfile(base, 'waveforms', ['*' target '*.wav']));
	files = {listing.name};
	note_names = extractBetween(files, 'ff.', '.');
	[~, index] = ismember(note_names, tuning.Note);
	F0s1 = tuning.Frequency(index);
	[F0s, ~] = sort(F0s1);
end

end
