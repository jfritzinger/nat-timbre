%% chirp_analysis_PDF.m
%
% Script that loads in the decomposed natural timbre stimuli and plots them
% into a PDF for easy analysis
%
% Author: J. Fritzinger
% Created: 2022-10-17; Last revision: 2024-10-17
%
% -------------------------------------------------------------------------
clear
import mlreportgen.dom.*
import mlreportgen.report.*

%% Get list of all timbre stimuli

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning 

target = 'Bassoon';
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
[~, order] = sort(pitch_order); % Sort freqs
F0s = pitch_order;
files = files(order);
F0s = F0s(order);

%% Set up report

% Initialize report
filename = 'Test';
images = {}; %hold all plots as images, need to delete when finished
datetime.setDefaultFormats('default','yyyy-MM-dd_hhmmss')
if ismac
	report_name = sprintf('%s%s_%s.pdf', '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/', datetime, filename);
else
	report_name = sprintf('%s%s_%s.pdf', 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\figures\', datetime, filename);
end
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.01in';
pm.PageMargins.Header = '0.01in';
pm.PageMargins.Bottom = '0.01in';
pm.PageMargins.Footer = '0.01in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';

%% Create PDF 

% Loop through all stimuli
for ind = 1 %:nfiles

	target = extractBefore(files{ind}, '.');
	target_F0 = F0s(ind);
	target_file = fullfile(fpath,'waveforms', files{ind});

	% Create plots 
	analyzeChirps(target,target_F0,target_file)
	plotChirpSpectrogram(target,target_F0,target_file)

	% Label
	h = Heading(2, sprintf("%s F0=%0.0f Hz\n", target, target_F0));
	b = Border();
	b.BottomStyle = 'single';
	b.BottomColor = 'LightGray';
	b.BottomWidth = '1pt';
	h.Style = [h.Style {Color('Black'), b}, {PageBreakBefore()}];
	append(rpt,h);


	% Add to PDF
	fig1 = gcf;
	[plt1, images] = addtoSTPDF(images, fig1, [num2str(round(target_F0)) '_1']);
	append(rpt, plt1);

	fig2 = gcf;
	[plt2, images] = addtoSTPDF(images, fig2, [num2str(round(target_F0)) '_2']);
	append(rpt, plt2);

end

%% Save 

% Closes and opens PDF to view
close(rpt);
for i = 1:length(images)
	delete(images{1,i}.Path);
end
rptview(rpt)


%% FUNCTIONS

function [img, images] = addtoSTPDF(images, fig, title)
import mlreportgen.dom.*

% Set figure size, recommended
if contains(title, '_1')
	values = [8, 4];
	fig.PaperSize = values;
	fig.PaperPosition = [0 0 values];
	fig.Units = 'inches';
	fig.Position(3:4) = values;
else
	values = [8, 6];
	fig.PaperSize = values;
	fig.PaperPosition = [0 0 values];
	fig.Units = 'inches';
	fig.Position(3:4) = values;
end


% Add the plot to the document
name = sprintf('%s.svg', title);
print(fig, name, '-dsvg');
img = Image(name);
delete(fig) %delete plot figure window
images = [images {img}];

end