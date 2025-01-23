%% PDF_AN_neurograms 
clear
import mlreportgen.dom.*
import mlreportgen.report.*

%% Load in natural timbre file names 

if ismac
	fpath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data';
else
	fpath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\';
end
tuning = readtable(fullfile(fpath, 'Tuning.xlsx')); % Load in tuning 
target = 'Bassoon';

% Get all .wav files containing the target instrument name
listing = dir(fullfile(fpath, 'waveforms', ['*' target '*.wav']));
files = {listing.name};

% Extract note names and find corresponding frequencies
note_names = extractBetween(files, 'ff.', '.');
[~, index] = ismember(note_names, tuning.Note);
F0s = tuning.Frequency(index);

% Sort files and frequencies by pitch
[F0s, order] = sort(F0s);
files = files(order);

% Initialize variables for later use (if needed)
nfiles = numel(files);
wav_npts = zeros(1, nfiles);
wav_data = cell(1, nfiles);

%% Set up report

% Initialize report
filename = 'Neurogram3';
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
pm.PageMargins.Left = '0.1in';
pm.PageMargins.Right = '0.1in';


%% Plot neurograms 

if ismac
	savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/data/AN_neurogram';
else
	savepath = 'C:\Users\jfritzinger\Box\02 - Code\Nat-Timbre\data\AN_neurogram';
end

fontsize = 10;
for ii = 1:nfiles

	% Get single stimulus
	target = extractBefore(files{ii}, '.');
	target_F0 = extractBetween(files{ii}, 'ff.', '.wav');
	savefile = sprintf('%s_F0_%s_AN.mat', target, target_F0{1});
	load(fullfile(savepath, savefile),'AN', 'model_params', 'params')
	F0 = F0s(ii);

	% Label
	h = Heading(2, sprintf("%s F0=%0.0f Hz\n", target, F0));
	b = Border();
	b.BottomStyle = 'single';
	b.BottomColor = 'LightGray';
	b.BottomWidth = '1pt';
	h.Style = [h.Style {Color('Black'), b}, {PageBreakBefore()}];
	append(rpt,h);

	% Load in stimulus 
	stim = params.stimulus(ii,:);

	% Analyze AN_sout
	an_sout = squeeze(AN.an_sout);
	CFs = AN.CFs;
	t = linspace(0, params.dur, size(an_sout,2));
	
	NHN = CFs/F0; % Calculate neural harmonic number (CF/F0) 
	limits = [min(CFs)/F0 30]; %max(CFs)/F0];
	period_lim = 1/F0*10;

	% Analyze spectrogram
	Fs = params.Fs;
    y2 = fft(stim);
    m = abs(y2);
    mdB = 20*log10(m);
    f = (0:length(y2)-1)*Fs/length(y2);
    mdB(mdB<0) = 0;
    f(f>max(CFs)) = [];
    mdB = mdB(1:length(f));
	NHN_stim = f/F0;

	% Plot MASD 
	spec_diff = diff(an_sout, 1);
	spec_abs = abs(spec_diff);
	%MASD = trapz(t, spec_abs, 2);
	MASD = trapz(t, spec_diff, 2);

	% Plot stimulus spectrogram
	figure('Position',[515,275,1082,573])
	tiledlayout(1, 7, 'TileSpacing','none')
	ax1 = nexttile([1,3]);
	plot(NHN_stim, mdB, 'k')
	view(90,90)
	set(gca, 'XDir','reverse')
	title('Stimulus Spectrogram')
	ylabel('Neural Harmonic Number (CF/F0)')
	xlim(limits)
	ylabel('Magnitude')
	set(gca, 'FontSize', fontsize)

	% Plot neurogram 
	ax2 = nexttile([1,3]);
	imagesc(NHN, t, an_sout')
	set(gca, 'YDir', 'normal');
	set(gca, 'XDir', 'reverse');
	view(90,90)
	ylim([0.05 period_lim+0.05])
	title('AN Neurogram')
	ylabel('Time (s)')
	xlim(limits)
	set(gca, 'FontSize', fontsize)

	% Plot MASD
	ax3 = nexttile;
	plot(NHN(1:end-1), MASD, 'k')
	set(gca, 'XDir','reverse')
	view(90,90)
	hold on
	xticks([])
	title('MASD')
	set(gca,'FontSize',fontsize)
	xlim(limits)
	xline(0)
	yline(0)

	linkaxes([ax1 ax2 ax3],'x')

	% Add to PDF
	fig1 = gcf;
	[plt1, images] = addtoSTPDF(images, fig1, num2str(round(F0)));
	append(rpt, plt1);

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

values = [8, 5];
fig.PaperSize = values;
fig.PaperPosition = [0 0 values];
fig.Units = 'inches';
fig.Position(3:4) = values;

% Add the plot to the document
name = sprintf('%s.svg', title);
print(fig, name, '-dsvg');
img = Image(name);
delete(fig) %delete plot figure window
images = [images {img}];

end