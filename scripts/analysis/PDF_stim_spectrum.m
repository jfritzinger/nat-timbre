% Saves all of the spectra for an instrument into a PDF file for reference
% J. Fritzinger, 7/19/2021

%% Load in data

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

% Load in stimulus
for ii = 1:size(files, 2)
    [y(ii,:), Fs] = audioread(fullfile(fpath,'waveforms', files{ii}));
end
t = linspace(0, length(y)/Fs, length(y));

%% Set up PDF report
import mlreportgen.dom.*
import mlreportgen.report.*

% Initialize report
images = {}; %hold all plots as images, need to delete when finished
report_name = 'Bassoon.pdf';
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.1in';
pm.PageMargins.Header = '0.1in';
pm.PageMargins.Bottom = '0.1in';
pm.PageMargins.Footer = '0.1in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';

%% Plot data in PDF format

for ii = 25 %1:size(files, 1)

	y2 = fft(y(ii, :));
    m = abs(y2);
    mdB = 20*log10(m);
    f = (0:length(y2)-1)*Fs/length(y2);
    mdB(mdB<0) = 0;
    f(f>Fs/2) = [];
    mdB = mdB(1:length(f));
    log_f = log10(f(:,2:end));
    center_of_mass(ii) = 10.^((sum(log_f.*mdB(:,2:end), 2))./sum(mdB(:,2:end), 2));
    mdB_all(ii,:) = mdB;
    f_all(ii, :) = f;

	% % Generate page heading
	% h = Heading(2, sprintf('%s', extractBefore(files(1).name, '.')));
	% b = Border();
	% b.BottomStyle = 'single';
	% b.BottomColor = 'LightGray';
	% b.BottomWidth = '1pt';
	% h.Style = [h.Style {Color('Black'), b}, {PageBreakBefore()}];
	% append(rpt,h);
     
    % Adds spectra 
    note = extractBetween(files{ii},'ff.','.wav');
    name = [extractBefore(files{ii}, '.') ', ' note{1} ', ' num2str(F0s(ii)) 'Hz'];
    [ssplt, images] = pdfStimulusSpectra(f_all(ii, :), mdB_all(ii, :),...
        center_of_mass(ii), ii, name, images);
	append(rpt, ssplt);

	% Spectrogram
	npts = round(0.003 * Fs); % 3ms 
	window = gausswin(npts);
	noverlap = round(npts * 0.9);
	[s, w, t] = spectrogram(y(ii,:), window, noverlap, [],Fs,'yaxis');


	figure
	spec = 20*log10(abs(s));
	pcolor(t, w, spec);
	shading interp
	set(gca, 'yscale', 'log')
	colormap(gray)
	colorbar
	yticks([100 200 500 1000 5000 10000])

	
end

close(rpt);
for i = 1:length(images)
    delete(images{1,i}.Path);
end
rptview(rpt)

%% Functions

function [img, images] = pdfStimulusSpectra(f, mdB, center_of_mass, num, name, images)
% Plots the natural timbre spectra
% J. Fritzinger, updated 7/19/2021
import mlreportgen.dom.*

hold on
semilogx(f,mdB);
%[pks, locs] = findpeaks(mdB, 'MinPeakDistance', dist);
%plot(f(locs),pks,'linewidth', 1.2, 'Color', [150 150 150]/255);
xline(center_of_mass, '--', 'linewidth', 1.5);
xlim([50 14000])
ylim([0 70])
grid on
set(gca, 'XScale', 'log')
legend(['SC = ' num2str(round(center_of_mass)) 'Hz']);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB SPL)')
xticks([100 200 500 1000 2000 5000 10000 20000])
box on
title(['Spectrum of ' name])

fig = gcf;
fig.PaperSize = [7.5 3];
fig.PaperPosition = [0 0 7.5 3];
fig.Units = 'inches';
fig.Position(3:4) = [7.5 3];

% Add the plot to the document
name = ['Spectra_',num2str(num)];
imgtype = '-dsvg';
imgname = [name '.svg'];
print(imgtype, imgname);
img = Image(imgname);
delete(gcf) %delete plot figure window
images = [images {img}];

end