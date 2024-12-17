% Saves all of the spectra for an instrument into a PDF file for reference
% J. Fritzinger, 7/19/2021

%% Load in data
clear all 

instrument = 'Bass';
fpath = 'C:\Users\jfritzinger\Box\02 - Physiology\07 - Natural Timbre\Cut Waveforms\';
files = dir([fpath instrument '.*']);
%files = dir([fpath '*.wav']);

% Load in tuning sheet
tuning_dir = 'C:\Users\jfritzinger\Box\02 - Physiology\07 - Natural Timbre\Cut Waveforms\';
tuning = readtable([tuning_dir 'Tuning.xlsx']);

% Sort timbre
k = [];
index = [];
for ii = 1:size(files, 1)
    for iii = 1:size(tuning)
        k = strfind(files(ii).name, tuning.Note{iii});
        if ~isnan(k)
            index(ii) = iii;
        end
    end
end
pitch_order = tuning.Frequency(index);
[pitch_order, order] = sort(pitch_order);

files = files(order);
for ii = 1:size(files, 1)
    [y(ii,:), Fs] = audioread([fpath files(ii).name]);
end
t = linspace(0, length(y)/Fs, length(y));



%% Calculate spectra

for ii = 1:size(files, 1)
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
end

%% Set up PDF report
import mlreportgen.dom.*
import mlreportgen.report.*

% Initialize report
images = {}; %hold all plots as images, need to delete when finished
report_name = [extractBefore(files(1).name, '.') '.pdf'];
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

% % Generate page heading
% h = Heading(2, sprintf('%s', extractBefore(files(1).name, '.')));
% b = Border();
% b.BottomStyle = 'single';
% b.BottomColor = 'LightGray';
% b.BottomWidth = '1pt';
% h.Style = [h.Style {Color('Black'), b}, {PageBreakBefore()}];
% append(rpt,h);

for pdf_ind = 1:size(files, 1)
     
    % Adds spectra 
    note = extractBetween(files(pdf_ind).name,'ff.','.wav');
    %string = extractBetween(files(pdf_ind).name,'ff.','.');
    name = [extractBefore(files(pdf_ind).name, '.') ', ' note{1} ', ' num2str(pitch_order(pdf_ind)) 'Hz'];
    [ssplt, images] = pdfStimulusSpectra(f_all(pdf_ind, :), mdB_all(pdf_ind, :),...
        center_of_mass(pdf_ind), pdf_ind, name, images);
    append(rpt, ssplt);
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