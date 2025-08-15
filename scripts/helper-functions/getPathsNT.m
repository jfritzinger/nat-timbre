function [base, datapath, savepath, ppi] = getPathsNT()

ppi = get(0, 'ScreenPixelsPerInch');
set(0,'DefaultAxesFontName', 'Arial')

% Set the 'base' filepath for creating all figures 
if ismac
	base = ['/Users/jfritzinger/Library/CloudStorage/Box-Box/02-projects/' ...
		'nat-timbre/data/'];
	savepath = ['/Users/jfritzinger/Library/CloudStorage/Box-Box/02-projects' ...
		'/nat-timbre/figures/'];
	%base = '/Volumes/DataFiles_JBF/Nat-Timbre/data';
else
	base = ['C:\Users\jfritzinger\Box\02-projects\nat-timbre\' ...
		'data'];
	savepath = ['C:\Users\jfritzinger\Box\02-projects\nat-timbre\' ...
		'figures\2025-manuscript'];
end

% Basic paths for loading data and saving figures 
datapath = fullfile(base, 'Neural_Data');
end