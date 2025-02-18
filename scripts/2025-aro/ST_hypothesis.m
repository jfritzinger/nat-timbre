%% ST_hypothesis
clear 

%% Set up figure
figure('Position',[25,476,620,670])
font_size = 16;

%% Parameters
CF = 1200; 

% Stimulus parameters
params.fpeak_mid = 1200;
params.Delta_F = 200;
params.num_harms = 11;
params.stp_otc = 1;
params.Fs = 100000;
params.mnrep = 1;
params.physio = 0;
params.dur = 0.3;
params.ramp_dur = 0.02;
params.spl = 70;
params.g = 24;
params = generate_ST(params);

% Model parameters
model_params.type = 'SFIE';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.BMF = 100;
model_params.CF_range = [125 10000];
model_params.num_CFs = 100;
model_params.CFs = logspace(log10(125), log10(10000), 100);
model_params.nAN_fibers_per_CF = 5;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 10; % how many times to run the AN model
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw 
% implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - 
% this is the 'noise' associated with spont. activity of AN fibers - 
% see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to 
% omit 1st 50 ms, use 0.050
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium 
% SR, 3 = high SR)
model_params.Fs = 100000;


%% Model 

%AN_HSR = modelAN(params, model_params); % HSR for IC input
%SFIE = wrapperIC(AN_HSR.an_sout, params, model_params); % SFIE output
%save('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Aim 2 - 
% Timbre/Data/Intro_ModelResponse.mat', 'AN_HSR', 'SFIE')

addpath '/Users/jfritzinger/Projects/synth-timbre/scripts/helper-functions'
[base, datapath, ~] = getPaths();
filename = 'Intro_ModelResponse.mat';
load(fullfile(datapath, filename), 'AN_HSR', 'SFIE')
avAN = AN_HSR.average_AN_sout;
avBE = SFIE.average_ic_sout_BE;
avBS = SFIE.average_ic_sout_BS;
CFs = AN_HSR.CFs;

%% Plot temporal 

% Figure parameters
CF_list = [875 1200 1550];
start_time = 0.1; % s - for plotting
start  = start_time * params.Fs;
stop_time = 0.115; % s
stop = stop_time * params.Fs;
b_LP = fir1(5000,100/(params.Fs/2),'low'); % LP at 100 Hz, for envelope plots

iplot = 0; % step thru CFs within each panel

for CF_plot = CF_list  % for TIN response to 1000 Hz tone
    iplot = iplot + 1;

	% AN model responses
	vihc = model_IHC(params.stim,CF_plot,model_params.nrep,1/params.Fs,...
		params.dur*1.2,model_params.cohc,model_params.cihc,model_params.species);
	[an_sout,~,~] = model_Synapse(vihc,CF_plot,model_params.nrep,1/params.Fs,...
		model_params.fiberType,model_params.noiseType,model_params.implnt); 
	an_sout_LP = conv(an_sout,b_LP,'Same'); 
	env_an_sout = envelope(an_sout_LP);
    
    % Plot AN Temporal
    t = (1:length(an_sout))/params.Fs; % time vector for plots
    
    h(iplot) = subplot(5, 5, 11+iplot);
    hold on
    set(gca,'fontsize',14, 'xtick',[100 125 150], 'YDir', 'reverse')
    
    % Duplicate the AN responses for Fluctuation profile plot
    plot(an_sout(start:stop), t(start:stop)*1e3,'Color', '#44AA99',...
		'linewidth',1.5)
    [peaks_vals,peak_locations] = findpeaks(an_sout(start:stop),params.Fs,...
		'MinPeakProminence',100);
    plot(peaks_vals, (peak_locations + start_time)*1e3,'Color', '#117733',...
		'linewidth',3) % add envelope
    
    ylim([start_time*1e3, stop_time*1e3])
    xlim([0 1000]); % freq
    
    if iplot == 1
        ylabel('Time (ms)')
        yticks([100 115])
        yticklabels([100 115])
    elseif iplot == 2
        yticks([100 115])
        yticklabels([])
        xlabel('Rate (sp/s)')
    else
        yticks([100 115])
        yticklabels([])
    end
    set(gca,'fontsize',font_size)
    xticks([0 1000])
    grid on
	xticklabels([])
    
end

%% Plot stimulus 

plot_range = [300 6000];
ticks = [200 500 1000 2000 5000];

% Stimulus Creation
h(4) = subplot(5,5,1:5);
hold on
harmonics = params.Delta_F:params.Delta_F:10000;
num_harmonics = length(harmonics);
npts = params.dur * params.Fs; % # pts in stimulus
t = (0:(npts-1))/params.Fs; % time vector
component_scales_linear = 10.^(-1*abs(log2(harmonics/params.fpeak_mid)*...
	params.g)/20);
stimulus = zeros(1,npts);
for iharm = 1:num_harmonics
    comp_freq = harmonics(iharm);
    component = component_scales_linear(iharm) * sin(2*pi*comp_freq*t);
    stimulus = stimulus + component;          %Add component to interval
end
Level_scale = 20e-6*10.^(params.spl/20) * (1/rms(stimulus));
component_scales_linear = Level_scale * component_scales_linear;

% Now, make the stimulus for this_fpeak
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
xline(CF_list(1), ':','Color', [0.4 0.4 0.4], 'linewidth', 2);
xline(CF_list(3), ':','Color', [0.4 0.4 0.4], 'linewidth', 2);

for iharm = 1:num_harmonics
    stim = component_scales_linear(iharm) * sin(2*pi*harmonics(iharm)*t);
    
    % Plot each stimulus
    y = fft(stim);
    mdB = 20*log10(abs(y));
    level(iharm) = findpeaks(mdB(1:length(mdB)/2), 'MinPeakProminence',200);
    stem(harmonics(iharm), level(iharm), 'Marker', 'none', 'LineWidth', ...
        3, 'Color', 'k');
end

% Plot envelope of stimulus
plot(harmonics, level, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':');
set(gca,'fontsize',font_size)
ylabel('Level (dB SPL)')
ylim([0 70])
grid on
set(gca, 'XScale', 'log')
hold off
xlim(plot_range)
xticks(ticks)
xlabel('Freq. (Hz)')


%% AN Plot
h(5) = subplot(5, 5, 6:10);
hold on
plot(CFs, avAN, 'linewidth', 2, 'color', '#117733');
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
xlim(plot_range)
ylim([0 300])
xticks(ticks)
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
xline(CF_list(1), ':','Color', [0.4 0.4 0.4], 'linewidth', 2);
xline(CF_list(3), ':','Color', [0.4 0.4 0.4], 'linewidth', 2);
xlabel('CF (Hz)')
ylabel('Avg. Rate (sp/s)')

%% IC BE Plot
h(6) = subplot(5, 5, 17:20);
hold on
plot(CFs, avBE, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
xlim(plot_range)
ylim([0 50])
xticks(ticks)
yticks([0 30 60])
xticklabels([])
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
xline(CF_list(1), ':','Color', [0.4 0.4 0.4], 'linewidth', 2);
xline(CF_list(3), ':','Color', [0.4 0.4 0.4], 'linewidth', 2);

%% IC BS Plot
h(7) = subplot(5, 5, 22:25);
hold on
plot(CFs, avBS, 'linewidth', 2, 'color', [0, 0.4470, 0.7410]);
set(gca,'fontsize',font_size)
grid on
set(gca, 'XScale', 'log')
xlim(plot_range)
ylim([0 50])
xticks(ticks)
yticks([0 30 60])
xlabel('CF (Hz)')
xline(CF, '--', 'Color', [0.4 0.4 0.4], 'linewidth', 2); % Add CF line
xline(CF_list(1), ':','Color', [0.4 0.4 0.4], 'linewidth', 2);
xline(CF_list(3), ':','Color', [0.4 0.4 0.4], 'linewidth', 2);

%% Plot MTFs
% clear params
%
% % Stimulus Parameters and Generation
% params.type = 'typMTFN';
% params.ramp_dur = 0.05;
% params.noise_state = 0;
% params.noise_band = [100, 10000];
% params.dur = 1; %1; % s
% params.reptim = 1.5;
% params.fms = [2, 600, 3]; % fm_lo, fm_hi, steps per octave
% params.mdepths = [0,0,1];
% params.binmode = 2;
% params.No = 30;
% params.spl = 30;
% params.all_mdepths = 0;
% params.Fs = 100000;
% params.nrep = 1;
% params.mnrep = 5;
% params.raised_sine = 1;
% params = generate_MTF(params);
% params.num_stim = size(params.stim, 1);
% 
% % Model Parameters 
% model_params.type = 'SFIE';
% model_params.range = 2; % 1 = population model, 2 = single cell model
% model_params.species = 1; % 1 = cat, 2 = human
% model_params.BMF = 100;
% model_params.CF_range = 1200;
% model_params.num_CFs = 1;
% model_params.CFs = 1200;
% model_params.nAN_fibers_per_CF = 10;
% model_params.cohc = 1; % (0-1 where 1 is normal)
% model_params.cihc = 1; % (0-1 where 1 is normal)
% model_params.nrep = 1; % how many times to run the AN model
% model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw 
% implementation(See Zilany etal., 2009)
% model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this 
% is the 'noise' associated with spont. activity of AN fibers - see 
% Zilany et al., 2009. "0" lets you "freeze" it.
% model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
% model_params.onsetWin = 0.020; % exclusion of onset response, e.g. 
% to omit 1st 50 ms, use 0.050
% model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium 
% SR, 3 = high SR)
% model_params.Fs = 100000;
% 
% % Run models 
% AN_HSR = modelAN(params, model_params); % HSR for IC input
% SFIE = wrapperIC(AN_HSR.an_sout, params, model_params); % SFIE output
% save('/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Aim 2 - 
% Timbre/Data/Intro_MTFModelResponse.mat', 'AN_HSR', 'SFIE', 'params')

filename = 'Intro_MTFModelResponse.mat';
load(fullfile(datapath, filename), 'AN_HSR', 'SFIE', 'params')

for ineuron = 1:2

	if ineuron == 1
		h(8) = subplot(5, 5, 16);
		model = SFIE.average_ic_sout_BE;
	else
		h(9) = subplot(5, 5, 21);
		model = SFIE.average_ic_sout_BS;
	end
	[~, avBE, stdBE] = plotMTF(params, model, 0);
	hold on
	line([1 params.all_fms(end)], [1 1]*avBE(1),'Color',[0.4 0.4 0.4], ...
		'LineWidth', 2);
	plot(params.all_fms,smooth(avBE),'-k', 'LineWidth', 2) % smoothed MTF
	hold off
	xtick = [1 2 5 10 20 50 100 200 500];
	xlim(xtick([1 end]))
	set(gca,'XTick',xtick,'XScale', 'log')
	grid on
	set(gca,'FontSize',font_size)
	xticklabels([])
	yticklabels([])
	if ineuron == 1
		ylim([22 32])
		ylabel('Avg. Rate (sp/s)                         ')
		hLegend = legend('Unmod.', 'Location','northwest', 'EdgeColor',...
			'none');
		hLegend.ItemTokenSize = [12,8];
		title('MTF')
	else
		ylim([15 34])
	end
end


%% Move positions

bottom = fliplr(linspace(0.07, 0.83, 5));
bottom(3) = 0.43;
bottom(4) = 0.23;
height = 0.135;

set(h(4), 'Position',[0.3 bottom(1) 0.655 height]);
set(h(5), 'Position',[0.3 bottom(2) 0.655 height]);

set(h(1), 'Position',[0.38 bottom(3) 0.13 height]);
set(h(2), 'Position',[0.55 bottom(3) 0.13 height]);
set(h(3), 'Position',[0.72 bottom(3) 0.13 height]);

set(h(8), 'Position',[0.14 bottom(4) 0.123 height]);
set(h(6), 'Position',[0.3 bottom(4) 0.655 height]);

set(h(9), 'Position',[0.14 bottom(5) 0.123 height]);
set(h(7), 'Position',[0.3 bottom(5) 0.655 height]);

%% Labels / Annotations / Positions 
titlesize = 22; 

annotation('textbox',[0.3 0.84 0.15 0.1], 'String','Stimulus',...
	'FontSize',titlesize,'EdgeColor','none', 'Rotation',90, ...
	'FitBoxToText','off', 'FontWeight','bold');
annotation('textbox', [0.3 0.66 0.15 0.1],'String','AN Rate',...
	'FontSize',titlesize,'FitBoxToText','off',...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
annotation('textbox',[0.3 0.43 0.2 0.1],'String','AN Temporal',...
	'FontSize',titlesize,'FitBoxToText','off',...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');
annotation('textbox',[0.17 0.11 0.15 0.1],'String','IC BS',...
	'FontSize',titlesize, 'FontWeight','bold',...
	'FitBoxToText','off','EdgeColor','none', 'Rotation',90);
annotation('textbox',[0.17 0.27 0.15 0.1],'String','IC BE',...
	'FontSize',titlesize,'FitBoxToText','off',...
	'EdgeColor','none', 'Rotation',90, 'FontWeight','bold');

% Create lines
annotation('line',[0.442 0.535],[0.43 0.361],'Color',[0.4 0.4 0.4],...
	'LineWidth',2,'LineStyle',':');
annotation('line',[0.442 0.53],[0.566 0.637],'Color',[0.4 0.4 0.4],...
	'LineWidth',2,'LineStyle',':');

annotation('line',[0.765 0.660],[0.43 0.361],'Color',[0.4 0.4 0.4],...
	'LineWidth',2,'LineStyle',':');
annotation('line',[0.765 0.660],[0.566 0.637],'Color',[0.4 0.4 0.4],...
	'LineWidth',2,'LineStyle',':');

annotation('line',[0.603 0.603],[0.566 0.637],'Color',[0.4 0.4 0.4],...
	'LineWidth',2,'LineStyle','--');
annotation('line',[0.603 0.603],[0.43 0.361],'Color',[0.4 0.4 0.4],...
	'LineWidth',2,'LineStyle','--');


%% Export 

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/2025-aro';
set(gcf, 'Renderer', 'painters')
print('-dsvg', '-vector', fullfile(savepath,'ST_hypothesis.svg'))


