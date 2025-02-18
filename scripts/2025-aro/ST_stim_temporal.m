%% ST_stim_temporal
clear

%% Set up figure

ppi = get(0, 'ScreenPixelsPerInch');

figure('Position',[50,50,400,130]);
colors = {'#000000', '#bdbdbd'};
CF_color = '#3690c0';
fontsize = 20;

%% Create Stimulus & plot 

% Parameters
params.Fs = 100000;
params.F0 = 200; % fundamental frequency 
params.Fc = 1400; % peak harmonic frequency 
params.fpeak_mid = 1400;
params.Delta_F = 200;
params.dur = 0.3; % duration
params.ramp_dur = 0.02; % ramp duration
params.stimdB = 70; % overall level
params.G = 24; % slope, in dB/oct 
%params.fpeaks = [850 1100 1400 1700 2050];
params.fpeaks = [900 1400 2250];
params.mnrep = 1;
params.num_harms = 13;
params.stp_otc = 1;
params.physio = 0;
params.g = 24;
params.spl = 70;
params = generate_ST(params);


%% Plot temporal response 

t = linspace(0, 300, params.Fs*params.dur);
stim = params.stim(1,:);
plot(t, stim, 'k')
xlim([100 130])
xlabel('Time (ms)')
ylabel('Amp.')
set(gca,'FontSize', fontsize)
yticklabels([])


%% Export 

savepath = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Nat-Timbre/figures/2025-aro';
set(gcf, 'Renderer', 'painters')
print('-dsvg', '-vector', fullfile(savepath,'ST_stim_temporal.svg'))
