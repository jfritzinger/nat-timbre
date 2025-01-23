function MASD = calculateMASD(Ftmp, Ttmp, this_spect, fi, target_F0, iplot)
% Calculates and optionally plots the Mean Absolute Spectral Difference (MASD)
%
% Inputs:
%   Ftmp        - Frequency vector
%   Ttmp        - Time vector
%   this_spect  - Spectrogram matrix
%   fi          - Frequency index limit
%   iplot       - Flag to plot results (1) or not (0)
%
% Output:
%   MASD        - Mean Absolute Spectral Difference vector
%
% The function computes MASD by:
%   1. Truncating input data to the specified frequency index
%   2. Calculating the difference between adjacent spectral frames
%   3. Taking the absolute value of these differences
%   4. Integrating over time using the trapezoidal method
%
% If iplot is set to 1, it also generates a plot of MASD vs. frequency,
% including harmonic lines based on a target fundamental frequency (F0).

% Plot MASD 
f = Ftmp(1:fi);
t = Ttmp;
spec = this_spect(1:fi,:);
spec_diff = diff(spec, 1);
%spec_abs = abs(spec_diff);
%MASD = trapz(t, spec_abs, 2);
MASD = trapz(t, spec_diff, 2);

% figure
% pcolor(1000*Ttmp,Ftmp(1:fi-1),spec_diff);
% shading interp
% ylabel('Frequency (Hz)')
% xlabel('Time (ms)')
% set(gca,'YTickLabel',[]);
% title('Synthesized Spectrogram')
% colormap(gray)
% ax(2).CLim = [max(max(this_spect(1:fi,:)))-20, max(max(this_spect(1:fi,:)))]; % min was 30
% set(gca,'FontSize',fontsize)
% set(gca,'box','off')
% hold on

if iplot == 1
	figure
	plot(f(1:end-1), MASD, 'k')
	set(gca, 'XDir','reverse')
	view(90,90)
	hold on
	num_harms = 25;
	for iharm = 1:num_harms
		xline(target_F0*iharm, ':');
	end
	xticks([])
	title('MASD')
	set(gca,'FontSize',16)
	set(gca,'box','off')
	ylim([-1.5 1.5])
end

end