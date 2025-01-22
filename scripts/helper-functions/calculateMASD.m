function MASD = calculateMASD(Ftmp, Ttmp, this_spect, fi)


% Plot MASD 
f = Ftmp(1:fi);
t = Ttmp;
spec = this_spect(1:fi,:);
spec_diff = diff(spec, 1);
spec_abs = abs(spec_diff);
MASD = trapz(t, spec_abs, 2);

if iplot == 1
	figure
	%ax = axes('Position',[0.83 0.14 0.12 0.80]);
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
	ylim([0 2])
end

end