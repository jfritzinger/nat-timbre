function save_figure(filename)

[~, ~, savepath, ~] = getPathsNT();
print(fullfile(savepath, '2025-manuscript', [filename '.png']),...
	'-dpng', '-r600')
% print('-vector', fullfile(savepath, '2025-manuscript', [filename '.tif']),...
% 	'-dtiff', '-r600')
% print('-vector', fullfile(savepath, '2025-manuscript', [filename '.eps']),...
% 	'-depsc', '-r600')


end