function save_figure_MARC(filename)

[~, ~, savepath, ~] = getPathsNT();
% print(fullfile(savepath, '2025-MARC', [filename '.png']),...
% 	'-dpng', '-r600')
% print('-vector', fullfile(savepath, '2025-MARC', [filename '.tif']),...
% 	'-dtiff', '-r600')
print('-vector', fullfile(savepath, '2025-MARC', [filename '.svg']),...
	'-dsvg', '-r600')

end