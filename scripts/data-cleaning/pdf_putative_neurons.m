%% Create PDF of all putative neurons for timbre
% J. Fritzinger, 8/27/24
clear 
import mlreportgen.dom.*
import mlreportgen.report.*

%% Load in spreadsheet 

[base, datapath, savepath, ppi] = getPaths();
sheetname = 'Data_Table.xlsx';
sessions = readtable(sheetname, 'PreserveVariableNames',true);
num_neurons = size(sessions, 1)-1;

%% Set up PDF 

% Identify current machine, which changes paths
[userid, base_dir, ~, report_path, data_path] = findPaths();
filename = 'PutativeTimbreResponses';

% Initialize report
images = {}; %hold all plots as images, need to delete when finished
datetime.setDefaultFormats('default','yyyy-MM-dd_hhmmss')
report_name = sprintf('%s%s_%s.pdf', report_path, datetime, filename);
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.01in';
pm.PageMargins.Header = '0.01in';
pm.PageMargins.Bottom = '0.01in';
pm.PageMargins.Footer = '0.01in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';

%% Plot all sessions 

% Only one MTF type & sort by CF 
% MTF_list = strcmp(sessions.MTF, 'BE');
% index = find(MTF_list);
% CF_list = sessions.CF(MTF_list);
% [~, order] = sort(CF_list);
% num_sessions = length(CF_list);


CF_list = sessions.CF(MTF_list);
[~, order] = sort(CF_list);
num_sessions = length(CF_list);

timerVal = tic;
for isession = 1:num_sessions

    % Label the session
    % 	putative = sessions.Putative_Units{ineuron}; % Grabs multi/good/excellent excel file
    % 	CF = sessions.CF(ineuron);
    % 	MTF_shape = sessions.MTF{ineuron};
    ineuron = index(order(isession));
    putative = sessions.Putative_Units{ineuron}; % Grabs multi/good/excellent excel file
    CF = sessions.CF(ineuron);
    MTF_shape = sessions.MTF{ineuron};

    % Set message
    fprintf('Creating plots... %s, CF = %0.0fHz, %s\n', putative, CF, MTF_shape);

    % Set heading
	h = Heading(2, sprintf('%s, CF = %0.0fHz, %s', putative, CF, MTF_shape));
	b = Border();
	b.BottomStyle = 'single';
	b.BottomColor = 'LightGray';
	b.BottomWidth = '1pt';
	h.Style = [h.Style {Color('Black'), b}, {PageBreakBefore()}];
	append(rpt,h);

	% Load in data 
	filename = sprintf('%s.mat', putative);
	load(fullfile(datapath,'neural_data', filename));

	% Plot BIN/RM/MTF
	BINplt = [];
	rmplt = [];
	mtfplt = [];
	if ~isempty(sessions.char_spl{ineuron})
		[fig, data_BIN] = plotPhysBIN(data{1,2}.cluster, data{1,2}, data{1,2}.stim);
		[BINplt, images] = addToPutativePDF(data{1,2}.cluster, data{1,2}, putative, images, fig);
	end
	if ~isempty(sessions.("type=RM"){ineuron})
		[fig, data_RM] = plotPhysRM(data{2,2}.cluster, data{2,2},data{2,2}.stims, CF);
		[rmplt, images] = addToPutativePDF(data{2,2}.cluster, data{2,2}, putative, images, fig);
		rmplt.Height = '2in';
		rmplt.Width = '3.5in';
	end
	if ~isempty(sessions.typMTFN{ineuron})
		[~,~,~, fig, data_MTF] = plotPhysMTF([], data{3,2},[]);
		[mtfplt, images] = addToPutativePDF(data{3,2}.cluster, data{3,2}, putative, images, fig);
	end
	addMultipleChar(rpt, BINplt, mtfplt, rmplt);

	% Plot RVF
	if ~isempty(sessions.RVF{ineuron})
		fig = plotPhysRVF([], data{5,1},[]);
		[rvfplt, images] = addToPutativePDF(data{5,1}.cluster, data{5,1}, putative, images, fig);
		append(rpt, rvfplt);
	end

	% Plot STRFs
	if ~isempty(sessions.STRF{ineuron})
		[~, fig, data_STRF] = plotPhysSTRF([],data{4,2}, []);
		data{4,2}.plot_type = 'STRF';
		[strfplt, images] = addToPutativePDF(data{4,2}.cluster, data{4,2},putative, images, fig);
		append(rpt, strfplt);
	end

	% Plot Natural Timbre
	data_NT = data(13:end, 2);
	data_NT = data_NT(~cellfun(@isempty, data_NT));
	if ~isempty(data_NT)
		[~, fig, data_NT] = plotPhysNT([], data_NT,[], 0, data_NT);
		data_NT{1}.plot_type = 'NatTimbre';
		[ntplt, images] = addToPutativePDF(data_NT{1,1}.cluster, data_NT, putative, images, fig);
		append(rpt, ntplt);
	end
end


%% Save PDF 

close(rpt);
for i = 1:length(images)
	delete(images{1,i}.Path);
end
rptview(rpt)

elapsedTime = toc(timerVal)/60;
disp(['This took ' num2str(elapsedTime) ' minutes'])

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function addMultipleChar(rpt, BINplt, mtfplt, rmplt)
import mlreportgen.dom.*
import mlreportgen.report.*

if ~isempty(BINplt) && ~isempty(mtfplt) && ~isempty(rmplt)
	Bin_table = Table({BINplt,rmplt,mtfplt});
	Bin_table.Style = {Width('100%'), ResizeToFitContents(false)};
	Bin_table.BorderColor = 'White';
	append(rpt, Bin_table);
elseif ~isempty(BINplt) && ~isempty(mtfplt) 
	Bin_ITD_table = Table({BINplt,mtfplt});
	Bin_ITD_table.Style = {Width('100%'), ResizeToFitContents(false)};
	Bin_ITD_table.BorderColor = 'White';
	append(rpt, Bin_ITD_table);
elseif ~isempty(BINplt) && ~isempty(rmplt) 
	Bin_ILD_table = Table({BINplt,rmplt});
	Bin_ILD_table.Style = {Width('100%'), ResizeToFitContents(false)};
	Bin_ILD_table.BorderColor = 'White';
	append(rpt, Bin_ILD_table);
elseif isempty(mtfplt) && ~isempty(rmplt)
	ILD_ITD_table = Table({rmplt,mtfplt});
	ILD_ITD_table.Style = {Width('100%'), ResizeToFitContents(false)};
	ILD_ITD_table.BorderColor = 'White';
	append(rpt, ILD_ITD_table);
elseif ~isempty(BINplt)
	append(rpt, BINplt);
elseif ~isempty(rmplt)
	append(rpt, rmplt);
elseif ~isempty(mtfplt) 
	append(rpt, mtfplt);
end
end
