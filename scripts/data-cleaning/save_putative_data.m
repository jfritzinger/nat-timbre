%% Fig0_SavingData.m
% J. Fritzinger, updated 10/23/23
%
% Saves all data for WB-TIN paper
clear

%% Load in PDF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spreadsheet_name = 'PutativeTable2.xlsx';
sessions = readtable(spreadsheet_name, 'PreserveVariableNames',true);
num_data = size(sessions, 1);

if ismac
	path = '/Users/jfritzinger/Library/CloudStorage/Box-Box/02 - Code/Synth-Timbre/data/manuscript/neural_data';
else
	path = 'C:\Users\jfritzinger\Box\02 - Code\Synth-Timbre\data\manuscript\neural_data';
end

for isesh = 1:num_data

	% Name of file to save
	filename = sessions.Putative_Units{isesh};

	%% Load in each session %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Create struct for each dataset
	dataset_types = {'Binaural Response', 'Response Map', 'MTFN', 'SCHR', 'STRF',...
		'RVF', 'Oboe', 'Bassoon', 'Other', 'RVF'};
	num_dataset_types = length(dataset_types);
	data = cell(num_dataset_types, 2);

	for ind = 10:34 

		% Get the path to each session
		full_name1 = sessions{isesh,ind};
		full_name = full_name1{1};

        if isempty(full_name)
            full_name1 = sessions{isesh,10};
    		full_name = full_name1{1};
            rabbit = str2num(full_name(3:4));
    		session = str2num(full_name(6:8));
    		tetrode = str2num(full_name(12));
    		neuron = str2num(full_name(15));
            [userid, base_dir, ~, report_path, data_path] = findPaths();
            session_name = sprintf('R%03dS%03d', rabbit, session);
            rab_num = num2str(rabbit);
            if ismac
                session_dir_name = fullfile(base_dir, ['R0' num2str(rab_num)]);
            else
                session_dir_name = base_dir{contains(base_dir, rab_num)};
            end
            session_dir = fullfile(session_dir_name, session_name);

            % Post_process
            [clusters, params, stims] = loadPhysiologySession(session_dir, session_name, userid);

            % Get data for tetrode/neuron of interest
            cluster = clusters([clusters.tetrode] == tetrode); % Select tetrode
            cluster = cluster([cluster.neuron] == neuron); % Select neuron

            % Check if spreadsheet matches data
            [bin, contra, F0_100, F0_200, level_43, level_63, level_73, level_83, even, ...
        		oboe, bassoon, other] = finddata(params);
            switch ind
                case 17
                    ST = cellfun(@(p) strcmp(p.type,'STRF')&&p.binmode==2,params);
                case 18
                    ST = cellfun(@(p) strcmp(p.type,'STRF')&&p.binmode==1,params);
                case 19
                    ST = cellfun(@(p)strcmp(p.type,'SCHR'),params);
                case 20
                    ST = cellfun(@(p)strcmp(p.type,'RVF'),params);
                case 21
            	    ST = any(oboe);
                case 22
            	    ST = any(bassoon);
                case 23
            	    ST = any(other);
            end
            if ind <17
                ST = 0;
            end
            if sum(ST)>0
               %error([full_name ': Forgot data for index ' num2str(ind)])
            end
            continue
        end

		rabbit = str2num(full_name(3:4));
		session = str2num(full_name(6:8));
		tetrode = str2num(full_name(12));
		neuron = str2num(full_name(15));

		[userid, base_dir, ~, report_path, data_path] = findPaths();
		session_name = sprintf('R%03dS%03d', rabbit, session);
		rab_num = num2str(rabbit);
		if ismac
			session_dir_name = fullfile(base_dir, ['R0' num2str(rab_num)]);
		else
			session_dir_name = base_dir{contains(base_dir, rab_num)};
		end
		session_dir = fullfile(session_dir_name, session_name);

		% Post_process
		[clusters, params, stims] = loadPhysiologySession(session_dir, session_name, userid);

		% Get data for tetrode/neuron of interest
		cluster = clusters([clusters.tetrode] == tetrode); % Select tetrode
		cluster = cluster([cluster.neuron] == neuron); % Select neuron

		%% Organize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% Find sessions of interest
		putative = sessions.Putative_Units{isesh};
		for ids = 1:length(params)
			params{ids}.putative = putative;
			params{ids}.tetrode = tetrode;
			params{ids}.neuron = neuron;
		end

		% Check if this session has an spl shift
		spl_shift = checkIfSPLShift(rabbit, session);

		switch ind
			case 10 %%%%%% Binaural noise response
				has_bin = cellfun(@(p)strcmp(p.type,'char_spl'), params);
				data_binnoise = params{has_bin};
				data_binnoise.spls = [data_binnoise.spls(1:2)+spl_shift data_binnoise.spls(3)];
				for ilist = 1:length(data_binnoise.list)
					data_binnoise.list(ilist).spl = data_binnoise.list(ilist).spl+spl_shift;
				end
				data_binnoise.stim = get_stims(data_binnoise, stims);
				data_binnoise.cluster = cluster;
				data{1,2} = data_binnoise;

			case 11 % ITD
				% Not saving
			case 12 % ILD
				% Not saving
			case 13 %%%%%% RM binaural
				has_rm = cellfun(@(p)strcmp(p.type,'type=RM')&&p.binmode==2, params);
				if any(has_rm)
					data_rm = params{has_rm};
					data_rm.spls = [data_rm.spls(1:2)+spl_shift data_rm.spls(3)];
					data_rm.all_spls = data_rm.all_spls+spl_shift;
					for ilist = 1:length(data_rm.list)
						data_rm.list(ilist).spl = data_rm.list(ilist).spl+spl_shift;
					end
					data_rm.stims = get_stims(data_rm, stims);
					data_rm.cluster = cluster;
					data{2,2} = data_rm;
				end

			case 14 %%%%%% RM contra
				has_rm = cellfun(@(p)strcmp(p.type,'type=RM')&&p.binmode==1, params);
				if any(has_rm)
					data_rm = params{has_rm};
					data_rm.spls = [data_rm.spls(1:2)+spl_shift data_rm.spls(3)];
					data_rm.all_spls = data_rm.all_spls+spl_shift;
					for ilist = 1:length(data_rm.list)
						data_rm.list(ilist).spl = data_rm.list(ilist).spl+spl_shift;
					end
					data_rm.stims = get_stims(data_rm, stims);
					data_rm.cluster = cluster;
					data{2,1} = data_rm;
				end

			case 15 %%%%%% MTFN binaural
				has_mtf = cellfun(@(p)strcmp(p.type,'typMTFN')&&(p.all_mdepths(1) == 0)&&...
					p.binmode==2, params);
				if any(has_mtf)
					data_mtf = params{has_mtf};
					data_mtf.spl = data_mtf.spl+spl_shift;
					data_mtf.stims = get_stims(data_mtf, stims);
					data_mtf.cluster = cluster;
					data{3,2} = data_mtf;
				end

			case 16 %%%%%% MTFN contra
				has_mtf = cellfun(@(p)strcmp(p.type,'typMTFN')&&(p.all_mdepths(1) == 0)&&...
					p.binmode==1, params);
				if any(has_mtf)
					data_mtf = params{has_mtf};
					data_mtf.spl = data_mtf.spl+spl_shift;
					data_mtf.stims = get_stims(data_mtf, stims);
					data_mtf.cluster = cluster;
					data{3,1} = data_mtf;
				end

			case 17 %%%%%% STRF binaural
				has_strf = cellfun(@(p) strcmp(p.type,'STRF')&&p.binmode==2,params);
				if any(has_strf)
					data_strf = params{has_strf};
					data_strf.spl = data_strf.spl + spl_shift;
					data_strf.stims = get_stims(data_strf, stims);
					data_strf.cluster = cluster;
					data{4,2} = data_strf;
				end

			case 18 %%%%%% STRF contra
				has_strf = cellfun(@(p) strcmp(p.type,'STRF')&&p.binmode==1,params);
				if any(has_strf)
					data_strf = params{has_strf};
					data_strf.spl = data_strf.spl + spl_shift;
					data_strf.stims = get_stims(data_strf, stims);
					data_strf.cluster = cluster;
					data{4,1} = data_strf;
				end

			case 19 %%%%%% SCHR
				has_schr = cellfun(@(p)strcmp(p.type,'SCHR'),params);
				if any(has_schr)
					data_schr = params{has_schr};
					data_schr.stimdB = data_schr.stimdB+spl_shift;
					data_schr.stim = get_stims(data_schr, stims);
					data_schr.cluster = cluster;
					data{5,2} = data_schr;
				end
			case 20 %%%%%% RVF
				has_rvf = cellfun(@(p)strcmp(p.type,'RVF'),params);
				if any(has_rvf)
					data_rvf = params{has_rvf};
					data_rvf.stim = get_stims(data_rvf, stims);
					data_rvf.cluster = cluster;
					data{5,1} = data_rvf;
				end
			case 21 %%%%%% Oboe
				has_oboe = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
					(isfield(p, 'target') && strcmp(p.target, 'Oboe')||...
					contains(p.list(1).wav_file, 'Oboe')), params);
				if any(has_oboe)
					data_oboe = params{has_oboe};
					data_oboe.spl = data_oboe.signal_spls+spl_shift;
					data_oboe.stims = get_stims(data_oboe, stims);
					data_oboe.cluster = cluster;
					data{6,2} = data_oboe;
				end

			case 22 %%%%%% Bassoon
				has_bassoon = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
					(isfield(p, 'target') && strcmp(p.target, 'Bassoon')||...
					contains(p.list(1).wav_file, 'Bassoon')), params);
				if any(has_bassoon)
					data_bassoon = params{has_bassoon};
					data_bassoon.spl = data_bassoon.signal_spls+spl_shift;
					data_bassoon.stims = get_stims(data_bassoon, stims);
					data_bassoon.cluster = cluster;
					data{7,2} = data_bassoon;
				end

			case 23 %%%%%% Other natural timbre
				has_nt = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre'), params);
				index = find(has_nt);

				if any(has_nt)
					for ii = 1:length(index)
						jj= index(ii);
						idata = 8+ii;

						data_nt = params{jj};
						data_nt.spl = data_nt.signal_spls+spl_shift;
						data_nt.stims = get_stims(data_nt, stims);
						data_nt.cluster = cluster;
						data{idata,2} = data_nt;
					end
				end
		end
	end

	% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	save(fullfile(path, filename), 'data')

end


%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to get the stimulus properties for each DSID
function [stim] = get_stims(data, stims)
ds = data.dsid;
stim = stims(ds);
end

function [oboe, bassoon, other] = finddata(population)

	% Find oboe 
	oboe = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
		(isfield(p, 'target') && strcmp(p.target, 'Oboe')||...
		contains(p.list(1).wav_file, 'Oboe')), population);

	% Find bassoon
	bassoon = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre') && ...
		(isfield(p, 'target') && strcmp(p.target, 'Bassoon')||...
		contains(p.list(1).wav_file, 'Bassoon')), population);

	% Find other
	other = cellfun(@(p) ~isempty(p) && strcmp(p.type, 'Natural_Timbre'), population)...
        & ~bassoon & ~oboe;

end