function T = getF0NeuronTable(nat_data, target, index, F0s, type, min_dis)

% switch type
% 	case 'Data'

switch target
	case 'Bassoon'
		h_all = [];
		for itarget = 1:length(F0s)

			if strcmp(type, 'Model')
				spikes = nat_data(index).bass_spikerate{itarget}*1000; % ms
			else
				spikes = nat_data(index).bass_spikerate{itarget}; % ms
			end
			spikereps = nat_data(index).bass_spikerep{itarget};
			
			edges = 0:min_dis:300;
			t = 0+min_dis/2:min_dis:300-min_dis/2;
			for irep = 1:20
				h_bass(irep, :) = histcounts(spikes(spikereps==irep), edges);
			end
			h_all = [h_all; h_bass];
		end
		T = array2table(h_all);
		response = reshape(repmat(F0s, 1, 20)', 1, []);
		response = response';
		T.Response = response;

	case 'Oboe'
		h_all = [];
		for itarget = 1:length(F0s)
			if strcmp(type, 'Model')
				spikes = nat_data(index).oboe_spikerate{itarget}*1000; % ms
			else
				spikes = nat_data(index).oboe_spikerate{itarget}; % ms
			end

			spikereps = nat_data(index).oboe_spikerep{itarget};
			edges = 0:min_dis:300;
			t = 0+min_dis/2:min_dis:300-min_dis/2;
			for irep = 1:20
				h_bass(irep, :) = histcounts(spikes(spikereps==irep), edges);
			end
			h_all = [h_all; h_bass];
		end
		T = array2table(h_all);
		response = reshape(repmat(F0s, 1, 20)', 1, []);
		response = response';
		T.Response = response;

	case 'Invariant'
		response = reshape(repmat(F0s, 1, 40)', 1, []);
		response = response';
		ind_b = 25:40;
		ind_o = [1 3:17];

		% Model including all F0s
		h_all = [];
		for itarget = 1:length(ind_b)
			if strcmp(type, 'Model')
				spikes_bass = nat_data(index).bass_spikerate{ind_b(itarget)}*1000; % ms
				spikereps_bass = nat_data(index).bass_spikerep{ind_b(itarget)};
				spikes_oboe = nat_data(index).oboe_spikerate{ind_o(itarget)}*1000; % ms
				spikereps_oboe = nat_data(index).oboe_spikerep{ind_o(itarget)};
			else
				spikes_bass = nat_data(index).bass_spikerate{ind_b(itarget)}; % ms
				spikereps_bass = nat_data(index).bass_spikerep{ind_b(itarget)};
				spikes_oboe = nat_data(index).oboe_spikerate{ind_o(itarget)}; % ms
				spikereps_oboe = nat_data(index).oboe_spikerep{ind_o(itarget)};
			end

			% Arrange data for SVM
			edges = 0:min_dis:300;
			t = 0+min_dis/2:min_dis:300-min_dis/2;
			for irep = 1:20
				h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
				h_oboe(irep, :) = histcounts(spikes_oboe(spikereps_oboe==irep), edges);
			end
			h_all = [h_all; h_bass; h_oboe];
		end
		T = array2table(h_all);
		T.Response = response;
end

% 	case 'Model'
% 		switch target
% 			case 'Bassoon'
% 				h_all = [];
% 				for itarget = 1:length(F0s)
% 					spikes = nat_data(index).bass_PSTH_all{itarget}; % ms
% 					h_all = [h_all; spikes];
% 				end
% 				T = array2table(h_all);
% 				response = reshape(repmat(F0s, 1, 20)', 1, []);
% 				T.Response = response';
%
% 			case 'Oboe'
% 				h_all = [];
% 				for itarget = 1:length(F0s)
% 					spikes = nat_data(index).oboe_PSTH_all{itarget}; % ms
% 					h_all = [h_all; spikes];
% 				end
% 				T = array2table(h_all);
% 				response = reshape(repmat(F0s, 1, 20)', 1, []);
% 				T.Response = response';
%
% 			case 'Invariant'
% 				response = reshape(repmat(F0s, 1, 40)', 1, []);
% 				ind_b = 25:40;
% 				ind_o = [1 3:17];
%
% 				h = [];
% 				for target = 1:16
%
% 					% Get data
% 					spikes_bass = nat_data(index).bass_PSTH_all{ind_b(target)}; % ms
% 					spikes_oboe = nat_data(index).oboe_PSTH_all{ind_o(target)};
%
% 					% Arrange data for SVM
% 					h_all = [spikes_bass; spikes_oboe];
% 					h = [h; h_all];
% 				end
%
% 				% Put data into table
% 				T = array2table(h);
% 				T.Response = response';
% 		end
% end

end
