function T = getTimbrePopTable(nat_data, type, sesh, num_data, type2)

ind_b = 25:40;
ind_o = [1 3:17];

if strcmp(type, 'Rate')

	% Model including all F0s
	data_mat = NaN(2*20, num_data);
	data_mat2 = NaN(2*20*16, num_data);
	for target = 1:16
		for ii = 1:num_data

			% Arrange data for SVM
			X1 = nat_data(sesh(ii)).bass_raterep(:,ind_b(target));
			X2 = nat_data(sesh(ii)).oboe_raterep(:,ind_o(target));
			X = [X1; X2];
			data_mat(:,ii) = X;
		end
		idx = (1:40) + 40*(target-1);
		data_mat2(idx, :) = data_mat;
	end

	% Put data into table
	T = array2table(data_mat2);
	T.Response = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

elseif strcmp(type, 'Timing')
	% 	if strcmp(type2, 'Data')

	h_all2 = [];
	for ineuron = 1:num_data
		sesh_current = sesh(ineuron);
		h_all = [];
		h = [];
		for target = 1:16

			% Get data
			spikes_bass = nat_data(sesh_current).bass_spikerate{ind_b(target)}; % ms
			spikereps_bass = nat_data(sesh_current).bass_spikerep{ind_b(target)};
			spikes_oboe = nat_data(sesh_current).oboe_spikerate{ind_o(target)};
			spikereps_oboe = nat_data(sesh_current).oboe_spikerep{ind_o(target)};
			if strcmp(type2, 'Model')
				spikes_bass = spikes_bass*1000;
				spikes_oboe = spikes_oboe*1000;
			end

			% Arrange data for SVM
			min_dis = 1;
			edges = 0:min_dis:300;
			t = 0+min_dis/2:min_dis:300-min_dis/2;
			for irep = 1:20
				h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
				h_oboe(irep, :) = histcounts(spikes_oboe(spikereps_oboe==irep), edges);
			end
			h_all = [h_bass; h_oboe];
			h = [h; h_all];
		end
		h_all2 = [h_all2, h];
	end

	% 	elseif strcmp(type2, 'Model')
	%
	% 		h_all2 = [];
	% 		for ineuron = 1:num_data
	% 			sesh_current = sesh(ineuron);
	% 			h_all = [];
	% 			h = [];
	% 			for target = 1:16
	%
	% 				% Get data
	% 				spikes_bass = nat_data(sesh_current).bass_PSTH_all{ind_b(target)}; % ms
	% 				spikes_oboe = nat_data(sesh_current).oboe_PSTH_all{ind_o(target)};
	%
	% 				% Arrange data for SVM
	% 				for irep = 1:20
	% 					h_bass(irep, :) = spikes_bass(irep,:);
	% 					h_oboe(irep, :) = spikes_oboe(irep,:);
	% 				end
	% 				h_all = [h_bass; h_oboe];
	% 				h = [h; h_all];
	% 			end
	% 			h_all2 = [h_all2, h];
	% 		end

	% Put data into table
	T = array2table(h_all2);
	response = repmat([ones(1, 20) repmat(2, 1, 20)]', 16, 1);
	T.Response = response;

end

end