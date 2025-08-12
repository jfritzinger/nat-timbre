function T = getAllPopTable(nat_data, type)


switch type % rate or timing
	case 'Rate'
		%%
		% Get stimulus
		F0s_b = getF0s('Bassoon');
		F0s_o = getF0s('Oboe');

		% Get into 'Response'
		response_b = cell(75,1);
		for ii = 1:length(F0s_b)
			response_b{ii} = ['B_' num2str(round(F0s_b(ii)))];
		end
		for ii = 1:length(F0s_o)
			response_b{ii+length(F0s_b)} = ['O_' num2str(round(F0s_o(ii)))];
		end
		response = reshape(repmat(response_b, 1, 20)', 1, []);
		response = response';

		% Get data
		[sesh, num_data] = getTimbreSessions(nat_data);

		data_mat1 = NaN(length(F0s_b)*20, num_data);
		for ii = 1:num_data
			X1 = nat_data(sesh(ii)).bass_raterep';
			X2 = reshape(X1', [], 1);
			data_mat1(:,ii) = X2;
		end

		data_mat2 = NaN(length(F0s_o)*20, num_data);
		for ii = 1:num_data
			X1 = nat_data(sesh(ii)).oboe_raterep';
			X2 = reshape(X1', [], 1);
			data_mat2(:,ii) = X2;
		end

		data_mat = [data_mat1; data_mat2];
		T = array2table(data_mat);
		T.Response = response;


	case 'Timing'
		%%
		% Get bassoon stimulus
		F0s_b = getF0s('Bassoon');
		F0s_o = getF0s('Oboe');

		% Get into 'Response'
		response_b = cell(75,1);
		for ii = 1:length(F0s_b)
			response_b{ii} = ['B_' num2str(round(F0s_b(ii)))];
		end
		for ii = 1:length(F0s_o)
			response_b{ii+length(F0s_b)} = ['O_' num2str(round(F0s_o(ii)))];
		end
		response = reshape(repmat(response_b, 1, 20)', 1, []);
		response = response';

		num_data = size(nat_data, 2);
		h_all2 = [];
		h_all = [];
		for ineuron = 1:num_data
			h_all_bass = [];
			for itarget = 1:length(F0s_b)

				% Arrange bassoon data for SVM
				spikes = nat_data(ineuron).bass_spikerate{itarget}/1000; % ms
				spikereps = nat_data(ineuron).bass_spikerep{itarget};
				min_dis = 0.25;
				edges = 0:min_dis:300;
				for irep = 1:20
					h_bass(irep, :) = histcounts(spikes(spikereps==irep), edges);
				end
				h_all_bass = [h_all_bass; h_bass];
			end

			h_all_oboe = [];
			for itarget = 1:length(F0s_o)
				spikes = nat_data(ineuron).oboe_spikerate{itarget}/1000; % ms
				spikereps = nat_data(ineuron).oboe_spikerep{itarget};
				for irep = 1:20
					h_oboe(irep, :) = histcounts(spikes(spikereps==irep), edges);
				end
				h_all_oboe = [h_all_oboe; h_oboe];
			end
			h_all = [h_all_bass; h_all_oboe];
			h_all2 = [h_all2, h_all];
		end
		T = array2table(h_all2);
		T.Response = response;

end