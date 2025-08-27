function T = getF0NeuronISITable(nat_data, target, index, F0s, type, min_dis)

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

			% Calculate ISI
			edges = 0:min_dis:20;
			for irep = 1:20
				isi = diff(spikes(spikereps==irep));
				h_bass(irep, :) = histcounts(isi, edges);
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
			% Calculate ISI
			edges = 0:min_dis:20;
			for irep = 1:20
				isi = diff(spikes(spikereps==irep));
				h_bass(irep, :) = histcounts(isi, edges);
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
			% Calculate ISI
			edges = 0:min_dis:20;
			for irep = 1:20
				isi = diff(spikes_bass(spikereps_bass==irep));
				h_bass(irep, :) = histcounts(isi, edges);
				isi = diff(spikes_oboe(spikereps_oboe==irep));
				h_oboe(irep, :) = histcounts(isi, edges);
			end
			h_all = [h_all; h_bass; h_oboe];
		end
		T = array2table(h_all);
		T.Response = response;
end
