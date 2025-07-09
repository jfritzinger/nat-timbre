function T = getTimbreNeuronTable(nat_data, index, type)
ind_b = 25:40;
ind_o = [1 3:17];
% switch type % data or model
% 	case 'Data'

h = [];
for target = 1:16

	% Get data
	if strcmp(type, 'Model')
		spikes_bass = 1000*nat_data(index).bass_spikerate{ind_b(target)}; % ms
		spikereps_bass = nat_data(index).bass_spikerep{ind_b(target)};
		spikes_oboe = 1000*nat_data(index).oboe_spikerate{ind_o(target)};
		spikereps_oboe = nat_data(index).oboe_spikerep{ind_o(target)};
	else
		spikes_bass = nat_data(index).bass_spikerate{ind_b(target)}; % ms
		spikereps_bass = nat_data(index).bass_spikerep{ind_b(target)};
		spikes_oboe = nat_data(index).oboe_spikerate{ind_o(target)};
		spikereps_oboe = nat_data(index).oboe_spikerep{ind_o(target)};
	end

	% Arrange data for SVM
	min_dis = 0.25;
	edges = 0:min_dis:300;
	t = 0+min_dis/2:min_dis:300-min_dis/2;
	for irep = 1:20
		h_bass1 = histcounts(spikes_bass(spikereps_bass==irep), edges);
		h_bass(irep, :) = h_bass1; %(randperm(length(h_bass1)));
		h_oboe1 = histcounts(spikes_oboe(spikereps_oboe==irep), edges);
		h_oboe(irep, :) = h_oboe1; %(randperm(length(h_bass1)));
	end
	h_all = [h_bass; h_oboe];
	h = [h; h_all];
end

% Put data into table
T = array2table(h);
T.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);

% 	case 'Model'
%
% 		h = [];
% 		for target = 1:16
%
% 			% Get data
% 			spikes_bass = nat_data(index).bass_PSTH_all{ind_b(target)}; % ms
% 			spikes_oboe = nat_data(index).oboe_PSTH_all{ind_o(target)};
%
% 			% Arrange data for SVM
% 			for irep = 1:20
% 				h_bass(irep, :) = spikes_bass(irep,:);
% 				h_oboe(irep, :) = spikes_oboe(irep,:);
% 			end
% 			h_all = [h_bass; h_oboe];
% 			h = [h; h_all];
% 		end
%
% 		% Put data into table
% 		T = array2table(h);
% 		T.Instrument = repmat([ones(20,1); ones(20, 1)*2], 16, 1);
% end

end
