function T = getF0PopTable(nat_data, target, sesh, F0s, num_data, type, type2, type3)

switch type2 % rate or timing
	case 'Rate'
		switch target
			case 'Bassoon'
				data_mat = NaN(length(F0s)*20, num_data);
				for ii = 1:num_data
					if strcmp(type, 'linear')
						X1 = nat_data(sesh(ii)).bass_norm_rep';
					else
						X1 = nat_data(sesh(ii)).bass_raterep';
					end
					X2 = reshape(X1', [], 1);
					data_mat(:,ii) = X2;
				end
				response = reshape(repmat(F0s, 1, 20)', 1, []);
				T = array2table(data_mat);
				T.Response = response';

			case 'Oboe'
				data_mat = NaN(length(F0s)*20, num_data);
				for ii = 1:num_data
					if strcmp(type, 'linear')
						X1 = nat_data(sesh(ii)).oboe_norm_rep';
					else
						X1 = nat_data(sesh(ii)).oboe_raterep';
					end
					X2 = reshape(X1', [], 1);
					data_mat(:,ii) = X2;
				end
				response = reshape(repmat(F0s, 1, 20)', 1, []);
				T = array2table(data_mat);
				T.Response = response';

			case 'Invariant'
				ind_b = 25:40;
				ind_o = [1 3:17];
				data_mat1 = NaN(2*20, num_data);
				data_mat = NaN(length(F0s)*40, num_data);
				for itarget = 1:16
					for ii = 1:num_data
						if strcmp(type, 'linear')
							X1 = nat_data(sesh(ii)).bass_norm_rep(:,ind_b(itarget));
							X2 = nat_data(sesh(ii)).oboe_norm_rep(:,ind_o(itarget));
						else
							X1 = nat_data(sesh(ii)).bass_raterep(:,ind_b(itarget));
							X2 = nat_data(sesh(ii)).oboe_raterep(:,ind_o(itarget));
						end
						X = [X1; X2];
						data_mat1(:,ii) = X;
					end
					idx = (1:40) + 40*(itarget-1);
					data_mat(idx, :) = data_mat1;
				end
				T = array2table(data_mat);
				B = repmat(F0s, 1, 40)';  % Repeat each number in a column 40 times
				T.Response = B(:);
		end
	case 'Timing'
		if strcmp(target, 'Bassoon') || strcmp(target, 'Oboe')

			h_all2 = [];
			for ineuron = 1:num_data
				h_all = [];
				for itarget = 1:length(F0s)
					if strcmp(target, 'Bassoon')
						spikes_bass = nat_data(sesh(ineuron)).bass_spikerate{itarget}; % ms
						spikereps_bass = nat_data(sesh(ineuron)).bass_spikerep{itarget};
					else
						spikes_bass = nat_data(sesh(ineuron)).oboe_spikerate{itarget}; % ms
						spikereps_bass = nat_data(sesh(ineuron)).oboe_spikerep{itarget};
					end
					if strcmp(type3, 'Model')
						spikes_bass = spikes_bass*1000;
					end


					% Arrange data for SVM
					min_dis = 1;
					edges = 0:min_dis:300;
					%t = 0+min_dis/2:min_dis:300-min_dis/2;
					for irep = 1:20
						h_bass(irep, :) = histcounts(spikes_bass(spikereps_bass==irep), edges);
					end
					h_all = [h_all; h_bass];
				end
				h_all2 = [h_all2, h_all];
			end
			response = reshape(repmat(F0s, 1, 20)', 1, []);
			T = array2table(h_all2);
			T.response = response';

		else
			ind_b = 25:40;
			ind_o = [1 3:17];
			h_all2 = [];
			for ineuron = 1:num_data
				h_all = [];
				for itarget = 1:length(ind_b)

					spikes_bass = nat_data(sesh(ineuron)).bass_spikerate{ind_b(itarget)}; % ms
					spikereps_bass = nat_data(sesh(ineuron)).bass_spikerep{ind_b(itarget)};
					spikes_oboe = nat_data(sesh(ineuron)).oboe_spikerate{ind_o(itarget)}; % ms
					spikereps_oboe = nat_data(sesh(ineuron)).oboe_spikerep{ind_o(itarget)};
					if strcmp(type3, 'Model')
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
					h_all = [h_all; h_bass; h_oboe];
				end
				h_all2 = [h_all2, h_all];
			end

			% Put data into table
			response = reshape(repmat(F0s, 1, 40)', 1, [])';
			T = array2table(h_all2);
			T.response = response;
		end
end
end