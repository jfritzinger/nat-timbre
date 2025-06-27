function T = getF0PopTable(nat_data, target, sesh, F0s, num_data, type)

if strcmp(target, 'Bassoon')
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

	% Create table for model
	
	response = reshape(repmat(F0s, 1, 20)', 1, []);
	T = array2table(data_mat);
	T.Response = response';


elseif strcmp(target, 'Oboe')
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

	% Create table for model
	response = reshape(repmat(F0s, 1, 20)', 1, []);
	T = array2table(data_mat);
	T.Response = response';


elseif strcmp(target, 'Invariant')

	% Model including all F0s
	ind_b = 25:40;
	ind_o = [1 3:17];
	data_mat = NaN(2*20, num_data);
	data_mat2 = NaN(length(F0s)*40, num_data);
	for itarget = 1:16
		for ii = 1:num_data

			% Arrange data for SVM
			if strcmp(type, 'linear')
				X1 = nat_data(sesh(ii)).bass_norm_rep(:,ind_b(itarget));
				X2 = nat_data(sesh(ii)).oboe_norm_rep(:,ind_o(itarget));
			else
				X1 = nat_data(sesh(ii)).bass_raterep(:,ind_b(itarget));
				X2 = nat_data(sesh(ii)).oboe_raterep(:,ind_o(itarget));
			end
			X = [X1; X2];
			data_mat(:,ii) = X;
		end
		idx = (1:40) + 40*(itarget-1);
		data_mat2(idx, :) = data_mat;
	end

	% Create table for model
	T = array2table(data_mat2);
	B = repmat(F0s, 1, 40)';  % Repeat each number in a column 40 times
	T.Response = B(:);
end

end