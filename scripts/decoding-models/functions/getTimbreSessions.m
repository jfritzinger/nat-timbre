function [sesh, num_data] = getTimbreSessions(nat_data)
	
	% Find all rows with bassoon and oboe in them
	has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
	has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
	sesh = find(has_bass & has_oboe);
	num_data = length(sesh);

end