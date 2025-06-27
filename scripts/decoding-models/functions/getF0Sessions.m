function [sesh, num_data] = getF0Sessions(nat_data, type)

if strcmp(type, 'Bassoon')
	sesh = find(~cellfun(@isempty, {nat_data.bass_rate}));
elseif strcmp(type, 'Oboe')
	sesh = find(~cellfun(@isempty, {nat_data.oboe_rate}));
elseif strcmp(type, 'Invariant')
	has_bass = ~cellfun(@isempty, {nat_data.bass_rate});
	has_oboe = ~cellfun(@isempty, {nat_data.oboe_rate});
	sesh = find(has_bass & has_oboe);
end
num_data = length(sesh);

end