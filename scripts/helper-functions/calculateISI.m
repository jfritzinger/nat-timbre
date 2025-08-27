function calculateISI()

	ISI = arrayfun(@(ii) diff(spike_times_on(reps==ii)), 1:nrep, 'UniformOutput', false);
	ISI_all = vertcat(ISI{:});

end