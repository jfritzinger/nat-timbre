function data = analyzeMTF(param)
% Calculate modulation transfer function (MTF) metrics from neural responses
%   data = analyzeMTF(param) analyzes neural data and returns an average
%   rate response to the MTF stimulus. Includes MTF shape classification.
%
%   Inputs:
%       param - Structure containing experimental parameters and response data
%
%   Outputs:
%       data - Structure containing MTF analysis results
%           Fields:
%           .fms        - Modulation frequencies analyzed
%           .rate       - Mean firing rate matrix (frequency × depth)
%           .rate_std   - Standard deviation of firing rates
%           .rlb        - Lower confidence bound of rates
%           .rub        - Upper confidence bound of rates
%           .BMF        - Best modulation frequency (Hz)
%           .WMF        - Worst modulation frequency (Hz)
%           .MTF_shape  - MTF classification ('BE, BS', etc.)
%           .at_100     - Response at 100 Hz modulation
%           .at_200     - Response at 200 Hz modulation
%           .rate_sm    - Smoothed rate matrix
%
%   See also: smooth_rates, accumstats, MTFclassification

this_ds = param.stims.dsid==param.dsid;

if param.dur > 10
	dur = param.dur/1000; % stimulus duration in seconds.
else
	dur = param.dur;
end

all_mod_depths = double([param.list.mdepth]).';
all_mod_depths(all_mod_depths < -80) = max(all_mod_depths);
[~,~,mdi] = unique(all_mod_depths);
[fms,~,fmi] = unique(double([param.list.fm]).');
if fms(1) == 0
	fms(1) = 1.2;
end

num_mod_freqs = length(fms);
num_depths = length(param.all_mdepths);
map_size = [num_mod_freqs, num_depths];

spike_rates = param.cluster.num_spikes_delayed(this_ds)/...
	(dur - param.onsetWin/1000);

if length(fmi)==length(spike_rates)
	[rate,rate_std, ~, rlb, rub]  = accumstats({fmi,mdi},spike_rates, map_size);
	rate_sm = smooth_rates(rate, rlb, rub, []);
	[BMF,WMF,MTF_shape, at_100, at_200] = MTFclassification(spike_rates,fms, fmi);
end

data.fms = fms;
data.rate = rate;
data.rate_std = rate_std;
data.rlb = rlb;
data.rub = rub;
data.BMF = BMF;
data.WMF = WMF;
data.MTF_shape = MTF_shape;
data.at_100 = at_100;
data.at_200 = at_200;
data.rate_sm = rate_sm;

end