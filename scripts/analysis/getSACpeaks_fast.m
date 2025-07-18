function [peakstats, peakstats_detail] = getSACpeaks_fast(spktrainset,CIin,analwindow,F0,minspkpersweep,plotdir,unitid)
% getSACpeaks Robust analysis of peaks in a single SAC with bootstrap statistics

% Handle plotting
% plotopt = false;
% if nargin>5 && ~isempty(plotdir)
%     saveplot = true;
% else
%     saveplot = false;
% end

% User progress message
fprintf('getSACpeaks: ');
if nargin>6
    fprintf(' (#%d) ',unitid);
end
fprintf('F0: %dHz. \n',F0);

duration = diff(analwindow);

% Params for peak picking
am_period = 1e3/F0;
wfr = 1/50; % Smoothing window fraction of period
maxlagfactor = 1.5;
maxlag_ms = maxlagfactor*1e3/F0;

% Initialize outputs
peakstats_detail = struct();
peakstats = struct();

% Defensive checks
if isempty(CIin)
    fprintf('Empty correlation index input. Stopping.\n');
    return;
end

allspks = [spktrainset{:}];
numSweeps = size(spktrainset,1);
spikesPerSweep = length(allspks)/numSweeps;
if spikesPerSweep<minspkpersweep
    fprintf('Not enough spikes for analysis. Stopping.\n');
    return;
end
if F0>1000
    fprintf('AM frequency >1kHz. Stopping.\n');
    return;
end

% Analysis window setup
wlen = round(1e3*wfr/(F0*(CIin.lags(2)-CIin.lags(1))))+1;
isel = find(abs(CIin.lags)<maxlag_ms*1.1);
lagaxis = CIin.lags(isel);

% Smoothing
if wlen>1
    smoothsac = filtfilt(ones(1,wlen)/wlen,1,CIin.unsmoothedSAC(isel));
else
    smoothsac = CIin.unsmoothedSAC(isel);
end

% --- Peak picking
[peaksets, trofsets, salsets] = findmainpeaks_fast(smoothsac);

% --- Bootstrap null distribution
randn = 500;
fprintf('\nComputing bootstrapped distributions of peak saliences:\n');

randpeaksets = cell(randn,1);
randtrofsets = cell(randn,1);
randsaliences = cell(randn,1);
allrandsmoothsacs = cell(randn,1);

for ri = 1:randn
    % Shuffle interspike intervals across all sweeps
    rand_spikesets = cell(length(spktrainset),1);
    isis = cell(length(spktrainset),1);

    for sweepi = 1:length(spktrainset)
        isis{sweepi} = diff(spktrainset{sweepi});
    end
    isiarray = [isis{:}];
    surrograteisis = isiarray(randperm(length(isiarray)));

    counter = 1;
    for sweepi = 1:length(spktrainset)
        inc = length(spktrainset{sweepi})-1;
        if isscalar(spktrainset{sweepi})
            rand_spikesets{sweepi} =  [spktrainset{sweepi}(1)];
        elseif length(spktrainset{sweepi})>1
            rand_spikesets{sweepi} =  cumsum([spktrainset{sweepi}(1) surrograteisis(counter:counter+inc-1)]);
            counter = counter+inc;
        end
    end

    % Compute the SAC
    [randh2plot, randbc2plot] = SPTCORR(rand_spikesets,'nodiag',1.1*maxlag_ms,.05,duration,'LouageNorm');

    if ri == 1
        wlen = round(1e3*wfr/(F0*(randbc2plot(2)-randbc2plot(1))))+1;
        %isel = find(abs(randbc2plot)<maxlag_ms*1.1);
    end

    if wlen>1
        randsmoothsac = filtfilt(ones(1,wlen)/wlen,1,randh2plot, "ctf");
    else
        randsmoothsac = randh2plot;
    end
    allrandsmoothsacs{ri} = randsmoothsac;

    if rem(ri,50)==0
        fprintf('.'); 
    end
end

for ri = 1:randn
    [randpeaksets{ri}, randtrofsets{ri}, randsaliences{ri}]=...
		findmainpeaks_fast(allrandsmoothsacs{ri});
end
randsmoothsaceg = allrandsmoothsacs{1};

% --- Robust null distribution stats (main fix here!)
maxnpeaks = max(cellfun(@length, randsaliences)); % determine largest possible set
salvals = cell(1, maxnpeaks); % ROBUST CHANGE

saltmp = [randsaliences{:}];
for i=1:length(saltmp)
    if isempty(saltmp{i}), continue; end
    sallen =  length(saltmp{i});
    if isempty(salvals{sallen})
        salvals{sallen} = [];
    end
    sval = min(saltmp{i});
    if isempty(sval) || isnan(sval), sval = NaN; end
    salvals{sallen} = [salvals{sallen} sval]; % ROBUST CHANGE
end

% --- Histogram stats/quantiles/pvals
salbins = 0:0.01:10;
salhist = zeros(length(salvals),length(salbins));
salpval = zeros(length(salvals),length(salbins));
salmean = nan(1,length(salvals));
salmedian = nan(1,length(salvals));
salsd = nan(1,length(salvals));
salquantiles = nan(length(salvals),5);

for i=1:length(salvals)
    if ~isempty(salvals{i})
        salhist(i,:) = hist(salvals{i},salbins);
        salpval(i,:) = 1- cumsum(salhist(i,:))/sum(salhist(i,:));
        salmean(i) = mean(salvals{i}, 'omitnan');
        salmedian(i) = median(salvals{i}, 'omitnan');
        salsd(i) = std(salvals{i}, 'omitnan');
        salquantiles(i,:) = quantile(salvals{i},[.025 .25 .50 .75 .975]);
    end
end

% --- Peak p-values: robust data structure and indexing
saliences = cell(1, length(peaksets));
peakpvalues = cell(1, length(peaksets));
for pi = 1:length(peaksets)
    setlen = sum(~isnan(peaksets{pi}));
    saliences{setlen} = salsets{pi};
    peakpvalues{setlen} = nan(1, setlen);
    for pii = 1:setlen
        pind = find(salbins<saliences{setlen}(pii),1,'last');
        if ~isempty(pind)
            peakpvalues{setlen}(pii) = salpval(setlen,pind);
        end
    end
end

% --- Assign outputs: all the structure assignments, as in your code
peakstats.wlen =wlen;
peakstats.amfreq = F0;
peakstats.am_period = am_period;
peakstats.peaksets =peaksets;
peakstats.trofsets = trofsets;
peakstats.salsets =salsets;
peakstats.maxlag_ms = maxlag_ms;
peakstats.smoothsac =smoothsac;
peakstats.randsmoothsaceg =randsmoothsaceg;
peakstats.peakpvalues = peakpvalues;
peakstats.lagaxis = lagaxis;
peakstats.spikespersweep = spikesPerSweep;
peakstats.randsals.mean = salmean;
peakstats.randsals.median = salmedian;
peakstats.randsals.sd = salsd;
peakstats.randsals.quantiles = salquantiles;
peakstats.randsals.quantile_values = [.025 .25 .50 .75 .975];

% Significance thresholds (as in your code, robust to nans)
% get_sigidx = @(pv,th) find(cellfun(@(x) max(x),pv,'Uniform',false) < th); 
get_sigidx = @(pv,th) find(cellfun(@(x) ~isempty(x) && isnumeric(x) && max(x) < th, pv));

peakstats.setofsigpeaks = [];
peakstats.setofsigpeaklags = [];
pt = get_sigidx(peakpvalues, .05);
if ~isempty(pt)
    peakstats.setofsigpeaks = peaksets{max(pt)}(~isnan(peaksets{max(pt)}));
end
if ~isempty(peakstats.setofsigpeaks)
    peakstats.setofsigpeaklags = randbc2plot(peakstats.setofsigpeaks);
end

peakstats.setofsigpeaks_p01 = [];
peakstats.setofsigpeaklags_p01 = [];
pt = get_sigidx(peakpvalues, .01);
if ~isempty(pt)
    peakstats.setofsigpeaks_p01 = peaksets{max(pt)}(~isnan(peaksets{max(pt)}));
end
if ~isempty(peakstats.setofsigpeaks_p01)
    peakstats.setofsigpeaklags_p01 = randbc2plot(peakstats.setofsigpeaks_p01);
end

peakstats.setofsigpeaks_p001 = [];
peakstats.setofsigpeaklags_p001 = [];
pt = get_sigidx(peakpvalues, .001);
if ~isempty(pt)
    peakstats.setofsigpeaks_p001 = peaksets{max(pt)}(~isnan(peaksets{max(pt)}));
end
if ~isempty(peakstats.setofsigpeaks_p001)
    peakstats.setofsigpeaklags_p001 = randbc2plot(peakstats.setofsigpeaks_p001);
end

peakstats.setofsigpeaks_p0001 = [];
peakstats.setofsigpeaklags_p0001 = [];
pt = get_sigidx(peakpvalues, .0001);
if ~isempty(pt)
    peakstats.setofsigpeaks_p0001 = peaksets{max(pt)}(~isnan(peaksets{max(pt)}));
end
if ~isempty(peakstats.setofsigpeaks_p0001)
    peakstats.setofsigpeaklags_p0001 = randbc2plot(peakstats.setofsigpeaks_p0001);
end

peakstats_detail.salhist = salhist;
peakstats_detail.salpval = salpval;
peakstats_detail.randpeaksets = randpeaksets;
peakstats_detail.randtrofsets = randtrofsets;
peakstats_detail.randsaliences = randsaliences;
peakstats_detail.allrandsmoothsacs = allrandsmoothsacs;

fprintf('\tdone.\n');

% --- Plotting (identical to your code, omitted here...) ---


end
