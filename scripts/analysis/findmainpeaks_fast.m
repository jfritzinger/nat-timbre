function [peaksets, trofsets, salsets] = findmainpeaks_fast(x)
% Fast version of findmainpeaks5
% Picks salient peaks and troughs using built-in findpeaks

minProminence = 1e-6; % Adjust as needed for salience

% Detect peaks (maxima) and troughs (minima) using findpeaks
[~, peakIdx] = findpeaks(x, 'MinPeakProminence', minProminence);          % Peaks
[~, trofIdx] = findpeaks(-x, 'MinPeakProminence', minProminence);         % Troughs
trofIdx = trofIdx(:)';
peakIdx = peakIdx(:)';

if isempty(peakIdx) || isempty(trofIdx)
    peaksets = {};
    trofsets = {};
    salsets = {};
    return
end

% First set: global max and min
[~, maxPeakIdx] = max(x(peakIdx));
[~, minTrofIdx] = min(x(trofIdx));
peaksets = {peakIdx(maxPeakIdx)};
trofsets = {trofIdx(minTrofIdx)};
salsets = {x(peakIdx(maxPeakIdx)) - x(trofIdx(minTrofIdx))};
minsals = salsets{1};

% Preallocate for typical set size up to N (can be tuned)
maxSets = min(length(peakIdx), length(trofIdx));
peaksets(1, maxSets) = {[]};
trofsets(1, maxSets) = {[]};
salsets(1, maxSets) = {[]};

k = 2;
while true
    peaksAvail = setdiff(peakIdx, [peaksets{1:k-1}]);
    trofsAvail = setdiff(trofIdx, [trofsets{1:k-1}]);
    if isempty(peaksAvail) || isempty(trofsAvail)
        break
    end
    [pGrid, tGrid] = meshgrid(peaksAvail, trofsAvail);
    if isempty(pGrid) || isempty(tGrid)
        break;
    end
    salCandidates = x(pGrid) - x(tGrid);
    if isempty(salCandidates) || all(isnan(salCandidates(:)))
        break;
    end
    [maxSal, idx] = max(salCandidates(:));
    if isempty(idx) || maxSal == -Inf
        break;
    end
    [row, col] = ind2sub(size(salCandidates), idx);
    if col > numel(peaksAvail) || row > numel(trofsAvail)
        break;
    end
    nextPeak = peaksAvail(col);
    nextTrof = trofsAvail(row);
    % Update sets
    peaksets{k} = sort([peaksets{k-1}, nextPeak]);
    trofsets{k} = sort([trofsets{k-1}, nextTrof]);
    salsets{k} = x(peaksets{k}) - x(trofsets{k});
    minsals(k) = min(salsets{k});
    k = k + 1;
end

lastValid = find(~cellfun(@isempty, peaksets), 1, 'last');
peaksets = peaksets(1:lastValid);
trofsets = trofsets(1:lastValid);
salsets = salsets(1:lastValid);

end
