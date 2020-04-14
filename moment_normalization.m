function theNorm = moment_normalization(momentData)

% momentData is an N x 1 array of 4Pi moment data
% theNorm is the corresponding N x 1 array of normalization values

% turn the moment data into a square wave with values +1 or -1
segments = momentData./abs(momentData);

% shift the square wave and subtract to put a peak at each zero crossing
segments = abs(segments(1:end-1) - segments(2:end));

% find the index of each zero crossing - these indexes correspond to
% the location that the moment data (sine wave) becomes positive or goes
% negative. We want to normalize each of these segments to the max (or min)
% value within that segment
segmentIndex = find(segments);
segmentIndex = [1; segmentIndex];

% preallocate the array where we store the normalization
theNorm = zeros(size(momentData, 1), 1);

for i = 1:size(segmentIndex, 1) - 1
    
    % check if the data is above zero (find a max) or below zero (find a min)
    signTest = mean(momentData(segmentIndex(i):segmentIndex(i + 1)));
    
    % find the max (or min) value of the current segment and normalize all
    % the data in that segment to that max or min
    if signTest > 0
        theNorm(segmentIndex(i):segmentIndex(i + 1), 1) = max(momentData(segmentIndex(i):segmentIndex(i + 1)));
    else
        theNorm(segmentIndex(i):segmentIndex(i + 1), 1) = (-1)*min(momentData(segmentIndex(i):segmentIndex(i + 1)));
    end
    
end

% loop above does not address the last value(s) depending on the last peak
% position. Set the final value(s) manually
if signTest > 0
    theNorm(segmentIndex(end):end) = max(momentData(segmentIndex(i):segmentIndex(i + 1)));
else
    theNorm(segmentIndex(end):end) = (-1)*min(momentData(segmentIndex(i):segmentIndex(i + 1)));
end