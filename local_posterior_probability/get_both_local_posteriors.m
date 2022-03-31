function [post, finalwindow] = get_both_local_posteriors (x, y, g, thrs, w, minpoints, gft, increment)

%% Function that actually calculates local posterior probabilities using a 
%% sliding window approach
% x = scores on ClinVar from a prediction algorithm
% y = class labels on ClinVar, corresponding to scores x
% g = scores from gnomAD data
% thrs = thresholds to be examined
% w = weight of one negative data point
% minpoints = min ClinVar points
% gft = gnomAD frequency threshold
% increment = by how much the halfwindow increases in each step
% post = posterior
% finalwindow = returns final halfwindow size for debugging purposes

% Initialize
post = zeros(1, length(thrs));         % posterior for pathogenicity
finalwindow = zeros(1, length(thrs));  % final window size for each threshold

maxthrs = max(thrs);
minthrs = min(thrs);
lengthgnomad = length(g);

% Loop through each unique threshold
for i = 1 : length(thrs)
    % half window initially takes a 0 value as some predictors can have a 
    % lot of variants with identical scores
    halfwindow = 0;

    % hi and lo boundaries change whenever halfwindow takes a new value
    lo = thrs(i) - halfwindow; 
    hi = thrs(i) + halfwindow;

    % now we need to find the window size such that enough points are in it
    while 1
        % number of positive and negative points in the window
        pos = length(find(x >= lo & x <= hi & y == 1));
        neg = length(find(x >= lo & x <= hi & y == 0));

        % for windows that are near score boundaries (max or min) we need
        % to reduce the effective number of points proportionally to the
        % reduced interval size. c <= 1, with c = 1 meaning that the full
        % window is considered
        if hi > maxthrs
            c = (maxthrs - lo) / (hi - lo);
        elseif lo < minthrs
            c = (hi - minthrs) / (hi - lo);
        else
            c = 1;
        end
        
        % check just in case whether c is indeed positive
        if c <= 0
            error('Problem wih computing c');
        end

        % if the number of + and - points is too small, increase the window
        % size by an increment (2 * increment for the full window)
        if pos + neg < c * minpoints
            halfwindow = halfwindow + increment;
            lo = thrs(i) - halfwindow;
            hi = thrs(i) + halfwindow;
            continue
        end

        % if the number of positives and negatives is large enough, we need
        % to check whether there are enough gnomAD points in that window
        if length(find(g >= lo & g <= hi)) < gft * c * lengthgnomad
            gnomad_condition = 0;
        else
            gnomad_condition = 1;
        end
        
        if pos + neg >= c * minpoints && gnomad_condition == 1
            break
        else
            halfwindow = halfwindow + increment;
            lo = thrs(i) - halfwindow;
            hi = thrs(i) + halfwindow;
        end
    end
    
    % compute the posterior for each threshold
    post(i) = pos / (pos + w * neg);

    % compute minimum delta for each threshold (mostly for debugging)
    finalwindow(i) = halfwindow;
end

return