function [thresh] = get_all_thresholds (posteriors, thrs, Post)

%% Function that selects the score threshold for each evidential strength
%% level
% Note: This gets either the point estimate of the threshold (if
%       'posteriors' is a vector) or a matrix of 'b' thresholds (if
%       'posteriors' is a matrix), where 'b' is the number of bootstrap
%       iterations

% posteriors = posterior probability
% thrs = unique thresholds
% Post = posterior probability thresholds for the 4 strengths
% thresh = selected threshold

% loop through - if posteriors is a vector, only one set of thresholds will
% be returned; if posteriors is a matrix with 'b' rows, then 'b' sets of 
% thresholds will be returned
for b = 1 : size(posteriors, 1)
    % loop through the 4 strength levels (posterior thresholds)
    for j = 1 : length(Post)
        % find first score threshold with posterior probability less than
        % the posterior probability threshold (stores the index)
        ind = min(find(posteriors(b, :) < Post(j))) - 1;
        if ind > 0
            thresh(b, j) = thrs(ind); % get corresponding score
        else
            thresh(b, j) = NaN; % does not meet the strength level
        end
    end
end

return