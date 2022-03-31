function [DiscountedThreshold] = get_discounted_thresholds (thresh, Post, B, discountonesided, type)

%% Function that selects the score threshold for each evidential strength
%% level using the one-sided confidence interval
% Note: This takes in the matrix of 'b' thresholds and returns a set of
% thresholds corresponding to the one-sided confidence bound (discount
% factor)

% thresh = selected threshold
% Post = posterior probability thresholds for the 4 strengths
% B = number of bootstrapping iterations
% discountonesided = the percentage for confidence bounds
% type = 'pathogenic' or 'benign'
% DiscountedThreshold = discounted thresholds

% loop through each posterior probability strength
for j = 1 : length(Post)
    % there will be some bootstrap samples for which a tool does not meet
    % the posterior probability threshold
    invalids = sum(isnan(thresh(:, j)));

    if invalids > discountonesided * B
        % if these exceed the confidence bound fraction, then the tool does
        % not meet the posterior probability threshold for the given strength
        DiscountedThreshold(j) = NaN;
    else
        % if not, then the 'invalids' must be excluded from the confidence
        % bound calculations
        if strcmp(type, 'pathogenic') == 1
            t = sort(thresh(find(isnan(thresh(:, j)) == 0), j), 'descend');
        elseif strcmp(type, 'benign') == 1
            t = sort(thresh(find(isnan(thresh(:, j)) == 0), j));
        end
      
        DiscountedThreshold(j) = t(floor(discountonesided * B) - invalids + 1);

    end
end

return