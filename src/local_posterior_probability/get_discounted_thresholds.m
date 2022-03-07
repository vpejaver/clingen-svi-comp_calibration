function [DiscountedThreshold] = get_discounted_thresholds (thresh, Post, B, discountonesided, type)

% determine thresholds and discounted thresholds
for j = 1 : length(Post)
    invalids = sum(isnan(thresh(:, j)));
    if invalids > discountonesided * B
        DiscountedThreshold(j) = NaN;
    else
        if strcmp(type, 'pathogenic') == 1
            t = sort(thresh(find(isnan(thresh(:, j)) == 0), j), 'descend');
        elseif strcmp(type, 'benign') == 1
            t = sort(thresh(find(isnan(thresh(:, j)) == 0), j));
        end
       
        DiscountedThreshold(j) = t(floor(discountonesided * B) - invalids + 1);
        
        %DiscountedThreshold(j) = t(ceil(discountonesided * B) - invalids);
        %DiscountedThreshold(j) = ...
        %    (t(floor(discountonesided * B) - invalids) + ...
        %     t(floor(discountonesided * B) - invalids + 1)) / 2;
    end
end

return