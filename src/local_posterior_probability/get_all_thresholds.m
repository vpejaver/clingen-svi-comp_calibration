function [thresh] = get_all_thresholds (posteriors, thrs, Post)

for b = 1 : size(posteriors, 1)
    for j = 1 : length(Post)
        ind = min(find(posteriors(b, :) < Post(j))) - 1;
        if ind > 0
            thresh(b, j) = thrs(ind);
        else
            thresh(b, j) = NaN;
        end
    end
end

return