function [postp, postb] = get_both_bootstrapped_posteriors (x, y, g, w, thrs, B, minpoints, gft, increment)

%% Function that calculates local posterior probabilities using a 
%% sliding window approach but on bootstrapped samples
% Note: this function essentially calls get_both_local_posteriors.m 'B' times

% x = scores on ClinVar from a prediction algorithm
% y = class labels on ClinVar, corresponding to scores x
% g = scores from gnomAD data
% w = weight of one negative data point
% thrs = thresholds to be examined
% B = number of bootstrap samples
% minpoints = min ClinVar points
% gft = gnomAD frequency threshold
% increment = by how much the halfwindow increases in each step
% postp = pathogenic posterior
% postb = benign posterior

% initialize
global to_parallelize;
global to_printout;

% prealocate memory for the posteriors
postp = zeros(B, length(thrs)); % pathogenic
postb = zeros(B, length(thrs)); % benign


if to_parallelize == false
    %
    % Boostrapping w/o parallelization; slower but uses less memory
    %
    
    for b = 1 : B
        % entertain the user, if the user wants it
        if to_printout == true && mod(b, 10) == 0
            fprintf(1, ['\n Iteration %d/%d'], b, B);
        end
        
        % boostrapped sample indices
        qx = randi(length(x), 1, length(x));
        qg = randi(length(g), 1, length(g));
        
        % compute posteriors
        postp(b, :) = get_both_local_posteriors(x(qx), y(qx), g(qg), thrs, ...
            w, minpoints, gft, increment);
        postb(b, :) = 1 - flip(postp(b, :)); % must reverse indices for benign
    end
else
    %
    % Boostrapping w/ parallelization; faster but uses more memory
    %
    
    % prealocate memory for bootstrap indices
    qx = zeros(B, length(x)); % for clinvar
    qg = zeros(B, length(g)); % for gnomad
    
    parfor b = 1 : B
        % entertain the user, if the user wants it
        if to_printout == true && mod(b, 10) == 0
            fprintf(1, ['\n Iteration %d/%d'], b, B);
        end
        
        % boostrapped sample indices
        qx(b, :) = randi(length(x), 1, length(x));
        qg(b, :) = randi(length(g), 1, length(g));
        
        % compute posteriors, pathogenic and benign
        postp(b, :) = get_both_local_posteriors(x(qx(b, :)), y(qx(b, :)), ...
            g(qg(b, :)), thrs, w, minpoints, gft, increment);
        postb(b, :) = 1 - flip(postp(b, :)); % must reverse indices for benign
    end
end

return