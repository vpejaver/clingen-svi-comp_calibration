function lgnd = plot_both_posteriors_pub(postp, postb, thrs, Post_p, Post_b, method, flab, eshade)

%% Function to plot local posterior probabilities
% Note: (1) based heavily on 'plot_both_posteriors.m' in the 
%           'local_posterior_probability' directory
%       (2) there is a lot of hardcoded positional information here
%           (change at your own peril!)

% postp = posterior for pathogenicity
% postb = posterior for benignity
% thrs = thresholds to be examined
% Post_p = pathogenic posterior probability thresholds corresponding to the 4 strengths
% Post_b = benign posterior probability thresholds
% method = method name
% flab = A-Z figure labels (for publication purposes)
% eshade = information for "error" bars (gray curves in paper)
% lgnd = legend data structure is returned

% initialize
offset = max(thrs) - min(thrs); % set offset value based on score range

% pathogenic
if ~isempty(Post_p) % this check differentiates pathogenic from benign

    % plot the point estimate (although it says 'mean' it is using the full 
    % data set)
    h = plot(thrs, mean(postp, 1), 'color', 'black', 'linewidth', 1);
    hold on;

    % plot confidence bound curve
    if ~isempty(eshade{1}) & ~isempty(eshade{2})
        uci = eshade{1};
        lci = eshade{2};
        plot(thrs, uci, '-', 'color', '#808080', 'linewidth', 1);
        plot(thrs, lci, '-', 'color', '#808080', 'linewidth', 1);
    elseif ~isempty(eshade{1})
        uci = eshade{1};
        plot(thrs, uci, '-', 'color', '#808080', 'linewidth', 1);
    end

    % plot lines corresponding to posterior probability strength thresholds
    plot(thrs, ones(1, length(thrs)) * Post_p(1), '-', 'color', 'red', 'linewidth', 2);
    plot(thrs, ones(1, length(thrs)) * Post_p(2), '-.', 'color', 'red', 'linewidth', 2);
    plot(thrs, ones(1, length(thrs)) * Post_p(3), '--', 'color', 'red', 'linewidth', 2);
    plot(thrs, ones(1, length(thrs)) * Post_p(4), ':', 'color', 'red', 'linewidth', 2);

    % set axis labels
    axis([min(thrs) max(thrs) 0 1]);
    if strcmp(method, 'SIFT') | strcmp(method, 'FATHMM')
        xlabel('-Score');
    else
    	xlabel('Score');
    end
    ylabel('Posterior');

    % set x-axis tick information
    set(gca, 'Xtick', min(thrs):offset/5:max(thrs));
    set(gca, 'XtickLabel', strsplit(num2str(min(thrs):offset/5:max(thrs), '%.1f '), ' ')); %min(thrs):offset/5:max(thrs));
    set(gca, 'XtickLabelRotation', 45);

    % set y-axis tick information
    set(gca, 'Ytick', 0:1/5:1);
    set(gca, 'YtickLabel', strsplit(num2str(0:1/5:1, '%.1f '), ' '));
    
    % set title and gridline
    subtitle(method);
    set(gca, 'TitleHorizontalAlignment', 'left');
    text(min(thrs)-1.5*offset/5, 1.2, sprintf('%s', flab), 'fontsize', 16); %, 'fontweight', 'Bold');
    grid on;
    
    % set legend but don't show here
    lgnd = legend({'Point estimate', '', 'PP3\_VeryStrong', 'PP3\_Strong', 'PP3\_Moderate', 'PP3\_Supporting'}, 'location', 'west');
    set(lgnd, 'visible', 'off')
    hold off;
else   
    % benign
    thrs = thrs(length(thrs) : -1 : 1); % invert order of thresholds

    % plot the point estimate (although it says 'mean' it is using the full 
    % data set)
    h = plot(thrs, mean(postb, 1), 'color', 'black', 'linewidth', 1);
    hold on;

    % plot confidence bound curve
    if ~isempty(eshade{1}) & ~isempty(eshade{2})
        uci = eshade{1};
    	lci = eshade{2};
    	plot(thrs, uci, '-', 'color', '#808080', 'linewidth', 1);
    	plot(thrs, lci, '-', 'color', '#808080', 'linewidth', 1);
    elseif ~isempty(eshade{2})
        lci = eshade{2};
    	plot(thrs, lci, '-', 'color', '#808080', 'linewidth', 1);
    end

    % plot lines corresponding to posterior probability strength thresholds
    plot(thrs, ones(1, length(thrs)) * Post_b(1), '-', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);
    plot(thrs, ones(1, length(thrs)) * Post_b(2), '-.', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);
    plot(thrs, ones(1, length(thrs)) * Post_b(3), '--', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);
    plot(thrs, ones(1, length(thrs)) * Post_b(4), ':', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);

    % set axis labels (y-axis not labeled for benigns)
    axis([min(thrs) max(thrs) 0.975 1]);
    if strcmp(method, 'SIFT') | strcmp(method, 'FATHMM')
        xlabel('-Score');
    else
    	xlabel('Score');
    end

    % set x-axis tick information
    set(gca, 'Xtick', min(thrs):offset/5:max(thrs));
    set(gca, 'XtickLabel', strsplit(num2str(min(thrs):offset/5:max(thrs), '%.1f '), ' ')); %min(thrs):offset/5:max(thrs));
    set(gca, 'XtickLabelRotation', 45);
    
    % set y-axis tick information
    set(gca, 'Ytick', 0.975:(1-0.975)/5:1);
    set(gca, 'YtickLabel', strsplit(num2str(0.975:(1-0.975)/5:1, '%.3f '), ' ')); %0.975:(1-0.975)/5:1);
    
    % set legend but don't show here
    grid on;
    lgnd = legend({'', 'One-sided confidence bound', 'BP4\_VeryStrong', 'BP4\_Strong', 'BP4\_Moderate', 'BP4\_Supporting'}, 'location', 'west');
    set(lgnd, 'visible', 'off');
    hold off;
end

return