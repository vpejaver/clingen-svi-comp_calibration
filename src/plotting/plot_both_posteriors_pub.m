function lgnd = plot_both_posteriors (postp, postb, thrs, Post_p, Post_b, method, flab, eshade)

offset = max(thrs) - min(thrs);
if ~isempty(Post_p) % pathogenic
    %subplot(2, 1, 1, 'Parent', P);
    % make a plot
    h = plot(thrs, mean(postp, 1), 'color', 'black', 'linewidth', 1);
    hold on;

    if ~isempty(eshade{1}) & ~isempty(eshade{2})
        uci = eshade{1};
	lci = eshade{2};
	plot(thrs, uci, '-', 'color', '#808080', 'linewidth', 1);
	plot(thrs, lci, '-', 'color', '#808080', 'linewidth', 1);
	%area(thrs, uci, 'facecolor', [0, 0, 0], 'facealpha', 0.05, 'edgecolor', 'black');
    	%area(thrs, lci, 'facecolor', [1, 1, 1], 'edgecolor', 'black');
    elseif ~isempty(eshade{1})
    	uci = eshade{1};
	plot(thrs, uci, '-', 'color', '#808080', 'linewidth', 1);
        %area(thrs, uci, 'facecolor', [0, 0, 0], 'facealpha', 0.05, 'edgecolor', 'black');
    end

    plot(thrs, ones(1, length(thrs)) * Post_p(1), '-', 'color', 'red', 'linewidth', 2);
    %text(min(thrs) +  0.01*offset, Post_p(1) - 0.02, 'VSt', 'fontsize', 8, 'color', 'red');
    plot(thrs, ones(1, length(thrs)) * Post_p(2), '-.', 'color', 'red', 'linewidth', 2);
    %text(min(thrs) +  0.01*offset, Post_p(2) - 0.02, 'St', 'fontsize', 8, 'color', 'red');
    plot(thrs, ones(1, length(thrs)) * Post_p(3), '--', 'color', 'red', 'linewidth', 2);
    %text(min(thrs) +  0.01*offset, Post_p(3) - 0.02, 'M', 'fontsize', 8, 'color', 'red');
    plot(thrs, ones(1, length(thrs)) * Post_p(4), ':', 'color', 'red', 'linewidth', 2);
    %text(min(thrs) +  0.01*offset, Post_p(4) - 0.02, 'Su', 'fontsize', 8, 'color', 'red');

    axis([min(thrs) max(thrs) 0 1]);
    if strcmp(method, 'SIFT') | strcmp(method, 'FATHMM')
        xlabel('-Score');
    else
	xlabel('Score');
    end
    ylabel('Posterior');
    set(gca, 'Xtick', min(thrs):offset/5:max(thrs));
    set(gca, 'XtickLabel', strsplit(num2str(min(thrs):offset/5:max(thrs), '%.1f '), ' ')); %min(thrs):offset/5:max(thrs));
    set(gca, 'XtickLabelRotation', 45);	
    set(gca, 'Ytick', 0:1/5:1);
    set(gca, 'YtickLabel', strsplit(num2str(0:1/5:1, '%.1f '), ' '));
    %pbaspect([1.5, 1, 1]);
    subtitle(method);
    set(gca, 'TitleHorizontalAlignment', 'left');
    text(min(thrs)-1.5*offset/5, 1.2, sprintf('%s', flab), 'fontsize', 16); %, 'fontweight', 'Bold');
    grid on;
    lgnd = legend({'Point estimate', '', 'PP3\_VeryStrong', 'PP3\_Strong', 'PP3\_Moderate', 'PP3\_Supporting'}, 'location', 'west');
    set(lgnd, 'visible', 'off')
    hold off;
else   % benign
    %subplot(2, 1, 2, 'Parent', P)
    thrs = thrs(length(thrs) : -1 : 1);

    % make a plot
    h = plot(thrs, mean(postb, 1), 'color', 'black', 'linewidth', 1);
    hold on;

    if ~isempty(eshade{1}) & ~isempty(eshade{2})
        uci = eshade{1};
	lci = eshade{2};
	plot(thrs, uci, '-', 'color', '#808080', 'linewidth', 1);
	plot(thrs, lci, '-', 'color', '#808080', 'linewidth', 1);
	%area(thrs, uci, 'facecolor', [0, 0, 0], 'facealpha', 0.05, 'edgecolor', 'black');
    	%area(thrs, lci, 'facecolor', [1, 1, 1], 'edgecolor', 'black');
    elseif ~isempty(eshade{2})
        lci = eshade{2};
	plot(thrs, lci, '-', 'color', '#808080', 'linewidth', 1);
    	%area(thrs, lci, 'facecolor', [1, 1, 1], 'edgecolor', 'black');
    end

    plot(thrs, ones(1, length(thrs)) * Post_b(1), '-', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);
    %text(max(thrs) - 0.07*offset, Post_b(1) - 0.0005, 'VSt', 'fontsize', 8, 'color', [0, 0.4470, 0.7410]);
    plot(thrs, ones(1, length(thrs)) * Post_b(2), '-.', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);
    %text(max(thrs) - 0.045*offset, Post_b(2) - 0.0005, 'St', 'fontsize', 8, 'color', [0, 0.4470, 0.7410]);
    plot(thrs, ones(1, length(thrs)) * Post_b(3), '--', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);
    %text(max(thrs) - 0.04*offset, Post_b(3) - 0.0005, 'M', 'fontsize', 8, 'color', [0, 0.4470, 0.7410]);
    plot(thrs, ones(1, length(thrs)) * Post_b(4), ':', 'color', [0, 0.4470, 0.7410], 'linewidth', 2);
    %text(max(thrs) - 0.055*offset, Post_b(4) - 0.0005, 'Su', 'fontsize', 8, 'color', [0, 0.4470, 0.7410]);

    axis([min(thrs) max(thrs) 0.975 1]);
    if strcmp(method, 'SIFT') | strcmp(method, 'FATHMM')
        xlabel('-Score');
    else
	xlabel('Score');
    end
    %ylabel('Posterior');
    set(gca, 'Xtick', min(thrs):offset/5:max(thrs));
    set(gca, 'XtickLabel', strsplit(num2str(min(thrs):offset/5:max(thrs), '%.1f '), ' ')); %min(thrs):offset/5:max(thrs));
    set(gca, 'XtickLabelRotation', 45);
    set(gca, 'Ytick', 0.975:(1-0.975)/5:1);
    set(gca, 'YtickLabel', strsplit(num2str(0.975:(1-0.975)/5:1, '%.3f '), ' ')); %0.975:(1-0.975)/5:1);
    %pbaspect([1.5, 1, 1]);
    grid on;
    lgnd = legend({'', 'One-sided confidence bound', 'BP4\_VeryStrong', 'BP4\_Strong', 'BP4\_Moderate', 'BP4\_Supporting'}, 'location', 'west');
    set(lgnd, 'visible', 'off');
    hold off;
end

return