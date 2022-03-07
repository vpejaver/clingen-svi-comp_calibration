% Vikas Pejaver and Pedja Radivojac
% December 2020
% University of Washington and Northeastern University
% Script to select thresholds based on posteriors and bootstrapped data

clear
clc

%% Constants and defaults
infiles = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/new_bootstrapping/final_10k_3prct_100pts*.mat';
out_file = 'posterior_plots_13tools.eps';
toolnames = {'BayesDel-noAF', 'CADDphred', 'EA1.0', 'FATHMM', 'GERP++', 'MPC', 'MutPred2.0', 'phyloP100way', 'PolyPhen2-HVAR', 'PrimateAI', 'REVEL', 'SIFT', 'VEST4'};
titles = {'BayesDel', 'CADD', 'EA', 'FATHMM', 'GERP++', 'MPC', 'MutPred2', 'PhyloP', 'PolyPhen-2', 'PrimateAI', 'REVEL', 'SIFT', 'VEST4'};
fig_labels = 'A':'Z';
fig_labels = fig_labels(1:length(titles));
%toolnames = {'BayesDel-noAF', 'CADDphred', 'EA1.0', 'FATHMM', 'GERP++', 'MPC', 'MutPred2.0', 'PolyPhen2-HVAR', 'PrimateAI', 'REVEL', 'SIFT', 'VEST4'};
%titles = {'BayesDel', 'CADD', 'EA', 'FATHMM', 'GERP++', 'MPC', 'MutPred2', 'PolyPhen-2', 'PrimateAI', 'REVEL', 'SIFT', 'VEST4'};
shaded = 'CI'; % SE, CI or None
nrows = 7; %7
ncols = 4;

%% Prep figure
f = figure('units', 'normalized', 'outerposition', [0 0 0.33 1]);
P = uipanel('Parent', f, 'BorderType', 'none');

%% Loop through each MAT file and get information
filenames = strcat({dir(infiles).folder}, '/', {dir(infiles).name})';
for n = 1:length(filenames)
    n
    % Load MAT file
    load(filenames{n});

    % Get indices of the methods that were done
    idx = find(~cellfun(@isempty, ThresholdP));

    % Plot posteriors
    for i = idx

    	if ispc()
            D = load(['C:\Users\Predrag Radivojac\iCloudDrive\ClinGen\New Predictions\' files{i}]);
    	else
	    D = load(['/Users/vikaspejaver/Desktop/ClinGen-SVI/data/new_predictions_for_pedja/' files{i}]);
    	end
    
	fprintf(1, '\n%s', files{i});
    	fprintf(1, '\n');
    
	x = D(:, 1);
    	y = D(:, 2);
    
	if ispc()
            file = strrep(['C:\Users\Predrag Radivojac\iCloudDrive\ClinGen\New Predictions\' files{i}], 'BLB', 'U');
    	else
	    file = strrep(['/Users/vikaspejaver/Desktop/ClinGen-SVI/data/new_predictions_for_pedja/' files{i}], 'BLB', 'U');
    	end
    
	D = load(file);
    	g = D(D(:, 2) == 0, 1);
    
	% negate the scores for methods where high prediction means benign
    	if strcmp(methods{i}, 'SIFT') == 1 || strcmp(methods{i}, 'FATHMM') == 1
            x = -x;
            g = -g;
    	end

    	% Must remove one gnomAD variant with score > 1, as per Panos
    	if strcmp(methods{i}, 'EA1.0') == 1
            g = g(g <= 1);
    	end

        % thresholds for posterior for pathogenicity, reverse sorted
        thrs = flip(unique([x; g; floor(min([x; g])); ceil(max([x; g]))]));

	if strcmp(shaded, 'SE')
	    perrors = {posteriors_p{i}(1, :) - std(posteriors_p{i}(2:end, :)), posteriors_p{i}(1, :) + std(posteriors_p{i}(2:end, :))};
	    berrors = {posteriors_b{i}(1, :) - std(posteriors_b{i}(2:end, :)), posteriors_b{i}(1, :) + std(posteriors_b{i}(2:end, :))};
	elseif strcmp(shaded, 'CI')
            perrors = {prctile(posteriors_p{i}(2:end, :), 5), []};
            berrors = {[], prctile(posteriors_b{i}(2:end, :), 5)};
	else
	    perrors = {[], []};
	    berrors = {[], []};
	end

	% plot
	inds = find(strcmp(toolnames, methods{i}));
	this_title = titles{inds};
	subplot(nrows, ncols, 2*inds-1, 'Parent', P);
	lp = plot_both_posteriors_pub(posteriors_p{i}(1, :), posteriors_b{i}(1, :), thrs, Post_p, [], this_title, fig_labels(inds), perrors);
	subplot(nrows, ncols, 2*inds, 'Parent', P);
	lb = plot_both_posteriors_pub(posteriors_p{i}(1, :), posteriors_b{i}(1, :), thrs, [], Post_b, this_title, '', berrors);
	%end
	%break;
    end
    %break;
end

% Add legends
set(lp, 'visible', 'on');
set(lp, 'position', [0.5240 0.131 0.1448 0.0570]);
set(lb, 'visible', 'on');
set(lb, 'position', [0.6870 0.131 0.2152 0.0570]);

% Print to file
exportgraphics(P, out_file);
