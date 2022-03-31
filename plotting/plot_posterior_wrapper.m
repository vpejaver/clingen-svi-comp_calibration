%% Vikas Pejaver and Predrag Radivojac
% Icahn School of Medicine at Mount Sinai, Univerity of Washington, and 
% Northeastern University
% 2021-2022

%% Wrapper script to plot publication-ready version of local posterior 
%% probability plots (Figure 3)
% Input files: (1) Path to MAT files output by main.m in 
%                  'local_posterior_probability' with the filename regular 
%                  expression (see below). Files located in
%                  'results -> bootstrapped'
%              (2) Directory that contains ClinVar 2019 and gnomAD data 
%                  sets, 'data -> data_combined_tools_split' (this is 
%                  mainly used to obtain all unique thresholds for a given
%                  tool as it is not stored in the output of main.m)
% Output file: A single EPS image file containing plots for all 13 pairs of
%              plots (output format can be changed as per MATLAB
%              documentation)
 
%% Initialize
% clear screen and any standing variables in MATLAB workspace
clear
clc

%% Constants and defaults
infiles = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/new_bootstrapping/final_10k_3prct_100pts*.mat';
indir = '/Users/vikaspejaver/Desktop/ClinGen-SVI/data/data_combined_tools_split/'; % change as needed
out_file = 'posterior_plots_13tools.eps';
toolnames = {'BayesDel-noAF', 'CADDphred', 'EA1.0', 'FATHMM', 'GERP++', 'MPC', 'MutPred2.0', 'phyloP100way', 'PolyPhen2-HVAR', 'PrimateAI', 'REVEL', 'SIFT', 'VEST4'};
titles = {'BayesDel', 'CADD', 'EA', 'FATHMM', 'GERP++', 'MPC', 'MutPred2', 'PhyloP', 'PolyPhen-2', 'PrimateAI', 'REVEL', 'SIFT', 'VEST4'};
fig_labels = 'A':'Z';
fig_labels = fig_labels(1:length(titles));
shaded = 'CI'; % 'SE' - standard errors, 'CI' - confidence interval-based, or 'None' - no gray line in the plot
onesided = 5; % If 'CI' is chosen, the percentage confidence interval that is calculated
nrows = 7;
ncols = 4;

%% Prep figure
f = figure('units', 'normalized', 'outerposition', [0 0 0.33 1]);
P = uipanel('Parent', f, 'BorderType', 'none');

%% Loop through each MAT file and plot
filenames = strcat({dir(infiles).folder}, '/', {dir(infiles).name})';
for n = 1:length(filenames)
    n
    % load MAT file (a single file may contain results for multiple tools)
    load(filenames{n});

    % get indices of the methods that were done (this was used only for
    % development purposes when results were partially done)
    idx = find(~cellfun(@isempty, ThresholdP));

    % loop through each tool 
    for i = idx
    	% read in ClinVar variants
        D = load([indir files{i}]);
	    fprintf(1, '\n%s', files{i});
    	fprintf(1, '\n');
    
	    x = D(:, 1);
    	y = D(:, 2);
    
        % read in gnomAD variants
	    file = strrep([indir files{i}], 'BLB', 'U');
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
        
        % calculate "error bars"; here represented as a grey line around
        % the point estimate
	    if strcmp(shaded, 'SE') % these will be symmetric
	        perrors = {posteriors_p{i}(1, :) - std(posteriors_p{i}(2:end, :)), posteriors_p{i}(1, :) + std(posteriors_p{i}(2:end, :))};
	        berrors = {posteriors_b{i}(1, :) - std(posteriors_b{i}(2:end, :)), posteriors_b{i}(1, :) + std(posteriors_b{i}(2:end, :))};
	    elseif strcmp(shaded, 'CI') % these are set so that the CI-based threshold is more stringent than the point estimate
            perrors = {prctile(posteriors_p{i}(2:end, :), onesided), []};
            berrors = {[], prctile(posteriors_b{i}(2:end, :), onesided)};
	    else
	        perrors = {[], []};
	        berrors = {[], []};
	    end

    	% get correct tool name (title)
    	inds = find(strcmp(toolnames, methods{i}));
    	this_title = titles{inds};

        % plot pathogenic
    	subplot(nrows, ncols, 2*inds-1, 'Parent', P);
    	lp = plot_both_posteriors_pub(posteriors_p{i}(1, :), posteriors_b{i}(1, :), thrs, Post_p, [], this_title, fig_labels(inds), perrors);
    	
        % plot benign
        subplot(nrows, ncols, 2*inds, 'Parent', P);
    	lb = plot_both_posteriors_pub(posteriors_p{i}(1, :), posteriors_b{i}(1, :), thrs, [], Post_b, this_title, '', berrors);
    end
end

% add legends (these are shown only once for the full figure)
set(lp, 'visible', 'on');
set(lp, 'position', [0.5240 0.131 0.1448 0.0570]);
set(lb, 'visible', 'on');
set(lb, 'position', [0.6870 0.131 0.2152 0.0570]);

% print to file
exportgraphics(P, out_file);
