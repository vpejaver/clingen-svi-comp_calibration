%% Vikas Pejaver
% Icahn School of Medicine at Mount Sinai and Univerity of Washington
% 2021-2022

%% Wrapper script to calculate local likelihood ratios for the default 
%% thresholds of three commonly used tools (Table 3 in paper)
% Input files: (1) List of MAT files output by main.m in 
%                  'local_posterior_probability'. Files located in
%                  'results -> bootstrapped'
%              (2) Directory that contains ClinVar 2019 and gnomAD data 
%                  sets, 'data -> data_combined_tools_split' (this is 
%                  mainly used to obtain all unique thresholds for a given
%                  tool as it is not stored in the output of main.m)
% Output file: None (prints to standard output)

%% Initialize
% clear screen and any standing variables in MATLAB workspace
clear
clc

%% Constants and defaults
indir = '/Users/vikaspejaver/Desktop/ClinGen-SVI/data/data_combined_tools_split/'; % change as needed
infiles = {'/Users/vikaspejaver/Desktop/ClinGen-SVI/results/new_bootstrapping/final_10k_3prct_100pts_cadd.mat',
	  '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/new_bootstrapping/final_10k_3prct_100pts_bd_pp2.mat',
	  '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/new_bootstrapping/final_10k_3prct_100pts.mat'};
toolnames = {'CADDphred', 'PolyPhen2-HVAR', 'SIFT'};
titles = {'CADD', 'PolyPhen-2', 'SIFT'};
def_thresh = [20, 0.902, 0.05];
priorp = 0.0441;

%% Loop through each MAT file and get information
fprintf(1, 'Method\tLR+\tCI\tFraction_gnomAD\n');
for n = 1:length(infiles)
    % load MAT file
    load(infiles{n});

    % get indices of the methods that were done (this was used only for
    % development purposes when results were partially done)
    i = find(strcmp(methods, toolnames{n}));
    
    % read in ClinVar variants
    D = load([indir files{i}]);
    x = D(:, 1);
    y = D(:, 2);
    
    % read in gnomAD variants
	file = strrep([indir files{i}], 'BLB', 'U');
    D = load(file);
    g = D(D(:, 2) == 0, 1);
    
    % thresholds for posterior for pathogenicity, reverse sorted
    thrs = flip(unique([x; g; floor(min([x; g])); ceil(max([x; g]))]));
    pcb = prctile(posteriors_p{i}(2:end, :), [2.5, 97.5]); % get CI

    % find threshold and local posterior probability value corresponding to
    % default threshold
    idx = find(thrs == def_thresh(n));
    this_post = posteriors_p{i}(1, idx); 
    
    % local likelihood ratio is the odds of posterior divided by odds of
    % prior
    lrpos = (this_post / (1 - this_post)) / (priorp / (1 - priorp));
    ci = (pcb(:, idx) / (1 - pcb(:, idx))) / (priorp / (1 - priorp));
    
    % get fraction of gnomAD variants that meet the default threshold for
    % pathogenicity
    if strcmp(methods{i}, 'SIFT')
        this_frac = length(find(g <= def_thresh(n))) / length(g);
    else
        this_frac = length(find(g >= def_thresh(n))) / length(g);
    end
    fprintf(1, '%s\t%.6f\t(%.6f, %.6f)\t%.6f\n', titles{n}, lrpos, ci(1), ci(2), this_frac);
end
