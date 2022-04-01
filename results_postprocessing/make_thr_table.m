%% Vikas Pejaver
% Icahn School of Medicine at Mount Sinai and Univerity of Washington
% 2021-2022

%% Wrapper script to extract threshold information and print formatted 
%% output (Table 2 and Supplemental Table S1 in the paper)
% Input files: Path to MAT files output by main.m in 
%              'local_posterior_probability' with the filename regular 
%              expression (see below). Files located in
%              'results -> bootstrapped'
% Output file: A tab-delimited TXT file containing thresholds formatted to 
%              look as close to these two tables in the paper

%% Initialize
% clear screen and any standing variables in MATLAB workspace
clear
clc

%% Constants and defaults
infiles = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/new_bootstrapping/1pc/final_10k_1prct_100pts*.mat';
out_file = 'new_thresholds_dcprior_nosmooth_bootstrap10000_1pc_posterior_output.txt';
choice = 'discounted'; % Either actual, mean, or discounted (only actual and discounted used in the paper)
toolnames = {'BayesDel-noAF', 'CADDphred', 'EA1.0', 'FATHMM', 'GERP++', 'MPC', 'MutPred2.0', 'phyloP100way', 'PolyPhen2-HVAR', 'PrimateAI', 'REVEL', 'SIFT', 'VEST4'};

%% Open output file
oid = fopen([choice, '_', out_file], 'w');
if oid == -1
    error('ERROR: cannot open output file!');
end

%% Loop through each MAT file and get information
thresholds = [];
tools = {};

filenames = strcat({dir(infiles).folder}, '/', {dir(infiles).name})';
for n = 1:length(filenames)
    n
    % load MAT file
    load(filenames{n});

    % get indices of the methods that were done (this was used only for
    % development purposes when results were partially done)
    idx = find(~cellfun(@isempty, ThresholdP));

    % put thresholds together
    for i = idx
        if strcmp(choice, 'mean')
    	    thresholds = [thresholds; [nanmean(pthresh{i}(2:end, :)) nanmean(bthresh{i}(2:end, :))]];
    	elseif strcmp(choice, 'actual')
    	    thresholds = [thresholds; [ThresholdP{i} ThresholdB{i}]];
    	else
    	    thresholds = [thresholds; [DiscountedThresholdP{i} DiscountedThresholdB{i}]];
    	end
    	tools = [tools; methods{i}];
    end
end

%% Reorder
[~, ordr] = ismember(toolnames, tools);
ordr(ordr == 0) = [];
tools = tools(ordr);
thresholds = thresholds(ordr, :);

%% Write to file
fprintf(oid, 'Method\tPP3_VeryStrong\t\tPP3_Strong\t\tPP3_Moderate\t\tPP3_Supporting\t\tBP4_VeryStrong\t\tBP4_Strong\t\tBP4_Moderate\t\tBP4_Supporting\n');
for i = 1:length(tools)
    values = regexprep(sprintf('%.3f\t', thresholds(i, :)), '\t$', '');
    fprintf(oid, '%s\t%s\n', tools{i}, values);
end
fclose(oid);
