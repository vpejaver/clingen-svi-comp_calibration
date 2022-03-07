% Vikas Pejaver and Pedja Radivojac
% December 2020
% University of Washington and Northeastern University
% Script to select thresholds based on posteriors and bootstrapped data

clear
clc

%% Constants and defaults
infiles = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/new_bootstrapping/1pc/final_10k_1prct_100pts*.mat';
out_file = 'new_thresholds_dcprior_nosmooth_bootstrap10000_1pc_posterior_output.txt';
choice = 'discounted'; % Either actual, mean, or discounted 
toolnames = {'BayesDel-noAF', 'CADDphred', 'EA1.0', 'FATHMM', 'GERP++', 'MPC', 'MutPred2.0', 'phyloP100way', 'PolyPhen2-HVAR', 'PrimateAI', 'REVEL', 'SIFT', 'VEST4'};

%% Open output file
%oid = 1;
oid = fopen([choice, '_', out_file], 'w');
if oid == -1
    error('ERROR: cannot open output file!');
end

%% Loop through each MAT file and get information
thresholds = [];
%post_point = [];
%post_lower = [];
%post_upper = [];
tools = {};

filenames = strcat({dir(infiles).folder}, '/', {dir(infiles).name})';
for n = 1:length(filenames)
    n
    % Load MAT file
    load(filenames{n});

    % Get indices of the methods that were done
    idx = find(~cellfun(@isempty, ThresholdP));

    % Put thresholds together
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
%return

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
