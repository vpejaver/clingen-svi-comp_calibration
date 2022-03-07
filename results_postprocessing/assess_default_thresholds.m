% Vikas Pejaver and Pedja Radivojac
% December 2020
% University of Washington and Northeastern University
% Script to select thresholds based on posteriors and bootstrapped data

clear
clc

%% Constants and defaults
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
    % Load MAT file
    load(infiles{n});

    % Get indices of the methods that were done
    i = find(strcmp(methods, toolnames{n}));

    % Plot posteriors
    if ispc()
        D = load(['C:\Users\Predrag Radivojac\iCloudDrive\ClinGen\New Predictions\' files{i}]);
    else
	D = load(['/Users/vikaspejaver/Desktop/ClinGen-SVI/data/new_predictions_for_pedja/' files{i}]);
    end
    
    x = D(:, 1);
    y = D(:, 2);
    
    if ispc()
        file = strrep(['C:\Users\Predrag Radivojac\iCloudDrive\ClinGen\New Predictions\' files{i}], 'BLB', 'U');
    else
	file = strrep(['/Users/vikaspejaver/Desktop/ClinGen-SVI/data/new_predictions_for_pedja/' files{i}], 'BLB', 'U');
    end
    
    D = load(file);
    g = D(D(:, 2) == 0, 1);
    
    % thresholds for posterior for pathogenicity, reverse sorted
    thrs = flip(unique([x; g; floor(min([x; g])); ceil(max([x; g]))]));
    pcb = prctile(posteriors_p{i}(2:end, :), [2.5, 97.5]);

    idx = find(thrs == def_thresh(n));
    this_post = posteriors_p{i}(1, idx); %mean(posteriors_p{i}(2:end, idx)); 
    lrpos = (this_post / (1 - this_post)) / (priorp / (1 - priorp));
    ci = (pcb(:, idx) / (1 - pcb(:, idx))) / (priorp / (1 - priorp));
    if strcmp(methods{i}, 'SIFT')
        this_frac = length(find(g <= def_thresh(n))) / length(g);
    else
        this_frac = length(find(g >= def_thresh(n))) / length(g);
    end
    fprintf(1, '%s\t%.6f\t(%.6f, %.6f)\t%.6f\n', titles{n}, lrpos, ci(1), ci(2), this_frac);
end
