%% Predrag Radivojac and Vikas Pejaver
% Northeastern University, Icahn School of Medicine at Mount Sinai and
% Univerity of Washington
% 2021-2022

%% Wrapper script to compute local posterior probabilities, 
%% plot rough versions of posterior probability plots, and 
%% estimate thresholds (full and bootstrapped)
% Input files: For each tool, the ClinVar 2019 P/LP+B/LB variants and 
%              gnomAD variants. Files located in 
%              'data -> data_combined_tools_split'
% Output file: A single MAT file for all 13 tools (or however many are
%              looped through below) containing all variables in this 
%              wrapper script. Files located in 'results -> bootstrapped'
 
%% Initialize
% clear screen and any standing variables in MATLAB workspace
clear
clc

% randomize seed for bootstrapping
rng('shuffle');

%% Output file name
filetosave = 'cadd_primateai_10k_3prct_100pts';

%% Parallel thread setup
% invokes parfor when true, otherwise the code is not parallelized
global to_parallelize;
to_parallelize = true;

% set the desired number of parallel threads
parallelthreads = 4;

% invokes printouts when true, otherwise the user is not entertained
global to_printout;
to_printout = true;

% enable parallel threads
if to_parallelize == true
    p = gcp('nocreate');
    if isempty(p)                  % if no parallel pool active
        parpool(parallelthreads);  % create one
    elseif p.NumWorkers ~= parallelthreads % if no. threads incorrect
        delete(gcp('nocreate'));           % delete the pool
        parpool(parallelthreads);          % then create one
        % error will be thrown if the no. threads is too large
    end
end


%% Files and tool names
indir = '/Users/vikaspejaver/Desktop/ClinGen-SVI/data/data_combined_tools_split/'; % change as needed
files = {'BayesDel_nsfp33a_noAF_PLP_BLB_predictions.txt', 'EA_1.0_PLP_BLB_predictions.txt', ...
    'MutPred2.0_score_PLP_BLB_predictions.txt', 'pph2_prob_PLP_BLB_predictions.txt', ...
    'REVEL_score_PLP_BLB_predictions.txt', 'VEST4_score_PLP_BLB_predictions.txt', ...
    'GERP++_RS_PLP_BLB_predictions.txt', 'phyloP100way_vertebrate_PLP_BLB_predictions.txt', ...
    'SIFT_score_PLP_BLB_predictions.txt', 'FATHMM_score_PLP_BLB_predictions.txt', ...
    'CADDv1.6_PHRED_PLP_BLB_predictions.txt', 'MPC_score_PLP_BLB_predictions.txt', ...
    'PrimateAI_score_PLP_BLB_predictions.txt'};

methods = {'BayesDel-noAF', 'EA1.0', 'MutPred2.0', ...
    'PolyPhen2-HVAR', 'REVEL', 'VEST4', 'GERP++', ...
    'phyloP100way', 'SIFT', 'FATHMM', 'CADDphred', ...
    'MPC', 'PrimateAI'};


%% Tuneable parameters
c = 1124;       % LR+ constant from Tavtigian et al. (changes w/ alpha)
alpha = 0.0441; % Prior probability of pathogenicity (changes w/ c)
B = 10000;      % Number of boostrapping iterations
discountonesided = 0.05;     % discount for one-sided 95th upper confidence interval
windowclinvarpoints = 100;   % minimum number of clinvar variants in a local window
windowgnomadfraction = 0.03; % minimum fraction of gnomad variants in a local window

% increment size for each method (default value is 0.001)
% for some methods we know what the increment should be, because we
% precomputed them; for others we have to assume something small
%
% this is needed because different tools output scores on different scale
% and with different granularity
increments = 0.001 * ones(1, length(methods));
for i = 1 : length(methods)
    if strcmp(methods{i}, 'BayesDel-noAF') == 1
        increments(i) = 0.01;
    elseif strcmp(methods{i}, 'CADDraw') == 1
        increments(i) = 0.01;
    elseif strcmp(methods{i}, 'CADDphred') == 1
        increments(i) = 0.1;
    elseif strcmp(methods{i}, 'EA1.0') == 1
        increments(i) = 0.001;
    elseif strcmp(methods{i}, 'MutPred') == 1
        increments(i) = 0.001;
    elseif strcmp(methods{i}, 'MutPred2.0') == 1
        increments(i) = 0.001;
    elseif strcmp(methods{i}, 'PolyPhen2-HVAR') == 1
        increments(i) = 0.001;
    elseif strcmp(methods{i}, 'REVEL') == 1
        increments(i) = 0.001;
    elseif strcmp(methods{i}, 'VEST4') == 1
        increments(i) = 0.001;
    elseif strcmp(methods{i}, 'GERP++') == 1
        increments(i) = 0.01;
    elseif strcmp(methods{i}, 'phyloP100way') == 1
        increments(i) = 0.01;
    elseif strcmp(methods{i}, 'SIFT') == 1
        increments(i) = 0.001;
    elseif strcmp(methods{i}, 'FATHMM') == 1
        increments(i) = 0.01;
    end
end

%
%
%
%% Below this part there's no parameters to set
%
%
%
%%
% Thresholds for posteriors, p and b, for alpha given above
% 1 = VS, 2 = ST, 3 = MO, 4 = SU
for j = 1 : 4
    Post_p(j) = c ^ (1 / 2 ^ (j - 1)) * alpha / ...
        ((c ^ (1 / 2 ^ (j - 1)) - 1) * alpha + 1);
    Post_b(j) = (c ^ (1 / 2 ^ (j - 1))) * (1 - alpha) / ...
        (((c ^ (1 / 2 ^ (j - 1))) - 1) * (1 - alpha) + 1);
end


%% Loop through and estimate intervals
for i = 1 : length(files) % does all 13 tools but change as needed
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

    % the worth of a negative example to satisfy the prior alpha
    w = (1 - alpha) * sum(y == 1) / (sum(y == 0) * alpha);
    
    % thresholds for posterior for pathogenicity, reverse sorted
    thrs = flip(unique([x; g; floor(min([x; g])); ceil(max([x; g]))]));
        
    tic
    % first row corresponds to original data
    posteriors_p{i}(1, :) = ...
        get_both_local_posteriors(x, y, g, thrs, w, ...
        windowclinvarpoints, windowgnomadfraction, increments(i)); 
    toc
    
    % benign posteriors are 1 - pathogenic posteriors, then reverse ordered
    posteriors_b{i}(1, :) = 1 - flip(posteriors_p{i}(1, :));
    
    % plot the posteriors to entertain the user, must negate thresholds for
    % methods where high scores mean benign and low scores mean pathogenic
    if strcmp(methods{i}, 'SIFT') == 1 || strcmp(methods{i}, 'FATHMM') == 1
        plot_both_posteriors(posteriors_p{i}(1, :), posteriors_b{i}(1, :), -thrs, Post_p, Post_b, methods{i});
    else
        plot_both_posteriors(posteriors_p{i}(1, :), posteriors_b{i}(1, :), thrs, Post_p, Post_b, methods{i});
    end
    
    % rows 2 to B+1 correspond to bootstrapped data posteriors
    [posteriors_p{i}(2 : B + 1, :), posteriors_b{i}(2 : B + 1, :)] = ...
        get_both_bootstrapped_posteriors(x, y, g, w, thrs, B, ...
        windowclinvarpoints, windowgnomadfraction, increments(i));
    
    % get thresholds for the regular & bootstrapped posteriors (pathogenic)
    pthresh{i} = get_all_thresholds(posteriors_p{i}, thrs, Post_p);
    
    % get thresholds for the regular & bootstrapped posteriors (benign)
    bthresh{i} = get_all_thresholds(posteriors_b{i}, flip(thrs), Post_b);
    
    % obtain nonbootstrapped thresholds (pathogenic, benign)
    ThresholdP{i} = pthresh{i}(1, :);
    ThresholdB{i} = bthresh{i}(1, :);
    
    % obtain discounted (one-sided confidence bound-based) thresholds (pathogenic, benign)
    DiscountedThresholdP{i} = ...
        get_discounted_thresholds(pthresh{i}(2 : B + 1, :), Post_p, B, discountonesided, 'pathogenic');
    DiscountedThresholdB{i} = ...
        get_discounted_thresholds(bthresh{i}(2 : B + 1, :), Post_b, B, discountonesided, 'benign');
    
    % unnegate thresholds for  methods where scores were negated
    if strcmp(methods{i}, 'SIFT') == 1 || strcmp(methods{i}, 'FATHMM') == 1
        pthresh{i} = -pthresh{i};
        bthresh{i} = -bthresh{i};
        ThresholdP{i} = -ThresholdP{i};
        ThresholdB{i} = -ThresholdB{i};
        DiscountedThresholdP{i} = -DiscountedThresholdP{i};
        DiscountedThresholdB{i} = -DiscountedThresholdB{i};
    end

    % print thresholds on the screen
    print_thresholds(ThresholdP{i}, ThresholdB{i}, DiscountedThresholdP{i}, DiscountedThresholdB{i});
    
    % give computer a 3-minute rest    
    pause(180);
    
    eval(['save ' filetosave ' -v7.3']);
end

fprintf(1, '\n');
fprintf(1, '\n');
