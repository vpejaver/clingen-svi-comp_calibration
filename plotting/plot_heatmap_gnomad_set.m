% Vikas Pejaver
% January 2021
% University of Washington
% Script to plot heatmap of gnomAD variant proportions belonging to
% each evidentiary criterion

clear
clc

%% Constants and defaults
% Models where HIGHER predictions indicate MORE pathogenic
in_dir = '/Users/vikaspejaver/Desktop/ClinGen-SVI/data/new_predictions_for_pedja/';
thr_file = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/discounted_new_thresholds_dcprior_nosmooth_bootstrap10000_3pc_posterior_output.txt';
out_file = 'white_discounted_new_notinalltrainingbdmpc_dcprior_nosmooth_bootstrap10000_gnomad_props.eps';
fold_unpred = 1; % 1 means unpredicted variants will be considered indeterminate as well
N = 363894; %363920; % Actual number of variants
to_skip = {'BayesDel-wiAF', 'hEAt1.0', 'hEAt2.0', 'MutPred', 'phastCons100way'};

files = {...
    'BayesDel_nsfp33a_noAF_PLP_U_predictions.txt', 'BayesDel_nsfp33a_wiAF_PLP_U_predictions.txt', ...
    'CADDv1.6_PHRED_PLP_U_predictions.txt', 'EA_1.0_PLP_U_predictions.txt', ...
    'hEAt_1.0_PLP_U_predictions.txt', 'hEAt_2.0_PLP_U_predictions.txt', ...
    'MutPred_score_PLP_U_predictions.txt', 'MutPred2.0_score_PLP_U_predictions.txt', ...
    'pph2_prob_PLP_U_predictions.txt', 'REVEL_score_PLP_U_predictions.txt', ...
    'VEST4_score_PLP_U_predictions.txt', 'GERP++_RS_PLP_U_predictions.txt', ...
    'phastCons100way_vertebrate_PLP_U_predictions.txt', ...
    'phyloP100way_vertebrate_PLP_U_predictions.txt', 'MPC_score_PLP_U_predictions.txt', 'PrimateAI_score_PLP_U_predictions.txt'};

methods = {'BayesDel-noAF', 'BayesDel-wiAF', ...
    'CADD', 'EA1.0', ...
    'hEAt1.0', 'hEAt2.0', ...
    'MutPred', 'MutPred2.0', ...
    'PolyPhen-2', 'REVEL', ...
    'VEST4', 'GERP++', ...
    'phastCons100way', 'phyloP100way', 'MPC', 'PrimateAI'};
%xlabels = {'>= PVst', '[PSt, PVst)', '[PM, PSt)', '[PSu, PM)', 'Indet.', '[BM, BSu)', '[BSt, BM)', '[BVst, BSt)', '<BVst', 'Not pred.'};
%xlabels = {'>= PP3_VeryStrong', '[PP3_Strong, PP3_VeryStrong)', '[PP3_Moderate, PP3_Strong)', '[PP3_Supporting, PP3_Moderate)', 'Indet.', '[BP4_Moderate, BP4_Supporting)', '[BP4_Strong, BP4_Moderate)', '[BP4_VeryStrong, BP4_Strong)', '< BP4_VeryStrong', 'Not pred.'};
%xlabels = {'PP3_VeryStrong', 'PP3_Strong', 'PP3_Moderate', 'PP3_Supporting', 'Indeterminate', 'BP4_Supporting', 'BP4_Moderate', 'BP4_Strong', 'BP4_VeryStrong', 'Not predicted'};
if fold_unpred
     xlabels = {'PP3\_VeryStrong', 'PP3\_Strong', 'PP3\_Moderate', 'PP3\_Supporting', 'Indeterminate', 'BP4\_Supporting', 'BP4\_Moderate', 'BP4\_Strong', 'BP4\_VeryStrong'};
else
     xlabels = {'PP3\_VeryStrong', 'PP3\_Strong', 'PP3\_Moderate', 'PP3\_Supporting', 'Indeterminate', 'BP4\_Supporting', 'BP4\_Moderate', 'BP4\_Strong', 'BP4\_VeryStrong', 'Not predicted'};
end

% Models where HIGHER predictions indicate LESS pathogenic
selif = {'SIFT_score_PLP_U_predictions.txt', 'FATHMM_score_PLP_U_predictions.txt'};

sdohtem = {'SIFT', 'FATHMM'};

% Read in thresholds file
fid = fopen(thr_file, 'r');
if fid == -1
    error('ERROR: cannot open thresholds file!');
end
tline = fgetl(fid);
tline = fgetl(fid); % skip header line
toolnames = {};
thrs = [];
while ischar(tline)
    tokens = regexp(tline, '\t', 'split');
    toolnames = [toolnames, tokens{1}];
    this_thr = cellfun(@str2num, tokens(2:end));
    thrs = [thrs; [this_thr(1:4) this_thr(end:-1:5)]];
    tline = fgetl(fid);
end
fclose(fid);
toolnames = strrep(toolnames, 'PolyPhen2-HVAR', 'PolyPhen-2');
toolnames = strrep(toolnames, 'CADDphred', 'CADD');

%% Models where HIGHER predictions indicate MORE pathogenic
proportions = repmat(NaN, length(toolnames), length(xlabels));
for i = 1 : length(toolnames)
    if ismember(toolnames{i}, to_skip)
        continue;
    end
    flipflag = 0;
    I = find(strcmp(methods, toolnames{i}));
    if isempty(I)
        I = find(strcmp(sdohtem, toolnames{i}));
	this_file = selif{I};
	flipflag = 1; 
    else
        this_file = files{I};
    end
    
    D = load([in_dir this_file]);
    
    fprintf(1, '\n%s', toolnames{i});
    fprintf(1, '\n');

    %if flipflag
        %x = -D(D(:, 2) == 0, 1); %HIGHER predictions indicate LESS pathogenic
    %else
        x = D(D(:, 2) == 0, 1);
    %end 

    start = find(~isnan(thrs(i, :)), 1, 'first');
    finish = find(~isnan(thrs(i, :)), 1, 'last');
    for j = start:finish
        if j == start
	    if flipflag %strcmp(toolnames{i}, 'SIFT') | strcmp(toolnames{i}, 'FATHMM')
	        count = length(find(x <= thrs(i, j)));
	    else
	        count = length(find(x >= thrs(i, j)));
	    end
	else
	    if flipflag %strcmp(toolnames{i}, 'SIFT') | strcmp(toolnames{i}, 'FATHMM')
	        count = length(find(x <= thrs(i, j) & x > thrs(i, j-1)));
	    else
	        count = length(find(x >= thrs(i, j) & x < thrs(i, j-1)));
	    end
	end
	proportions(i, j) = (count * 100) / N;
    end
    if flipflag %strcmp(toolnames{i}, 'SIFT') | strcmp(toolnames{i}, 'FATHMM')
        count = length(find(x > thrs(i, finish)));
    else
        count = length(find(x < thrs(i, finish)));
    end
    
    proportions(i, finish+1) = (count * 100) / N;
    idx = find(strcmp(xlabels, 'Indeterminate'));
    if fold_unpred
	proportions(i, idx) = proportions(i, idx) + (((N - length(x)) * 100) / N);
    else
        proportions(i, end) = ((N - length(x)) * 100) / N;
    end
end

% Postprocess
inds = find(any(~isnan(proportions), 2));
proportions = proportions(inds, :);
toolnames = toolnames(inds);
toolnames = strrep(toolnames, 'BayesDel-noAF', 'BayesDel');
toolnames = strrep(toolnames, 'MutPred2.0', 'MutPred2');
toolnames = strrep(toolnames, 'phyloP100way', 'PhyloP');

% Added a 'flip' to match edits in the revised version of the paper
proportions = fliplr(proportions);
xlabels = fliplr(xlabels);

% Plot heatmap
f = figure('units','normalized','outerposition',[0 0 1 1])
P = uipanel('Parent',f,'BorderType','none');
%p.Title = 'Proportion of gnomAD variants in interval (%)';
%p.TitlePosition = 'lefttop';
%p.FontSize = 20;
%p.FontWeight = 'bold';

% [0.00, 0.45, 0.74]
p(1) = subplot(1, 3, 1, 'Parent', P);
cmap = colormap(p(1), [linspace(1, 0)' linspace(1, 0.45)' linspace(1, 0.74)']);
tmp = proportions(:, 1:idx-1);
%tmp(:, 6:9) = -tmp(:, 6:9);
h = heatmap(xlabels(1:idx-1), toolnames, tmp, 'FontSize', 40, 'FontName', 'Arial'); %proportions);
h.CellLabelFormat = strrep('%.2f', '-', '');
h.ColorLimits = [0, 75];
h.Colormap = cmap;
h.MissingDataColor = [1, 1, 1];
colorbar off;
%set(gca, 'FontSize', 20);
pos = get(gca, 'Position');
pos = [pos(1)+0.065 pos(2)+0.03 pos(3:4)]
set(gca, 'Position', pos);

%[0.5, 0.5, 0.5]
p(2) = subplot(1, 3, 2, 'Parent', P);
cmap = colormap(p(2), repmat(linspace(1, 0.5)', 1, 3));
tmp = proportions(:, idx);
%tmp(:, 6:9) = -tmp(:, 6:9);
h = heatmap(xlabels(idx), toolnames, tmp, 'FontSize', 40, 'FontName', 'Arial'); %proportions);
%h.Title = 'Proportion of gnomAD variants in interval (%)';
h.CellLabelFormat = strrep('%.2f', '-', '');
h.YDisplayLabels = repmat({''}, length(toolnames), 1);
h.ColorLimits = [0, 75];
h.Colormap = cmap;
colorbar off;
%set(gca, 'FontSize', 20);
pos = get(gca, 'Position');
pos = [pos(1) pos(2)+0.03 pos(3)-0.16 pos(4)]
set(gca, 'Position', pos);

% [0.64, 0.08, 0.18]
p(3) = subplot(1, 3, 3, 'Parent', P);
cmap = colormap(p(3), [linspace(1, 0.64)' linspace(1, 0.08)' linspace(1, 0.18)']);
tmp = proportions(:, idx+1:end);
%tmp(:, 6:9) = -tmp(:, 6:9);
h = heatmap(xlabels(idx+1:end), toolnames, tmp, 'FontSize', 40, 'FontName', 'Arial'); %proportions);
h.CellLabelFormat = strrep('%.2f', '-', '');
h.YDisplayLabels = repmat({''}, length(toolnames), 1);
h.ColorLimits = [0, 75];
h.Colormap = cmap;
h.MissingDataColor = [1, 1, 1];
colorbar off;
%set(gca, 'FontSize', 20);
pos = get(gca, 'Position');
pos = [pos(1)-0.225 pos(2)+0.03 pos(3:4)]
set(gca, 'Position', pos);

% Save plot
exportgraphics(P, out_file, 'Resolution', 300);
%set(0,'DefaultFigureColor', [1 1 1]);

%eval(['print -dpng ' out_file]);
%close(gcf);
