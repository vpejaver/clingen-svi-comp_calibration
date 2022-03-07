% Vikas Pejaver
% February 2021
% University of Washington
% Script to plot heatmap of ClinVar test set varian proportions for
% each evidentiary criterion (one plot per tool)

clear
clc

%% Constants and defaults
in_file = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/testset/wmp2_wpph2w_wcadd1.6_wBD_wEA_wnewpreds_rare2_wafs_wlabels_canonical_mappings_testset_notinalltrainingbdmpc_notinhgmd2019_diseaseonlyPgenelist_wvus_no0star_missense_clinvar_20201219.txt'; %'/Users/vikaspejaver/Desktop/ClinGen-SVI/results/final_files/wnewpreds_wpph2w_wcadd1.6_wmp2_wBD_wEA_wpreds_rare2_wafs_wlabels_canonical_mappings_defnotinalltrainingmpcbd_notinhgmd2019_diseaseonlyPgenelist_novus_no0star_missense_clinvar_20191202.txt'; '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/testset/wmp2_wpph2w_wcadd1.6_wBD_wEA_wnewpreds_rare2_wafs_wlabels_canonical_mappings_testset_notinalltrainingbdmpc_notinhgmd2019_diseaseonlyPgenelist_wvus_no0star_missense_clinvar_20201219.txt'; %'../../../results/testset/COL4A5_rare2.txt';
thr_file = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/discounted_new_thresholds_dcprior_nosmooth_bootstrap10000_3pc_posterior_output.txt';
%gene_file = '';
out_file = 'white_discounted_new_notinalltrainingbdmpc_dcprior_nosmooth_bootstrap10000_testset_lr.eps'; % 'discounted_new_notinalltrainingbdmpc_dcprior_nosmooth_bootstrap10000_trainingset_lr.eps'; %'COL4A5_nosmooth_bootstrap10_test_props_';
to_skip = {'BayesDel-wiAF', 'hEAt1.0', 'hEAt2.0', 'MutPred', 'phastCons100way'};

%classes = {'Pathogenic', 'Likely_pathogenic', 'Uncertain_significance', 'Likely_benign', 'Benign'};
classes = {'Pathogenic', 'Likely_pathogenic', 'Likely_benign', 'Benign'};

methods = {'BayesDel-noAF', 'BayesDel-wiAF', 'CADD', 'EA1.0', 'hEAt1.0', ...
    'hEAt2.0', 'MutPred', 'MutPred2.0', ...
    'PolyPhen-2', 'REVEL', 'VEST4', 'GERP++', ...
    'phastCons100way', 'phyloP100way', 'SIFT', 'FATHMM', 'MPC', 'PrimateAI'};
columns = {'BayesDel_nsfp33a_noAF', '', 'CADDv10x2E6_PHRED', 'EA_10x2E0', '', '', '', 'MutPred20x2E0_score', ...
    'pph2_prob', 'REVEL_score', 'VEST4_score', 'GERP0x2B0x2B_RS', ...
    '', 'phyloP100way_vertebrate', 'SIFT_score', 'FATHMM_score', 'MPC_score', 'PrimateAI_score'};

%xlabels = {'>= PVst', '[PSt, PVst)', '[PM, PSt)', '[PSu, PM)', 'Indet.', '[BM, BSu)', '[BSt, BM)', '[BVst, BSt)', '<BVst', 'Not pred.'};
%xlabels = {'>= PP3_VeryStrong', '[PP3_Strong, PP3_VeryStrong)', '[PP3_Moderate, PP3_Strong)', '[PP3_Supporting, PP3_Moderate)', 'Indet.', '[BP4_Moderate, BP4_Supporting)', '[BP4_Strong, BP4_Moderate)', '[BP4_VeryStrong, BP4_Strong)', '< BP4_VeryStrong', 'Not pred.'};
xlabels = {'PP3\_VeryStrong', 'PP3\_Strong', 'PP3\_Moderate', 'PP3\_Supporting', 'BP4\_Supporting', 'BP4\_Moderate', 'BP4\_Strong', 'BP4\_VeryStrong'};


% Read in data file
D = tdfread(in_file, '\t');

% Prepare labels
gene_labels = cellstr(D.Ensembl_geneid);
labels = cellstr(D.clnsig);

% Replace  problematic labels
labels = strrep(labels, 'Likely_benign,_drug_response', 'Likely_benign');
labels = strrep(labels, 'Benign/Likely_benign', 'Benign');
labels = strrep(labels, 'Pathogenic/Likely_pathogenic', 'Pathogenic');
labels = strrep(labels, 'Pathogenic,_risk_factor', 'Pathogenic');
labels = strrep(labels, 'Pathogenic/Likely_pathogenic,_risk_factor', 'Pathogenic');
labels = strrep(labels, 'Uncertain_significance,_other', 'Uncertain_significance');

% Get label counts
for i = 1:length(classes)
    N(i) = length(find(strcmp(labels, classes{i})));
end

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

% Loop through each predictor
k = 1;
to_remove = [];
lr = repmat([NaN], length(toolnames), length(xlabels));
cidx = find(strcmp(xlabels, 'PP3\_Supporting'));
for n = 1 : length(toolnames)
    if ismember(toolnames{n}, to_skip)
        continue;
    end
    fprintf(1, '\n%s', toolnames{n});
    fprintf(1, '\n');

    flipflag = 0;
    I = find(strcmp(methods, toolnames{n}));
    if strcmp(toolnames{n}, 'SIFT') | strcmp(toolnames{n}, 'FATHMM')
    	flipflag = 1;
    end

    if strcmp(columns{I}, '')
        to_remove = [to_remove n];
        continue;
    end
    if isfloat(D.(columns{I}))
        scores = D.(columns{I});
    else
        this_list = cellstr(D.(columns{I}));
    	this_list(strcmp(this_list, '.')) = {'NaN'};
    	num_list = cellfun(@str2num, this_list, 'uniformoutput', false);
    	num_list(cellfun(@isempty, num_list)) = {[NaN]};
    	scores = cell2mat(num_list);
    end

    inds = find(~isnan(scores));
    %inds = [1:length(scores)];
    scores = scores(inds);
    cats = labels(inds);
    genes = gene_labels(inds);

    start = find(~isnan(thrs(n, :)), 1, 'first');
    finish = find(~isnan(thrs(n, :)), 1, 'last');
    dpos = length(find(strcmp(cats, 'Pathogenic') | strcmp(cats, 'Likely_pathogenic')));
    dneg = length(find(strcmp(cats, 'Benign') | strcmp(cats, 'Likely_benign')));
    % PP3
    for j = start:cidx
        if j == start
    	    if flipflag %strcmp(toolnames{i}, 'SIFT') | strcmp(toolnames{i}, 'FATHMM')
    	        idx = find(scores <= thrs(n, j));
    	    else
    	        idx = find(scores >= thrs(n, j));
            end
        else
    	    if flipflag %strcmp(toolnames{i}, 'SIFT') | strcmp(toolnames{i}, 'FATHMM')
    	        idx = find(scores <= thrs(n, j) & scores > thrs(n, j-1));
    	    else
    	        idx = find(scores >= thrs(n, j) & scores < thrs(n, j-1));
    	    end
    	end
    	x = scores(idx);
    	g = genes(idx);
    	c = cats(idx);
    	npos = length(find(strcmp(c, 'Pathogenic') | strcmp(c, 'Likely_pathogenic')));
    	nneg = length(find(strcmp(c, 'Benign') | strcmp(c, 'Likely_benign')));
    	tpr = npos / dpos;
    	fpr = nneg / dneg;
    	lr(k, j) = tpr / fpr;
    end
    % BP4
    for j = finish:-1:cidx+1 
        if j == finish
            if flipflag %strcmp(toolnames{i}, 'SIFT') | strcmp(toolnames{i}, 'FATHMM')
    	        idx = find(scores >= thrs(n, j));
    	    else
    	        idx = find(scores <= thrs(n, j));
            end
        else
            if flipflag %strcmp(toolnames{i}, 'SIFT') | strcmp(toolnames{i}, 'FATHMM')
    	        idx = find(scores >= thrs(n, j) & scores < thrs(n, j+1));
    	    else
    	        idx = find(scores <= thrs(n, j) & scores > thrs(n, j+1));
    	    end
        end
        x = scores(idx);
    	g = genes(idx);
    	c = cats(idx);
    	npos = length(find(strcmp(c, 'Pathogenic') | strcmp(c, 'Likely_pathogenic')));
    	nneg = length(find(strcmp(c, 'Benign') | strcmp(c, 'Likely_benign')));
    	tpr = npos / dpos;
    	fpr = nneg / dneg;
    	lr(k, j) = tpr / fpr;
    end
    k = k + 1;
    %sum(sum(count(:, [5,10]))) / 222098
    %proportions = count ./ nansum(count, 2);

end

% Process names
toolnames(to_remove) = [];
toolnames = strrep(toolnames, 'BayesDel-noAF', 'BayesDel');
toolnames = strrep(toolnames, 'MutPred2.0', 'MutPred2');
toolnames = strrep(toolnames, 'phyloP100way', 'PhyloP');

% Plot heatmap
f = figure('units', 'normalized', 'outerposition', [0 0 1 1])
P = uipanel('Parent', f, 'BorderType', 'none');

% [0.64, 0.08, 0.18]
p(1) = subplot(1, 3, 1, 'Parent', P);
cmap = colormap(p(1), [linspace(1, 0.64)' linspace(1, 0.08)' linspace(1, 0.18)']);
tmp = lr(:, 1:cidx);
tmp(all(isnan(lr), 2), :) = [];
%tmp(:, 6:9) = -tmp(:, 6:9);
h = heatmap(xlabels(1:cidx), toolnames, tmp, 'FontSize', 20, 'FontName', 'Arial'); %proportions);
%h.Title = 'Likelihood ratio (odds of pathogenicity) in interval';
h.CellLabelFormat = strrep('%.2f', '-', '');
h.ColorLimits = [1, 100];
h.Colormap = cmap;
h.MissingDataColor = [1, 1, 1];
colorbar off;
%set(gca, 'FontSize', 20);
%set(gca, 'FontName', 'Arial');
pos = get(gca, 'Position');
pos = [pos(1)+0.065 pos(2)+0.03 pos(3:4)]
set(gca, 'Position', pos);

% [0.5, 0.5, 0.5]
p(2) = subplot(1, 3, 2, 'Parent', P);
cmap = colormap(p(2), repmat(linspace(1, 0.5)', 1, 3));
h = heatmap(repmat(0.5, length(toolnames), 1), 'CellLabelColor', 'none', 'GridVisible', 'off', 'FontSize', 20, 'FontName', 'Arial');
h.XDisplayLabels = '';
h.YDisplayLabels = repmat({''}, length(toolnames), 1);
h.ColorLimits = [0, 1];
h.Colormap = cmap;
colorbar off;
%set(gca, 'FontSize', 20);
%set(gca, 'FontName', 'Arial');
%h = heatmap(xlabels(idx), toolnames, tmp); %proportions);
%h.Title = 'Proportion of gnomAD variants in interval (%)';
%h.CellLabelFormat = strrep('%.2f', '-', '');
%h.YDisplayLabels = repmat({''}, length(toolnames), 1);
%h.ColorLimits = [0, 75];
%h.Colormap = cmap;
%colorbar off;
%set(gca, 'FontSize', 20);
pos = get(gca, 'Position');
pos = [pos(1) pos(2)+0.03 pos(3)-0.16 pos(4)]
set(gca, 'Position', pos);

% [0, 0.45, 0.74]
p(3) = subplot(1, 3, 3, 'Parent', P);
cmap = flipud(colormap(p(3), [linspace(1, 0)' linspace(1, 0.45)' linspace(1, 0.74)']));
tmp = lr(:, cidx+1:end);
tmp(all(isnan(lr), 2), :) = [];
%tmp(:, 6:9) = -tmp(:, 6:9);
h = heatmap(xlabels(cidx+1:end), toolnames, tmp, 'FontSize', 20, 'FontName', 'Arial'); %proportions);
h.CellLabelFormat = strrep('%.2f', '-', '');
h.YDisplayLabels = repmat({''}, length(toolnames), 1);
h.ColorLimits = [0, 1];
h.Colormap = cmap;
h.MissingDataColor = [1, 1, 1];
colorbar off;
%set(gca, 'FontSize', 20);
%set(gca, 'FontName', 'Arial');
pos = get(gca, 'Position');
pos = [pos(1)-0.225 pos(2)+0.03 pos(3:4)]
set(gca, 'Position', pos);

% Save plot
exportgraphics(P, out_file, 'Resolution', 1000);
%set(0,'DefaultFigureColor', [1 1 1]);
%set(gcf, 'Color', 'None');

%eval(['print -dpng ' out_file]);
%close(gcf);
