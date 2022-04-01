%% Vikas Pejaver
% Icahn School of Medicine at Mount Sinai and Univerity of Washington
% 2021-2022

%% Wrapper script to calculate fractions of variants without prediction
%% scores for each tool (Supplemental Table S2 in paper)
% Input files: List of tab-delimited TXT files containing variants from
%              each of the three data sets. Files located in
%              'data -> ClinVar2019 [AND] ClinVar2020 [AND] gnomAD'
% Output file: Tab-delimited TXT file containing fractions of variants 
%              without predictions for each tool (rows) and in each data 
%              set (columns)

%% Initialize
% clear screen and any standing variables in MATLAB workspace
clear
clc

%% Constants and defaults
in_files = {'/Users/vikaspejaver/Desktop/ClinGen-SVI/results/final_files/wnewpreds_wpph2w_wcadd1.6_wmp2_wBD_wEA_wpreds_rare2_wafs_wlabels_canonical_mappings_defnotinalltrainingmpcbd_notinhgmd2019_diseaseonlyPgenelist_novus_no0star_missense_clinvar_20191202.txt', ...
    '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/final_files/wnewpreds_wpph2w_wcadd1.6_wmp2_wBD_wEA_downsampled_wpreds_rare2_wafs_canonical_mappings_defnotinalltrainingmpcbd_filthgmd2019_clinvar20191202_nolcr_merged_gnomad_only_DPGQfilt_splitmulti_missense_diseaseonlyPgenelist.txt', ...
    '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/testset/wmp2_wpph2w_wcadd1.6_wBD_wEA_wnewpreds_rare2_wafs_wlabels_canonical_mappings_testset_notinalltrainingbdmpc_notinhgmd2019_diseaseonlyPgenelist_novus_no0star_missense_clinvar_20201219.txt'};
out_file = 'data_coverage_table.txt';
methods = {'BayesDel-noAF', 'BayesDel-wiAF', 'CADD', 'EA1.0', 'hEAt1.0', ...
    'hEAt2.0', 'MutPred', 'MutPred2.0', ...
    'PolyPhen-2', 'REVEL', 'VEST4', 'GERP++', ...
    'phastCons100way', 'phyloP100way', 'SIFT', 'FATHMM', 'MPC', 'PrimateAI'};
columns = {'BayesDel_nsfp33a_noAF', '', 'CADDv10x2E6_PHRED', 'EA_10x2E0', '', '', '', 'MutPred20x2E0_score', ...
    'pph2_prob', 'REVEL_score', 'VEST4_score', 'GERP0x2B0x2B_RS', ...
    '', 'phyloP100way_vertebrate', 'SIFT_score', 'FATHMM_score', 'MPC_score', 'PrimateAI_score'};

%% Open output file
oid = fopen(out_file, 'w');
if oid == -1
    error('Cannot open output file!');
end

%% Loop through files and process
for n = 1:length(in_files)
    n

    % read in data file
    D = tdfread(in_files{n}, '\t');

    % prepare labels
    gene_labels = cellstr(D.Ensembl_geneid);

    if n ~= 2
        % replace problematic labels
        labels = cellstr(D.clnsig);
        labels = strrep(labels, 'Likely_benign,_drug_response', 'Likely_benign');
        labels = strrep(labels, 'Benign/Likely_benign', 'Benign');
       	labels = strrep(labels, 'Pathogenic/Likely_pathogenic', 'Pathogenic');
        labels = strrep(labels, 'Pathogenic,_risk_factor', 'Pathogenic');
        labels = strrep(labels, 'Pathogenic/Likely_pathogenic,_risk_factor', 'Pathogenic');
        labels = strrep(labels, 'Uncertain_significance,_other', 'Uncertain_significance');
    end

    % loop through tools
    for i = 1:length(columns)
        if ~strcmp(columns{i}, '')
    	    if isfloat(D.(columns{i})) % this is easy as data type is numeric
                scores = D.(columns{i});
            else % format scores if not numeric data type
        		this_list = cellstr(D.(columns{i}));
        		this_list(strcmp(this_list, '.')) = {'NaN'};
        		num_list = cellfun(@str2num, this_list, 'uniformoutput', false);
        		num_list(cellfun(@isempty, num_list)) = {[NaN]};
        		scores = cell2mat(num_list);
            end
    	    denom = length(scores);

            % find all variants with NaN (no scores)
    	    inds = find(isnan(scores));
            scores = scores(inds);
    	    if n ~= 2
     	       cats = labels(inds);
    	    end
    	    genes = gene_labels(inds);
    	    num = length(scores);
            
            % calculate percentage of missing variants
            cov(i, n) = num/denom * 100;
        end
    end
end

%% Print to output
fprintf(oid, 'Method\tClinVar 2019\tgnomAD\tClinVar 2020\n');
for i = 1:length(columns)
    if ~strcmp(columns{i}, '')
        fprintf(oid, '%s\t%.1f\t%.1f\t%.1f\n', methods{i}, cov(i, 1), cov(i, 2), cov(i, 3));
    end
end
fclose(oid);
