% Vikas Pejaver
% February 2021
% University of Washington
% Script to plot heatmap of ClinVar test set varian proportions for
% each evidentiary criterion (one plot per tool)

clear
clc

%% Constants and defaults
in_file = '/Users/vikaspejaver/Desktop/ClinGen-SVI/results/correlations_gnomad.txt';
out_file = 'supp_correlation_gnomad.eps';
methods = {'BayesDel', 'CADD', 'EA', 'FATHMM', 'GERP++', 'MPC', 'MutPred2', ...
	    'PhyloP', 'PolyPhen-2', 'PrimateAI', 'REVEL', 'SIFT', 'VEST4'};

% Read in data file
[D, c] = tblread(in_file, '\t');
colnames = cellstr(c);

% Rorder
[~, ordr] = ismember(methods, colnames);
colnames = colnames(ordr);
M = triu(D) + tril(D');
for i = 1:size(M, 2)
    M(i, i) = 1;
end
M = M(ordr, ordr);
X = abs(tril(M));
X(X == 0) = NaN;

% Plot clustergram
f = figure('units', 'normalized', 'outerposition', [0 0 0.5 1])
P = uipanel('Parent', f, 'BorderType', 'none');
cmap = (flipud(bone));

h = heatmap(colnames, colnames, X, 'FontSize', 20, 'FontName', 'Arial', 'Parent', P);
h.CellLabelFormat = strrep('%.2f', '-', '');
%h.YDisplayLabels = repmat({''}, length(toolnames), 1);
h.ColorLimits = [0, 1];
h.Colormap = cmap;
h.MissingDataColor = [1, 1, 1];
colorbar off;

% Save plot
exportgraphics(P, out_file, 'Resolution', 300);
