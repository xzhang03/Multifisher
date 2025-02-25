function [pvalues, ncomps] = multifisher(data, comps)
% Fisher's exact test corrected for multiple comparisons
% [pvalues, ncomps] = multifisher(data, comps)
%
% data -  n x 2 matrix where the rows are conditions, each row is entered 
%         as [hits, totals]
% comps - 0     (1D): all rows compared to all other rows
%         1-x   (1D): compare all rows to row 1-x
%         2D:   each row specify which two conditions are compared: 
%               e.g., [Condtion1, Condition2; Condition1, Condition 3]
% pvalues m x 3 matrix, where each row specify a fisher's exact test
%         e.g., [Condition1, Condition2, p-value (inflated with Bonferroni]
% ncomps  Number of comparisons

if nargin < 2
    comps = 0;
end

% Convert to table
datatable = array2table(data);

% Independent variables
nvar = size(datatable,1);

% Comparason table
if size(comps,2) == 2
    % 2D
    pvalues = comps;
    pvalues(:,3) = 0;
elseif size(comps,2) == 1
    % 1D
    switch comps
        case 0
            pvalues = zeros(nchoosek(nvar,2), 3);
            ind = 0;
            for i = 1 : nvar
                for j = i + 1 : nvar
                    ind = ind + 1;
                    pvalues(ind,1) = i;
                    pvalues(ind,2) = j;
                end
            end
        otherwise
            pvalues = (1 : nvar)';
            pvalues(pvalues == comps) = [];
            pvalues = cat(2, ones(nvar-1,1) * comps, pvalues, zeros(nvar-1,1));
    end
end
ncomps = size(pvalues, 1);

% Doing the comparisons
for i = 1 : ncomps
    datatable_subset = datatable(pvalues(i,[1 2]),:);
    [~,p,~] = fishertest(datatable_subset);
    pvalues(i, 3) = p;
end
pvalues(:,3) = pvalues(:,3) * ncomps;


end
