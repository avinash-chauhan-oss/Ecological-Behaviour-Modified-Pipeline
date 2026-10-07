function model = applyDiet(model, dietTable, compartment, dietName)
% APPLYDIET Enforces dietary exchange constraints on a metabolic model.
%
% INPUTS:
%   model       - COBRA model structure.
%   dietTable   - Struct with fields .rxns (cell/string array) and .lbs (numeric).
%   compartment - External compartment suffix (e.g., '(e)', '[u]').
%   dietName    - Diet condition name ('High_Fiber', 'Western', 'No_Diet').
%
% OUTPUTS:
%   model       - COBRA model with updated exchange bounds.

exIdx = startsWith(model.rxns, 'EX_');

if isempty(dietTable) || (nargin > 3 && strcmpi(dietName, 'No_Diet'))
    model.lb(exIdx) = -1000;
    model.ub(exIdx) = 1000;
    return;
end

model.lb(exIdx) = 0;
if ~exist('compartment', 'var') || isempty(compartment), compartment = '(e)'; end

dietTable.rxns = regexprep(dietTable.rxns, 'EX_H2\b', 'EX_h2');
dietTable.rxns = regexprep(dietTable.rxns, 'EX_glc\b', 'EX_glc_D');

dietRxns = strrep(dietTable.rxns, '(e)', compartment);
cleanDiet = regexprep(lower(string(dietRxns)), '\s+', '');
cleanModel = regexprep(lower(string(model.rxns)), '\s+', '');

[uniqueDietRxns, firstIdx] = unique(cleanDiet, 'stable');
if numel(uniqueDietRxns) < numel(cleanDiet)
    cleanDiet = cleanDiet(firstIdx);
    dietLbs = dietTable.lbs(firstIdx);
else
    dietLbs = dietTable.lbs;
end

[tf, loc] = ismember(cleanDiet, cleanModel);

model.lb(loc(tf)) = dietLbs(tf);
model.ub(loc(tf)) = 1000;

end