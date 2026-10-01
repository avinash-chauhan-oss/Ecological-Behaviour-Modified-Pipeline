function model = applyDiet(model, dietTable, compartment, dietName)
% APPLYDIET Enforces specified dietary exchange constraints on a metabolic model.
% 
% INPUTS:
%   model       - COBRA model structure.
%   dietTable   - Table containing diet reactions and corresponding lower bounds.
%   compartment - String specifying the external compartment (e.g., '(e)').
%   dietName    - String specifying the diet condition name.
%
% OUTPUTS:
%   model       - COBRA model structure with updated exchange bounds.

exIdx = startsWith(model.rxns, 'EX_');

% Apply unconstrained baseline bounds if no diet is specified
if isempty(dietTable) || (nargin > 3 && strcmpi(dietName, 'No_Diet'))
    model.lb(exIdx) = -1000;
    model.ub(exIdx) = 1000;
    return;
end

% Constrain all external exchange reactions prior to applying diet specifics
model.lb(exIdx) = 0; 
if ~exist('compartment', 'var'), compartment = '(e)'; end

% Format reaction identifiers for matching
dietRxns = strrep(dietTable.rxns, '(e)', compartment);
cleanDiet = regexprep(lower(string(dietRxns)), '\s+', '');
cleanModel = regexprep(lower(string(model.rxns)), '\s+', '');

% Resolve duplicate reaction entries in the diet table
[uniqueDietRxns, firstIdx] = unique(cleanDiet, 'stable');
if numel(uniqueDietRxns) < numel(cleanDiet)
    cleanDiet = cleanDiet(firstIdx); 
    dietLbs = dietTable.lbs(firstIdx);
else
    dietLbs = dietTable.lbs;
end

[tf, loc] = ismember(cleanDiet, cleanModel);

% Apply mapped dietary lower bounds (VMH diets are already in mmol/gDW/h)
model.lb(loc(tf)) = dietLbs(tf);
model.ub(loc(tf)) = 1000;

end