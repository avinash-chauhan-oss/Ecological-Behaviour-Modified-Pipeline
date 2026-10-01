% MINIMAL_REACTOMES
% 
% Executes parallel extraction of minimal reactomes for a batch of metabolic 
% models across specified dietary conditions. Viability failures and 
% extinctions are recorded as NaN.

OverallTimer = tic;
fprintf('\n-----------------------------------------------------\n');
fprintf('GENERATING MINIMAL REACTOMES\n');
fprintf('-----------------------------------------------------\n');

base_dir = pwd;
diet_dir = fullfile(base_dir, 'Diets');
output_dir = fullfile(base_dir, 'minReactModels');

if ~exist(output_dir, 'dir')
    mkdir(output_dir); 
end

changeCobraSolver('gurobi', 'LP', 0);
changeCobraSolver('gurobi', 'MILP', 0);

model_names = {
    'Acidaminococcus_fermentans_DSM_20731.mat'; 'Akkermansia_muciniphila_ATCC_BAA_835.mat';
    'Alistipes_finegoldii_DSM_17242.mat'; 'Alistipes_putredinis_DSM_17216.mat'; 'Alistipes_shahii_WAL_8301.mat';
    'Bacteroides_caccae_ATCC_43185.mat'; 'Bacteroides_coprocola_M16_DSM_17136.mat'; 'Bacteroides_eggerthii_DSM_20697.mat';
    'Bacteroides_fragilis_NCTC_9343.mat'; 'Bacteroides_intestinalis_341_DSM_17393.mat'; 'Bacteroides_stercoris_ATCC_43183.mat';
    'Bacteroides_thetaiotaomicron_VPI_5482.mat'; 'Bacteroides_uniformis_ATCC_8492.mat'; 'Bifidobacterium_adolescentis_ATCC_15703.mat';
    'Bifidobacterium_bifidum_PRL2010.mat'; 'Bifidobacterium_longum_NCC2705.mat'; 'Bilophila_wadsworthia_3_1_6.mat';
    'Blautia_wexlerae_DSM_19850.mat'; 'Citrobacter_amalonaticus_Y19.mat'; 'Clostridium_clostridioforme_CM201.mat';
    'Collinsella_tanakaei_YIT_12063.mat'; 'Coprococcus_catus_GD_7.mat'; 'Desulfovibrio_desulfuricans_subsp_desulfuricans_DSM_642.mat';
    'Desulfovibrio_piger_ATCC_29098.mat'; 'Dialister_invisus_DSM_15470.mat'; 'Dorea_longicatena_DSM_13814.mat';
    'Enterobacter_cloacae_EcWSU1.mat'; 'Enterococcus_faecium_TX1330.mat'; 'Escherichia_coli_str_K_12_substr_MG1655.mat';
    'Eubacterium_rectale_M104_1.mat'; 'Faecalibacterium_prausnitzii_L2_6.mat'; 'Fusobacterium_varium_ATCC_27725.mat';
    'Haemophilus_parainfluenzae_T3T1.mat'; 'Helicobacter_pylori_26695.mat'; 'Klebsiella_pneumoniae_pneumoniae_MGH78578.mat';
    'Lactobacillus_mucosae_LM1.mat'; 'Lactobacillus_reuteri_SD2112_ATCC_55730.mat'; 'Megasphaera_elsdenii_DSM_20460.mat';
    'Parabacteroides_distasonis_ATCC_8503.mat'; 'Parabacteroides_johnsonii_DSM_18315.mat'; 'Phascolarctobacterium_succinatutens_YIT_12067.mat';
    'Prevotella_ruminicola_23.mat'; 'Pseudoflavonifractor_capillosus_strain_ATCC_29799.mat'; 'Roseburia_intestinalis_L1_82.mat';
    'Roseburia_inulinivorans_DSM_16841.mat'; 'Ruminococcus_bromii_L2_63.mat'; 'Ruminococcus_callidus_ATCC_2776001.mat';
    'Ruminococcus_torques_ATCC_27756.mat'; 'Shigella_flexneri_2002017.mat'; 'Streptococcus_salivarius_JIM8777.mat';
    'Subdoligranulum_variabile_DSM_15176.mat'; 'Veillonella_atypica_ACS_049_V_Sch6.mat'
};
model_files = struct('name', model_names, 'isdir', num2cell(false(size(model_names))));
%
nModels = numel(model_files);

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
nDiets = numel(dietNames);
DietTables = cell(nDiets-1, 1);

% Load dietary constraints from text files
for d = 1:nDiets-1
    file = fullfile(diet_dir, [dietNames{d} '.txt']);
    T = readtable(file, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);
    DietTables{d}.rxns = strrep(strtrim(string(T{:,1})), '[u]', '(e)');
    DietTables{d}.lbs  = T{:,2};
end

% Initialize output directories for each diet condition
for d = 1:nDiets
    folder = fullfile(output_dir, dietNames{d}); 
    if ~exist(folder, 'dir')
        mkdir(folder); 
    end
end

if isempty(gcp('nocreate'))
    parpool('local'); 
end

ParpoolResults = cell(nModels, 1);
fprintf('\nExecuting parallel reduction across %d models...\n', nModels);

parfor m = 1:nModels
    changeCobraSolver('gurobi', 'all', 0);
    ModelTimer = tic;
    file = model_files(m).name; 
    [~, modelName] = fileparts(file);
    localResults = cell(nDiets, 7); 

    try
        S = load(fullfile(base_dir, file)); 
        fn = fieldnames(S); 
        model = [];
        for k = 1:numel(fn)
            tmp = S.(fn{k}); 
            if isstruct(tmp) && isfield(tmp,'rxns')
                model = tmp; 
                break; 
            end
        end
        if iscell(model)
            model = model{1}; 
        end
        model.rxns = strtrim(model.rxns);
    catch 
        continue;
    end

    for d = 1:nDiets
        dietName = dietNames{d};
        saveFile = fullfile(output_dir, dietName, [modelName '.mat']);
        currentModel = model; 
        
        try
            % Enforce dietary constraints
            if d <= numel(DietTables)
                currentModel = applyDiet(currentModel, DietTables{d}, '(e)', dietName);
            else
                currentModel = applyDiet(currentModel, [], '(e)', dietName);
            end
        catch 
            localResults(d,:) = {modelName, dietName, NaN, NaN, NaN, NaN, 'DIET_FAIL'};
            continue;
        end
        
        WTsol = optimizeCbModel(currentModel, 'max', 'one');
        
        % Assess baseline viability; record extinction failures as NaN
        if isempty(WTsol) || WTsol.stat ~= 1 || WTsol.f < 1e-6
            localResults(d,:) = {modelName, dietName, NaN, NaN, numel(currentModel.rxns), numel(currentModel.rxns), 'WT_Extinction'};
            continue;
        end
        
        WTgrowth = WTsol.f; 
        WTRxns = numel(currentModel.rxns); 
        st = 0.05; 
        
        % Identify Non-Growth Associated Maintenance (NGAM) reaction
        atp_candidates = {'ATPM', 'rxn00062', 'NGAM', 'ATPM_c', 'maintenance', 'DM_atp_c_'};
        found_atpm = intersect(atp_candidates, currentModel.rxns);
        atpRxn = '';
        if ~isempty(found_atpm)
            atpRxn = found_atpm{1}; 
        end
        
        Jmin = [];
        try
            if ~isempty(atpRxn)
                [Jmin_temp, ~] = minReact(currentModel, st, 1e-8, {atpRxn});
            else
                [Jmin_temp, ~] = minReact(currentModel, st, 1e-8); 
            end
            Jmin = double(Jmin_temp);
        catch 
        end
        
        viableFound = false; 
        minimalModel = currentModel; 
        MRMgrowth = WTgrowth; 
        MRMRxns = WTRxns; 
        status = 'Unreducible';
        best_f = -Inf; 
        best_model = []; 
        best_rxns = NaN;

        % Evaluate candidate minimal reactomes — select the SMALLEST viable one
        if ~isempty(Jmin)
            for alt = 1:size(Jmin,1)
                candidateModel = removeRxns(currentModel, currentModel.rxns(~logical(Jmin(alt,:))'));
                if isempty(candidateModel.rxns)
                    continue; 
                end
                
                cSol = optimizeCbModel(candidateModel, 'max', 'one');

                if ~isempty(cSol) && cSol.stat == 1
                    if cSol.f > best_f
                        best_f = cSol.f; 
                        best_model = candidateModel; 
                        best_rxns = numel(candidateModel.rxns);
                    end
                    
                    % Confirm candidate viability against growth cutoff
                    min_threshold = max(st * WTgrowth, 1e-6);
                    
                    if cSol.f >= min_threshold - 1e-9
                        % Accept if viable AND smaller than current best
                        if ~viableFound || numel(candidateModel.rxns) < MRMRxns
                            viableFound = true; 
                            minimalModel = candidateModel;
                            MRMgrowth = cSol.f; 
                            MRMRxns = numel(candidateModel.rxns); 
                            status = 'PASS';
                        end
                        % Do NOT break — continue evaluating remaining candidates
                    end
                end
            end
        end

        if ~viableFound && ~isempty(best_model) && best_f >= 1e-6
            minimalModel = best_model; 
            MRMgrowth = best_f; 
            MRMRxns = best_rxns;
            status = 'SuboptimalMRM';
        end
        
        saveStruct = struct('minimalModel', minimalModel, 'WTgrowth', WTgrowth, 'MRMgrowth', MRMgrowth, ...
            'WTRxns', WTRxns, 'MRMRxns', MRMRxns, 'modelName', modelName, 'dietName', dietName);
        parsave(saveFile, saveStruct);
        
        localResults(d,:) = {modelName, dietName, WTgrowth, MRMgrowth, WTRxns, MRMRxns, status};
    end
    fprintf('%s processed in %.2f mins.\n', modelName, toc(ModelTimer)/60);
    ParpoolResults{m} = localResults;
end

MasterResults = vertcat(ParpoolResults{:});
Elapsed = toc(OverallTimer);

fprintf('\n------------------------------------------------\n');
fprintf('ALL MODELS COMPLETED\nTotal Time : %.2f minutes\n', Elapsed/60);

% Aggregate and export analytical results
if ~isempty(MasterResults)
    ResultTable = cell2table(MasterResults, 'VariableNames', {'Model','Diet','WTGrowth','MRMGrowth','WTRxns','MRMRxns','Status'});
    writetable(ResultTable, fullfile(output_dir, 'minReactModels_Summary.csv'));
    
    models = unique(string(ResultTable.Model), 'stable');
    WideTable = table(models(:), 'VariableNames', {'Model'});
    dietShort = {'HF','MED','UNH','VEG','WES','NOD'};
    
    for d = 1:numel(dietNames)
        wt = nan(height(WideTable),1); 
        mrm = nan(height(WideTable),1);
        for i = 1:height(WideTable)
            idx = strcmp(string(ResultTable.Model), string(WideTable.Model(i))) & strcmp(string(ResultTable.Diet), string(dietNames{d}));
            if any(idx)
                wt(i)  = ResultTable.WTGrowth(idx);
                mrm(i) = ResultTable.MRMGrowth(idx);
            end
        end
        WideTable.([dietShort{d} '_WT'])  = wt;
        WideTable.([dietShort{d} '_MRM']) = mrm;
    end
    wideFile = fullfile(output_dir, 'minReactModels_WideSummary.csv');
    writetable(WideTable, wideFile);
    fprintf('Wide summary saved to %s\n', wideFile);
end

function parsave(fname, dataStruct)
    save(fname, '-struct', 'dataStruct', '-v7.3');
end