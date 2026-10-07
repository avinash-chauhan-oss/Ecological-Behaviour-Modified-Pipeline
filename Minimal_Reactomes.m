% MINIMAL_REACTOMES
%
% Extracts minimal reactomes for 52 AGORA2 gut microbiome models under two
% dietary conditions (High Fiber, Western). For each model-diet pair, three
% reduction scenarios are evaluated:
%   Scenario 1: Strict reduction (GrowthRateCutoff = 1.0).
%   Scenario 2: Relaxed reduction (GrowthRateCutoff = 0.5).
%   Scenario 3: Relaxed reduction preserving exchange and transport reactions.
%
% Non-growing models are supplemented before reduction.

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
nModels = numel(model_files);

% Diet file names must match the actual .txt files in Diets/
dietNames  = {'high_fiber_vmh', 'Western_VMH'};
dietLabels = {'High_Fiber',     'Western'};
nDiets = numel(dietNames);
DietTables = cell(nDiets, 1);

for d = 1:nDiets
    file = fullfile(diet_dir, [dietNames{d} '.txt']);
    T = readtable(file, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);
    rxns = strrep(strtrim(string(T{:,1})), '[u]', '(e)');
    rxns = strrep(rxns, 'EX_H2(e)', 'EX_h2(e)');
    rxns = strrep(rxns, 'EX_glc(e)', 'EX_glc_D(e)');
    DietTables{d}.rxns = rxns;
    DietTables{d}.lbs  = T{:,2};
end

for d = 1:nDiets
    folder = fullfile(output_dir, dietLabels{d});
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end

if isempty(gcp('nocreate'))
    parpool('local');
end

% Output columns: Model, Diet, WTGrowth, WTRxns,
%   MRM_cut1_Growth, MRM_cut1_Rxns,
%   MRM_cut05_Growth, MRM_cut05_Rxns,
%   MRM_ext_Growth, MRM_ext_Rxns,
%   Status_cut1, Status_cut05, Status_ext
nCols = 13;
ParpoolResults = cell(nModels, 1);
fprintf('\nExecuting parallel reduction across %d models...\n', nModels);

parfor m = 1:nModels
    changeCobraSolver('gurobi', 'all', 0);
    ModelTimer = tic;
    file = model_files(m).name;
    [~, modelName] = fileparts(file);
    localResults = cell(nDiets, nCols);

    try
        S = load(fullfile(base_dir, file));
        fn = fieldnames(S);
        model = [];
        for k = 1:numel(fn)
            tmp = S.(fn{k});
            if isstruct(tmp) && isfield(tmp, 'rxns')
                model = tmp;
                break;
            end
        end
        if iscell(model), model = model{1}; end
        model.rxns = strtrim(model.rxns);
    catch ME
        fprintf('Error loading %s: %s\n', file, ME.message);
        continue;
    end

    for d = 1:nDiets
        dietLabel = dietLabels{d};
        saveFile = fullfile(output_dir, dietLabel, [modelName '.mat']);
        currentModel = model;

        try
            currentModel = applyDiet(currentModel, DietTables{d}, '(e)', dietLabel);
        catch ME
            fprintf('Diet error for %s/%s: %s\n', modelName, dietLabel, ME.message);
            localResults(d,:) = {modelName, dietLabel, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'DIET_FAIL', 'DIET_FAIL', 'DIET_FAIL'};
            continue;
        end

        WTsol = optimizeCbModel(currentModel, 'max', 'one');

        status_init = 'PASS';
        if isempty(WTsol) || WTsol.stat ~= 1 || WTsol.f < 1e-6
            tempModel = currentModel;
            ex_idx = find(startsWith(tempModel.rxns, 'EX_'));
            tempModel.lb(ex_idx) = -1000;

            try
                obj_idx = find(tempModel.c == 1);
                if isempty(obj_idx)
                    obj_idx = find(tempModel.c ~= 0);
                end
                if isempty(obj_idx)
                    error('No objective function found in model');
                end
                bioRxn = tempModel.rxns{obj_idx(1)};
                [minMed, ~] = minimalMedium(tempModel, bioRxn, 0.05);

                if ~isempty(minMed)
                    if isstruct(minMed)
                        missing_ex = fieldnames(minMed);
                    elseif iscell(minMed)
                        missing_ex = minMed;
                    else
                        missing_ex = {};
                    end
                    missing_ex = missing_ex(startsWith(missing_ex, 'EX_'));
                    currentModel = changeRxnBounds(currentModel, missing_ex, -10, 'l');
                    WTsol = optimizeCbModel(currentModel, 'max', 'one');
                    status_init = 'Supplemented';
                end
            catch ME
                fprintf('minimalMedium error for %s/%s: %s\n', modelName, dietLabel, ME.message);
            end

            if isempty(WTsol) || WTsol.stat ~= 1 || WTsol.f < 1e-6
                status_init = 'WT_Extinction';
            end
        end

        if strcmp(status_init, 'WT_Extinction')
            localResults(d,:) = {modelName, dietLabel, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'WT_Extinction', 'WT_Extinction', 'WT_Extinction'};
            continue;
        end

        WTgrowth = WTsol.f;
        WTRxns = numel(currentModel.rxns);

        atp_candidates = {'ATPM', 'rxn00062', 'NGAM', 'ATPM_c', 'maintenance', 'DM_atp_c_'};
        found_atpm = intersect(atp_candidates, currentModel.rxns, 'stable');
        eliList = {};
        if ~isempty(found_atpm)
            eliList = found_atpm(:);
        end

        exRxns = currentModel.rxns(startsWith(currentModel.rxns, 'EX_'));
        transIdx = false(numel(currentModel.rxns), 1);
        for r = 1:numel(currentModel.rxns)
            mets_in = find(currentModel.S(:, r) ~= 0);
            comps = regexp(currentModel.mets(mets_in), '\[(\w+)\]|\((\w+)\)', 'tokens');
            comp_ids = {};
            for ci = 1:numel(comps)
                if ~isempty(comps{ci})
                    comp_ids{end+1} = comps{ci}{1}{find(~cellfun(@isempty, comps{ci}{1}), 1)};
                end
            end
            if numel(unique(comp_ids)) > 1
                transIdx(r) = true;
            end
        end
        transRxns = currentModel.rxns(transIdx);
        eliList_extended = unique([eliList; exRxns; transRxns]);

        [mrm_cut1_model, mrm_cut1_growth, mrm_cut1_rxns, mrm_cut1_status, Jmin_cut1] = ...
            get_minimal_reactome(currentModel, 1.0, WTgrowth, eliList);

        [mrm_cut05_model, mrm_cut05_growth, mrm_cut05_rxns, mrm_cut05_status, Jmin_cut05] = ...
            get_minimal_reactome(currentModel, 0.5, WTgrowth, eliList);

        [mrm_ext_model, mrm_ext_growth, mrm_ext_rxns, mrm_ext_status, Jmin_ext] = ...
            get_minimal_reactome(currentModel, 0.5, WTgrowth, eliList_extended);
        if strcmp(status_init, 'Supplemented')
            pass_states = {'PASS', 'SuboptimalMRM', 'Unreducible'};
            if ismember(mrm_cut1_status, pass_states)
                mrm_cut1_status = ['Supplemented_' mrm_cut1_status];
            end
            if ismember(mrm_cut05_status, pass_states)
                mrm_cut05_status = ['Supplemented_' mrm_cut05_status];
            end
            if ismember(mrm_ext_status, pass_states)
                mrm_ext_status = ['Supplemented_' mrm_ext_status];
            end
        end

        saveStruct = struct( ...
            'minimalModel_cut1', mrm_cut1_model, ...
            'minimalModel_cut05', mrm_cut05_model, ...
            'minimalModel_ext', mrm_ext_model, ...
            'WTgrowth', WTgrowth, ...
            'MRMgrowth_cut1', mrm_cut1_growth, ...
            'MRMgrowth_cut05', mrm_cut05_growth, ...
            'MRMgrowth_ext', mrm_ext_growth, ...
            'WTRxns', WTRxns, ...
            'MRMRxns_cut1', mrm_cut1_rxns, ...
            'MRMRxns_cut05', mrm_cut05_rxns, ...
            'MRMRxns_ext', mrm_ext_rxns, ...
            'modelName', modelName, ...
            'dietName', dietLabel, ...
            'status_cut1', mrm_cut1_status, ...
            'status_cut05', mrm_cut05_status, ...
            'status_ext', mrm_ext_status);
        parsave(saveFile, saveStruct);

        localResults(d,:) = {modelName, dietLabel, WTgrowth, WTRxns, ...
            mrm_cut1_growth, mrm_cut1_rxns, mrm_cut05_growth, mrm_cut05_rxns, ...
            mrm_ext_growth, mrm_ext_rxns, mrm_cut1_status, mrm_cut05_status, mrm_ext_status};
    end
    fprintf('%s processed in %.2f mins.\n', modelName, toc(ModelTimer)/60);
    ParpoolResults{m} = localResults;
end

MasterResults = vertcat(ParpoolResults{:});
Elapsed = toc(OverallTimer);

fprintf('\n------------------------------------------------\n');
fprintf('ALL MODELS COMPLETED\nTotal Time : %.2f minutes\n', Elapsed/60);

if ~isempty(MasterResults)
    ResultTable = cell2table(MasterResults, 'VariableNames', ...
        {'Model','Diet','WTGrowth','WTRxns', ...
         'MRM_cut1_Growth','MRM_cut1_Rxns','MRM_cut05_Growth','MRM_cut05_Rxns', ...
         'MRM_ext_Growth','MRM_ext_Rxns','Status_cut1','Status_cut05','Status_ext'});
    writetable(ResultTable, fullfile(output_dir, 'minReactModels_Summary.csv'));

    models = unique(string(ResultTable.Model), 'stable');
    WideTable = table(models(:), 'VariableNames', {'Model'});
    dietShort = {'HF', 'WES'};

    for d = 1:numel(dietLabels)
        wt       = nan(height(WideTable), 1);
        wt_rxns  = nan(height(WideTable), 1);
        mrm_c1   = nan(height(WideTable), 1);
        mrm_c05  = nan(height(WideTable), 1);
        mrm_ext  = nan(height(WideTable), 1);
        rxn_c1   = nan(height(WideTable), 1);
        rxn_c05  = nan(height(WideTable), 1);
        rxn_ext  = nan(height(WideTable), 1);
        for i = 1:height(WideTable)
            idx = strcmp(string(ResultTable.Model), string(WideTable.Model(i))) & ...
                  strcmp(string(ResultTable.Diet), string(dietLabels{d}));
            if any(idx)
                wt(i)      = ResultTable.WTGrowth(idx);
                wt_rxns(i) = ResultTable.WTRxns(idx);
                mrm_c1(i)  = ResultTable.MRM_cut1_Growth(idx);
                mrm_c05(i) = ResultTable.MRM_cut05_Growth(idx);
                mrm_ext(i) = ResultTable.MRM_ext_Growth(idx);
                rxn_c1(i)  = ResultTable.MRM_cut1_Rxns(idx);
                rxn_c05(i) = ResultTable.MRM_cut05_Rxns(idx);
                rxn_ext(i) = ResultTable.MRM_ext_Rxns(idx);
            end
        end
        WideTable.([dietShort{d} '_WT_Growth'])  = wt;
        WideTable.([dietShort{d} '_WT_Rxns'])    = wt_rxns;
        WideTable.([dietShort{d} '_MRM_cut1'])   = mrm_c1;
        WideTable.([dietShort{d} '_Rxns_cut1'])  = rxn_c1;
        WideTable.([dietShort{d} '_MRM_cut05'])  = mrm_c05;
        WideTable.([dietShort{d} '_Rxns_cut05']) = rxn_c05;
        WideTable.([dietShort{d} '_MRM_ext'])    = mrm_ext;
        WideTable.([dietShort{d} '_Rxns_ext'])   = rxn_ext;
    end
    wideFile = fullfile(output_dir, 'minReactModels_WideSummary.csv');
    writetable(WideTable, wideFile);
    fprintf('Wide summary saved to %s\n', wideFile);
end

function [minimalModel, MRMgrowth, MRMRxns, status, Jmin] = get_minimal_reactome(currentModel, st, WTgrowth, eliList)
    Jmin = [];
    minimalModel = currentModel;
    MRMgrowth = WTgrowth;
    MRMRxns = numel(currentModel.rxns);
    status = 'Unreducible';

    try
        if ~isempty(eliList)
            [Jmin_temp, ~] = minReact(currentModel, st, 1e-8, eliList);
        else
            [Jmin_temp, ~] = minReact(currentModel, st, 1e-8);
        end
        Jmin = double(Jmin_temp);
    catch ME
        status = 'MINREACT_FAILED';
        fprintf('minReact error: %s\n', ME.message);
        return;
    end

    viableFound = false;
    best_f = -Inf;
    best_model = [];
    best_rxns = NaN;

    if ~isempty(Jmin)
        for alt = 1:size(Jmin, 1)
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

                min_threshold = max(st * WTgrowth, 1e-6);
                if cSol.f >= min_threshold - 1e-9
                    if ~viableFound || numel(candidateModel.rxns) < MRMRxns
                        viableFound = true;
                        minimalModel = candidateModel;
                        MRMgrowth = cSol.f;
                        MRMRxns = numel(candidateModel.rxns);
                        status = 'PASS';
                    end
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
end

function parsave(fname, dataStruct)
    save(fname, '-struct', 'dataStruct', '-v7.3');
end