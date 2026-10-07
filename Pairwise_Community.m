% PAIRWISE COMMUNITY ANALYSIS (WT & MRM)
%
% Evaluates pairwise community viability for 52 AGORA2 models under two
% dietary conditions (High Fiber, Western). Communities are grown on the
% ORIGINAL (unsupplemented) diets. Three MRM variants are evaluated:
%   1) Strict reduction (GrowthRateCutoff = 1.0).
%   2) Relaxed reduction (GrowthRateCutoff = 0.5).
%   3) Relaxed reduction preserving exchange and transport reactions.
%
% FVA-based growth ranges are computed for each member at the community
% optimum to enable robust interaction classification.

OverallTimer = tic;
fprintf('\n-----------------------------------------------------\n');
fprintf('INITIATING PAIRWISE COMMUNITY ANALYSIS\n');
fprintf('-----------------------------------------------------\n');

base_dir = pwd;
diet_dir = fullfile(base_dir, 'Diets');
mt_dir   = fullfile(base_dir, 'minReactModels');
output_csv = fullfile(base_dir, 'Community_Results_Wide_52Models.csv');

changeCobraSolver('gurobi', 'LP', 0);

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

Pairs = nchoosek(1:nModels, 2);
nPairs = size(Pairs, 1);

dietFileNames = {'high_fiber_vmh', 'Western_VMH'};
dietLabels = {'High_Fiber', 'Western'};
nDiets = numel(dietLabels);
DietTables = cell(nDiets, 1);

for d = 1:nDiets
    file = fullfile(diet_dir, [dietFileNames{d} '.txt']);
    T = readtable(file, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);
    rxns = strrep(strtrim(string(T{:,1})), '[u]', '(e)');
    rxns = strrep(rxns, 'EX_H2(e)', 'EX_h2(e)');
    rxns = strrep(rxns, 'EX_glc(e)', 'EX_glc_D(e)');
    DietTables{d}.rxns = rxns;
    DietTables{d}.lbs  = T{:,2};
end

mrmCuts = {'cut1', 'cut05', 'ext'};
nCuts = numel(mrmCuts);

colsPerDiet = 6 + 6 * nCuts;
headers = {'Community'};
for d = 1:nDiets
    dn = dietLabels{d};
    headers{end+1} = sprintf('WT_Rxns_%s', dn);
    headers{end+1} = sprintf('WT_Growth_%s', dn);
    headers{end+1} = sprintf('WT_M1_Min_%s', dn);
    headers{end+1} = sprintf('WT_M1_Max_%s', dn);
    headers{end+1} = sprintf('WT_M2_Min_%s', dn);
    headers{end+1} = sprintf('WT_M2_Max_%s', dn);
    for c = 1:nCuts
        cn = mrmCuts{c};
        headers{end+1} = sprintf('MRM_%s_Rxns_%s', cn, dn);
        headers{end+1} = sprintf('MRM_%s_Growth_%s', cn, dn);
        headers{end+1} = sprintf('MRM_%s_M1_Min_%s', cn, dn);
        headers{end+1} = sprintf('MRM_%s_M1_Max_%s', cn, dn);
        headers{end+1} = sprintf('MRM_%s_M2_Min_%s', cn, dn);
        headers{end+1} = sprintf('MRM_%s_M2_Max_%s', cn, dn);
    end
end

if isempty(gcp('nocreate'))
    parpool('local');
end

WorkerResults = cell(nPairs, 1);
fprintf('Evaluating %d pairwise communities...\n', nPairs);

parfor p = 1:nPairs
    changeCobraSolver('gurobi', 'LP', 0);
    idx1 = Pairs(p,1);
    idx2 = Pairs(p,2);
    [~,name1] = fileparts(model_files(idx1).name);
    [~,name2] = fileparts(model_files(idx2).name);
    communityName = [name1 '_AND_' name2];

    LocalRow = cell(1, length(headers));
    LocalRow{1} = communityName;

    try
        S1 = load(fullfile(base_dir, model_files(idx1).name));
        fn1 = fieldnames(S1); model1_WT = [];
        for k=1:numel(fn1), tmp=S1.(fn1{k}); if isstruct(tmp) && isfield(tmp,'rxns') && isfield(tmp,'mets'), model1_WT=tmp; break; end; end

        S2 = load(fullfile(base_dir, model_files(idx2).name));
        fn2 = fieldnames(S2); model2_WT = [];
        for k=1:numel(fn2), tmp=S2.(fn2{k}); if isstruct(tmp) && isfield(tmp,'rxns') && isfield(tmp,'mets'), model2_WT=tmp; break; end; end

        if iscell(model1_WT), model1_WT = model1_WT{1}; end
        if iscell(model2_WT), model2_WT = model2_WT{1}; end

        model1_WT.rxns = strtrim(model1_WT.rxns);
        model2_WT.rxns = strtrim(model2_WT.rxns);

        if ~isfield(model1_WT, 'csense')
            model1_WT.csense = repmat('E', size(model1_WT.S, 1), 1);
        end
        if ~isfield(model2_WT, 'csense')
            model2_WT.csense = repmat('E', size(model2_WT.S, 1), 1);
        end

        ex1 = startsWith(model1_WT.rxns, 'EX_');
        model1_WT.lb(ex1) = -1000; model1_WT.ub(ex1) = 1000;
        ex2 = startsWith(model2_WT.rxns, 'EX_');
        model2_WT.lb(ex2) = -1000; model2_WT.ub(ex2) = 1000;

        bio1_WT = model1_WT.rxns(model1_WT.c ~= 0);
        bio2_WT = model2_WT.rxns(model2_WT.c ~= 0);

        if isempty(bio1_WT) || isempty(bio2_WT)
            WorkerResults{p} = LocalRow;
            continue;
        end

        removeFields = {'osense','osenseStr','metCharge','metCharges','metFormulas','metNames','metCHEBIID','metHMDBID','metKEGGID','metPubChemID','metSmile','metSmiles','metInChIString','metInchiString','metSEEDID','rules','grRules','rxnGeneMat','genes','geneNames','subSystems','proteinClasses','comments','rxnConfidenceScores','citations','ecNumbers'};
        for f = 1:numel(removeFields)
            if isfield(model1_WT,removeFields{f}), model1_WT = rmfield(model1_WT,removeFields{f}); end
            if isfield(model2_WT,removeFields{f}), model2_WT = rmfield(model2_WT,removeFields{f}); end
        end

        WT_Community = createMultipleSpeciesModel({model1_WT; model2_WT}, {bio1_WT{1}; bio2_WT{1}}, 'mergeGenesFlag', false);

        WT_Community.c(:) = 0;
        bioIdx1_WT = find(strcmp(WT_Community.rxns, ['model1_' bio1_WT{1}]), 1);
        bioIdx2_WT = find(strcmp(WT_Community.rxns, ['model2_' bio2_WT{1}]), 1);

        if isempty(bioIdx1_WT), bioIdx1_WT = find(endsWith(WT_Community.rxns, ['model1_' bio1_WT{1}]), 1); end
        if isempty(bioIdx2_WT), bioIdx2_WT = find(endsWith(WT_Community.rxns, ['model2_' bio2_WT{1}]), 1); end

        if isempty(bioIdx1_WT) || isempty(bioIdx2_WT)
            WorkerResults{p} = LocalRow;
            continue;
        end

        WT_Community.c(bioIdx1_WT) = 1;
        WT_Community.c(bioIdx2_WT) = 1;
        WT_Community.lb(bioIdx1_WT) = 0.001;
        WT_Community.lb(bioIdx2_WT) = 0.001;
        WT_Community.rxns = strtrim(WT_Community.rxns);

        exBioIdxWT = contains(WT_Community.rxns, 'EX_') & contains(WT_Community.rxns, 'biomass') & contains(WT_Community.rxns, '[u]');
        WT_Community.lb(exBioIdxWT) = 0;
        WT_Community.ub(exBioIdxWT) = 0;

    catch ME
        WorkerResults{p} = LocalRow;
        continue;
    end

    for d = 1:nDiets
        col = 2 + (d-1) * colsPerDiet;
        dietLabel = dietLabels{d};

        % --- WT COMMUNITY ---
        WT_Rxns = NaN; WT_G = NaN;
        WT_M1_Min = NaN; WT_M1_Max = NaN;
        WT_M2_Min = NaN; WT_M2_Max = NaN;

        try
            WT_Current = WT_Community;
            WT_Current = applyDiet(WT_Current, DietTables{d}, '[u]', dietLabel);

            solWT = optimizeCbModel(WT_Current, 'max', 'one');
            WT_Rxns = numel(WT_Current.rxns);

            if ~isempty(solWT) && solWT.stat == 1 && solWT.f >= 0.001
                WT_G = solWT.f;

                [WT_M1_Min, WT_M1_Max, WT_M2_Min, WT_M2_Max] = ...
                    community_fva(WT_Current, bioIdx1_WT, bioIdx2_WT, solWT.f);
            end
        catch
        end

        LocalRow{col}   = WT_Rxns;
        LocalRow{col+1} = WT_G;
        LocalRow{col+2} = WT_M1_Min;
        LocalRow{col+3} = WT_M1_Max;
        LocalRow{col+4} = WT_M2_Min;
        LocalRow{col+5} = WT_M2_Max;

        % --- MRM COMMUNITIES (one block per cut) ---
        for c = 1:nCuts
            cutName = mrmCuts{c};
            mCol = col + 6 + (c-1) * 6;

            MT_Rxns = NaN; MT_G = NaN;
            MT_M1_Min = NaN; MT_M1_Max = NaN;
            MT_M2_Min = NaN; MT_M2_Max = NaN;

            try
                mtFile1 = fullfile(mt_dir, dietLabel, [name1 '.mat']);
                mtFile2 = fullfile(mt_dir, dietLabel, [name2 '.mat']);

                mt_valid = exist(mtFile1, 'file') && exist(mtFile2, 'file');

                if mt_valid
                    mrmField = ['minimalModel_' cutName];
                    growthField = ['MRMgrowth_' cutName];

                    MS1 = load(mtFile1);
                    MS2 = load(mtFile2);

                    if isfield(MS1, growthField) && MS1.(growthField) >= 1e-6 && isfield(MS1, mrmField)
                        model1_MT = MS1.(mrmField);
                    else
                        mt_valid = false;
                    end

                    if mt_valid && isfield(MS2, growthField) && MS2.(growthField) >= 1e-6 && isfield(MS2, mrmField)
                        model2_MT = MS2.(mrmField);
                    else
                        mt_valid = false;
                    end
                end

                if mt_valid
                    while iscell(model1_MT), model1_MT = model1_MT{1}; end
                    while iscell(model2_MT), model2_MT = model2_MT{1}; end

                    model1_MT.rxns = strtrim(model1_MT.rxns);
                    model2_MT.rxns = strtrim(model2_MT.rxns);

                    if ~isfield(model1_MT, 'csense')
                        model1_MT.csense = repmat('E', size(model1_MT.S, 1), 1);
                    end
                    if ~isfield(model2_MT, 'csense')
                        model2_MT.csense = repmat('E', size(model2_MT.S, 1), 1);
                    end

                    ex1_MT = startsWith(model1_MT.rxns, 'EX_');
                    model1_MT.lb(ex1_MT) = -1000; model1_MT.ub(ex1_MT) = 1000;
                    ex2_MT = startsWith(model2_MT.rxns, 'EX_');
                    model2_MT.lb(ex2_MT) = -1000; model2_MT.ub(ex2_MT) = 1000;

                    bio1_MT = find(model1_MT.c ~= 0, 1);
                    bio2_MT = find(model2_MT.c ~= 0, 1);

                    if ~isempty(bio1_MT) && ~isempty(bio2_MT)
                        bioRxn1 = model1_MT.rxns{bio1_MT};
                        bioRxn2 = model2_MT.rxns{bio2_MT};

                        coreFields = {'rxns','mets','S','lb','ub','c','b','csense','description'};
                        allF = fieldnames(model1_MT);
                        for k = 1:numel(allF)
                            if ~ismember(allF{k}, coreFields)
                                if isfield(model1_MT, allF{k}), model1_MT = rmfield(model1_MT, allF{k}); end
                                if isfield(model2_MT, allF{k}), model2_MT = rmfield(model2_MT, allF{k}); end
                            end
                        end

                        MT_Community = createMultipleSpeciesModel({model1_MT; model2_MT}, {bioRxn1; bioRxn2}, 'mergeGenesFlag', false);

                        MT_Community.c(:) = 0;
                        bioIdx1_MT = find(strcmp(MT_Community.rxns, ['model1_' bioRxn1]), 1);
                        bioIdx2_MT = find(strcmp(MT_Community.rxns, ['model2_' bioRxn2]), 1);

                        if isempty(bioIdx1_MT), bioIdx1_MT = find(endsWith(MT_Community.rxns, ['model1_' bioRxn1]), 1); end
                        if isempty(bioIdx2_MT), bioIdx2_MT = find(endsWith(MT_Community.rxns, ['model2_' bioRxn2]), 1); end

                        if ~isempty(bioIdx1_MT) && ~isempty(bioIdx2_MT)
                            MT_Community.c(bioIdx1_MT) = 1;
                            MT_Community.c(bioIdx2_MT) = 1;
                            MT_Community.lb(bioIdx1_MT) = 0.001;
                            MT_Community.lb(bioIdx2_MT) = 0.001;
                            MT_Community.rxns = strtrim(MT_Community.rxns);

                            exBioIdxMT = contains(MT_Community.rxns, 'EX_') & contains(MT_Community.rxns, 'biomass') & contains(MT_Community.rxns, '[u]');
                            MT_Community.lb(exBioIdxMT) = 0;
                            MT_Community.ub(exBioIdxMT) = 0;

                            MT_Current = MT_Community;
                            MT_Current = applyDiet(MT_Current, DietTables{d}, '[u]', dietLabel);

                            solMT = optimizeCbModel(MT_Current, 'max', 'one');
                            MT_Rxns = numel(MT_Current.rxns);

                            if ~isempty(solMT) && solMT.stat == 1 && solMT.f >= 0.001
                                MT_G = solMT.f;

                                [MT_M1_Min, MT_M1_Max, MT_M2_Min, MT_M2_Max] = ...
                                    community_fva(MT_Current, bioIdx1_MT, bioIdx2_MT, solMT.f);
                            end
                        end
                    end
                end
            catch
            end

            LocalRow{mCol}   = MT_Rxns;
            LocalRow{mCol+1} = MT_G;
            LocalRow{mCol+2} = MT_M1_Min;
            LocalRow{mCol+3} = MT_M1_Max;
            LocalRow{mCol+4} = MT_M2_Min;
            LocalRow{mCol+5} = MT_M2_Max;
        end
    end

    WorkerResults{p} = LocalRow;
end

Results = vertcat(WorkerResults{:});
ResultTable = cell2table(Results, 'VariableNames', headers);
writetable(ResultTable, output_csv);

fprintf('\n-----------------------------------------------------\n');
fprintf('ANALYSIS COMPLETE\nTotal Time : %.2f minutes\n', toc(OverallTimer)/60);
fprintf('Data saved to: %s\n', output_csv);
fprintf('-----------------------------------------------------\n');


function [M1_Min, M1_Max, M2_Min, M2_Max] = community_fva(model, bioIdx1, bioIdx2, totalGrowth)
    M1_Min = NaN; M1_Max = NaN;
    M2_Min = NaN; M2_Max = NaN;

    [nMets, ~] = size(model.S);
    model.S(nMets+1, bioIdx1) = 1;
    model.S(nMets+1, bioIdx2) = 1;
    model.b(nMets+1) = totalGrowth;
    model.mets{nMets+1} = 'CommunityBiomassFix';
    if isfield(model, 'csense')
        model.csense = [model.csense(:); 'E'];
    end

    model.c(:) = 0;
    model.c(bioIdx1) = 1;
    solMin = optimizeCbModel(model, 'min');
    if ~isempty(solMin) && solMin.stat == 1, M1_Min = solMin.f; end
    solMax = optimizeCbModel(model, 'max');
    if ~isempty(solMax) && solMax.stat == 1, M1_Max = solMax.f; end

    model.c(:) = 0;
    model.c(bioIdx2) = 1;
    solMin = optimizeCbModel(model, 'min');
    if ~isempty(solMin) && solMin.stat == 1, M2_Min = solMin.f; end
    solMax = optimizeCbModel(model, 'max');
    if ~isempty(solMax) && solMax.stat == 1, M2_Max = solMax.f; end
end