% =========================================================================
% PAIRWISE COMMUNITY ANALYSIS (WT & MRM) 
% 
% DESCRIPTION:
% Evaluates pairwise community viability across dietary conditions.
% Maximizes the total overall community growth rate (Model1 + Model2) 
% Enforces a biological 
% coexistence floor (1e-6) to prevent mathematical solver artifacts from 
% artificially starving one community member. 
% =========================================================================

OverallTimer = tic;
fprintf('\n-----------------------------------------------------\n');
fprintf('INITIATING PAIRWISE COMMUNITY ANALYSIS\n');
fprintf('-----------------------------------------------------\n');

base_dir = pwd;
diet_dir = fullfile(base_dir, 'Diets');
mt_dir   = fullfile(base_dir, 'minReactModels'); % Pulling from initial phase directory
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
%
nModels = numel(model_files);

Pairs = nchoosek(1:nModels, 2);
nPairs = size(Pairs, 1);

dietNames = {'high_fiber_vmh', 'Mediterranian', 'Unhealthy', 'vegetarian_diet', 'Western_VMH', 'No_Diet'};
nDiets = numel(dietNames);
DietTables = cell(nDiets, 1);

% Load Diet Constraints
for d = 1:nDiets
    if strcmp(dietNames{d}, 'No_Diet')
        DietTables{d}.rxns = [];
        DietTables{d}.lbs  = [];
        continue;
    end
    
    file = fullfile(diet_dir, [dietNames{d} '.txt']);
    T = readtable(file, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);
    rxns = strrep(strtrim(string(T{:,1})), '[u]', '(e)');
    DietTables{d}.rxns = rxns;
    DietTables{d}.lbs  = T{:,2};
end

% Construct Output Headers
headers = {'Community'};
for d = 1:nDiets
    dn = dietNames{d};
    headers{end+1} = sprintf('WT_Rxns_%s', dn);
    headers{end+1} = sprintf('WT_Growth_%s', dn);
    headers{end+1} = sprintf('WT_M1_Growth_%s', dn);
    headers{end+1} = sprintf('WT_M2_Growth_%s', dn);
    headers{end+1} = sprintf('MRM_Rxns_%s', dn);
    headers{end+1} = sprintf('MRM_Growth_%s', dn);
    headers{end+1} = sprintf('MRM_M1_Growth_%s', dn);
    headers{end+1} = sprintf('MRM_M2_Growth_%s', dn);
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
        % --- LOAD WILD-TYPE MODELS ---
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
        
        bio1_WT = model1_WT.rxns(model1_WT.c ~= 0);
        bio2_WT = model2_WT.rxns(model2_WT.c ~= 0);

        if isempty(bio1_WT) || isempty(bio2_WT)
            WorkerResults{p} = LocalRow;
            continue;
        end
        
        removeFields = {'C','ctrs','d','dsense','osense','osenseStr','metCharge','metCharges','metFormulas','metNames','metCHEBIID','metHMDBID','metKEGGID','metPubChemID','metSmile','metSmiles','metInChIString','metInchiString','metSEEDID','rules','grRules','rxnGeneMat','genes','geneNames','subSystems','proteinClasses','comments','rxnConfidenceScores','citations','ecNumbers'};
        for f = 1:numel(removeFields)
            if isfield(model1_WT,removeFields{f}), model1_WT = rmfield(model1_WT,removeFields{f}); end
            if isfield(model2_WT,removeFields{f}), model2_WT = rmfield(model2_WT,removeFields{f}); end
        end
        
        % Build WT Community Model
        WT_Community = createMultipleSpeciesModel({model1_WT; model2_WT}, {bio1_WT{1}; bio2_WT{1}}, 'mergeGenesFlag', false);
        
        WT_Community.c(:) = 0;
        bioIdx1_WT = find(strcmp(WT_Community.rxns,['model1_' bio1_WT{1}]), 1);
        bioIdx2_WT = find(strcmp(WT_Community.rxns,['model2_' bio2_WT{1}]), 1);
        
        if isempty(bioIdx1_WT), bioIdx1_WT = find(endsWith(WT_Community.rxns,['model1_' bio1_WT{1}]), 1); end
        if isempty(bioIdx2_WT), bioIdx2_WT = find(endsWith(WT_Community.rxns,['model2_' bio2_WT{1}]), 1); end

        if isempty(bioIdx1_WT) || isempty(bioIdx2_WT)
            WorkerResults{p} = LocalRow;
            continue;
        end

        % Maximize Community Growth (Sum of both models)
        WT_Community.c(bioIdx1_WT) = 1;
        WT_Community.c(bioIdx2_WT) = 1;
        
        % BIOLOGICAL COEXISTENCE & RATIO COUPLING
        % Allow flexibility but prevent the "Winner-Takes-All" LP artifact
        % by constraining the growth ratio to max 10:1 in either direction.
        WT_Community.lb(bioIdx1_WT) = 1e-6;
        WT_Community.lb(bioIdx2_WT) = 1e-6;
        
        [m_wt, ~] = size(WT_Community.S);
        WT_Community.S(m_wt+1, bioIdx1_WT) = 1;
        WT_Community.S(m_wt+1, bioIdx2_WT) = -10;
        WT_Community.b(m_wt+1) = 0;
        WT_Community.mets{m_wt+1} = 'CouplingConstraint_1to2';
        if isfield(WT_Community, 'csense'), WT_Community.csense(m_wt+1) = 'L'; end
        
        WT_Community.S(m_wt+2, bioIdx1_WT) = -10;
        WT_Community.S(m_wt+2, bioIdx2_WT) = 1;
        WT_Community.b(m_wt+2) = 0;
        WT_Community.mets{m_wt+2} = 'CouplingConstraint_2to1';
        if isfield(WT_Community, 'csense'), WT_Community.csense(m_wt+2) = 'L'; end
        
        WT_Community.rxns = strtrim(WT_Community.rxns);
       
    catch ME
        WorkerResults{p} = LocalRow;
        continue;
    end

    for d = 1:nDiets
        col = 2 + (d-1)*8; 
        dietName = dietNames{d};
        
        WT_Rxns = NaN; WT_G = NaN; WT_M1 = NaN; WT_M2 = NaN;
        MT_Rxns = NaN; MT_G = NaN; MT_M1 = NaN; MT_M2 = NaN;
        
        % --- EVALUATE WILD-TYPE COMMUNITY ---
        try
            WT_Current = WT_Community;
            % Apply diet using [u] for the universal community compartment
            WT_Current = applyDiet(WT_Current, DietTables{d}, '[u]', dietName);
            
            solWT = optimizeCbModel(WT_Current, 'max', 'one');
            WT_Rxns = numel(WT_Current.rxns);
            
            if ~isempty(solWT) && solWT.stat == 1 && solWT.f >= 1e-6
                WT_G  = solWT.f;
                WT_M1 = solWT.x(bioIdx1_WT);
                WT_M2 = solWT.x(bioIdx2_WT);
            end
        catch ME
        end
        
        % --- EVALUATE MINIMAL REACTOME (MRM) COMMUNITY ---
        try
            % Load MT models directly from their specific diet folders
            mtFile1 = fullfile(mt_dir, dietName, [name1 '.mat']);
            mtFile2 = fullfile(mt_dir, dietName, [name2 '.mat']);
            
            mt_valid = exist(mtFile1,'file') && exist(mtFile2,'file');
            
            if mt_valid
                S1 = load(mtFile1);
                S2 = load(mtFile2);
                
                % STRICT VALIDATION: Ensure the MT models actually survived the minimal reactome phase
                if isfield(S1, 'MRMgrowth') && S1.MRMgrowth >= 1e-6 && isfield(S1,'minimalModel')
                    model1_MT = S1.minimalModel; 
                else
                    mt_valid = false; 
                end
                
                if isfield(S2, 'MRMgrowth') && S2.MRMgrowth >= 1e-6 && isfield(S2,'minimalModel')
                    model2_MT = S2.minimalModel; 
                else
                    mt_valid = false; 
                end
                
                if mt_valid
                    while iscell(model1_MT), model1_MT = model1_MT{1}; end
                    while iscell(model2_MT), model2_MT = model2_MT{1}; end
                    
                    model1_MT.rxns = strtrim(model1_MT.rxns);
                    model2_MT.rxns = strtrim(model2_MT.rxns);
                    
                    bio1_MT = find(model1_MT.c~=0, 1);
                    bio2_MT = find(model2_MT.c~=0, 1);

                    if isempty(bio1_MT) || isempty(bio2_MT)
                        mt_valid = false;
                    else
                        bioRxn1 = model1_MT.rxns{bio1_MT};
                        bioRxn2 = model2_MT.rxns{bio2_MT};
                    end
                end
            end
            
            if mt_valid
                coreFields = {'rxns','mets','S','lb','ub','c','b','csense','description'};
                allF = fieldnames(model1_MT);
                for k = 1:numel(allF)
                    if ~ismember(allF{k},coreFields)
                        if isfield(model1_MT,allF{k}), model1_MT = rmfield(model1_MT,allF{k}); end
                        if isfield(model2_MT,allF{k}), model2_MT = rmfield(model2_MT,allF{k}); end
                    end
                end
                
                % Build MRM Community Model
                MT_Community = createMultipleSpeciesModel({model1_MT; model2_MT}, {bioRxn1; bioRxn2}, 'mergeGenesFlag', false);
                
                MT_Community.c(:) = 0;
                bioIdx1_MT = find(strcmp(MT_Community.rxns,['model1_' bioRxn1]), 1);
                bioIdx2_MT = find(strcmp(MT_Community.rxns,['model2_' bioRxn2]), 1);

                if isempty(bioIdx1_MT) || isempty(bioIdx2_MT)
                    bioIdx1_MT = find(endsWith(MT_Community.rxns,['model1_' bioRxn1]), 1);
                    bioIdx2_MT = find(endsWith(MT_Community.rxns,['model2_' bioRxn2]), 1);
                end

                if ~isempty(bioIdx1_MT) && ~isempty(bioIdx2_MT)
                    
                    % Maximize Community Growth (Sum of both models)
                    MT_Community.c(bioIdx1_MT) = 1;
                    MT_Community.c(bioIdx2_MT) = 1;
                    
                    % BIOLOGICAL COEXISTENCE FLOOR & RATIO COUPLING for MRM
                    MT_Community.lb(bioIdx1_MT) = 1e-6;
                    MT_Community.lb(bioIdx2_MT) = 1e-6;
                    
                    [m_mt, ~] = size(MT_Community.S);
                    MT_Community.S(m_mt+1, bioIdx1_MT) = 1;
                    MT_Community.S(m_mt+1, bioIdx2_MT) = -10;
                    MT_Community.b(m_mt+1) = 0;
                    MT_Community.mets{m_mt+1} = 'CouplingConstraint_1to2';
                    if isfield(MT_Community, 'csense'), MT_Community.csense(m_mt+1) = 'L'; end
                    
                    MT_Community.S(m_mt+2, bioIdx1_MT) = -10;
                    MT_Community.S(m_mt+2, bioIdx2_MT) = 1;
                    MT_Community.b(m_mt+2) = 0;
                    MT_Community.mets{m_mt+2} = 'CouplingConstraint_2to1';
                    if isfield(MT_Community, 'csense'), MT_Community.csense(m_mt+2) = 'L'; end
                    
                    MT_Community.rxns = strtrim(MT_Community.rxns);
                    
                    MT_Current = MT_Community;
                    MT_Current = applyDiet(MT_Current, DietTables{d}, '[u]', dietName);
                    
                    solMT = optimizeCbModel(MT_Current, 'max', 'one');
                    MT_Rxns = numel(MT_Current.rxns);
                    
                    if ~isempty(solMT) && solMT.stat == 1 && solMT.f >= 1e-6
                        MT_G  = solMT.f;
                        MT_M1 = solMT.x(bioIdx1_MT);
                        MT_M2 = solMT.x(bioIdx2_MT);
                    end
                end
            end
        catch ME
        end
        
        LocalRow{col}   = WT_Rxns;
        LocalRow{col+1} = WT_G;
        LocalRow{col+2} = WT_M1;
        LocalRow{col+3} = WT_M2;
        LocalRow{col+4} = MT_Rxns;
        LocalRow{col+5} = MT_G;
        LocalRow{col+6} = MT_M1;
        LocalRow{col+7} = MT_M2;
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