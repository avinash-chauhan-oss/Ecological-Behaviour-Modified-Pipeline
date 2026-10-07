% Plot6.m - Phylum-Pair Interaction Composition (Modified Pipeline)
% Shows composition of interactions by phylum pair for WT and MRM scenarios.

%% Configuration
commFile = 'Community_Results_Wide_52Models.csv';
monoFile = 'minReactModels/minReactModels_WideSummary.csv';

if ~isfile(commFile) || ~isfile(monoFile)
    error('Required CSV files not found. Run the pipeline first.');
end

commData = readtable(commFile, 'PreserveVariableNames', true);
monoData = readtable(monoFile, 'PreserveVariableNames', true);

dietLabels = {'High_Fiber', 'Western'};
dietPrefixes = {'HF', 'WES'};
dietTitles = {'High Fiber Diet', 'Western Diet'};
mrmCuts = {'cut1', 'cut05', 'ext'};
mrmTitles = {'MRM (Cutoff 1.0)', 'MRM (Cutoff 0.5)', 'MRM (Cutoff 0.5 + EXT)'};

taxMap = getPhylumMapping();
phyla = unique(values(taxMap));
numPhyla = length(phyla);
phylumIdx = containers.Map(phyla, 1:numPhyla);

threshold = 0.1;

% 1: Mutualism, 2: Commensalism, 3: Amensalism, 4: Parasitism, 5: Competition, 6: Neutral
colorMap = [0 1 0; 0.5 1 0.5; 1 0.5 0.5; 1 0.5 0; 1 0 0; 0.8 0.8 0.8];

%% Processing
for d = 1:length(dietLabels)
    diet = dietLabels{d};
    pref = dietPrefixes{d};
    
    intData = cell(4, 1);
    for k = 1:4
        intData{k} = zeros(numPhyla, numPhyla, 6);
    end
    
    for i = 1:height(commData)
        names = strsplit(commData.Community{i}, '_AND_');
        sp1 = names{1};
        sp2 = names{2};
        
        p1 = taxMap(sp1);
        p2 = taxMap(sp2);
        idx1 = phylumIdx(p1);
        idx2 = phylumIdx(p2);
        
        if idx1 > idx2
            temp = idx1; idx1 = idx2; idx2 = temp;
        end
        
        wt_m1_mono = monoData{strcmp(monoData.Model, sp1), [pref '_WT_Growth']};
        wt_m2_mono = monoData{strcmp(monoData.Model, sp2), [pref '_WT_Growth']};
        wt_m1_comm = commData{i, ['WT_M1_Max_' diet]};
        wt_m2_comm = commData{i, ['WT_M2_Max_' diet]};
        
        alpha1WT = calcAlpha(wt_m1_comm, wt_m1_mono);
        alpha2WT = calcAlpha(wt_m2_comm, wt_m2_mono);
        typeWT = classifyInteraction(alpha1WT, alpha2WT, threshold);
        intData{1}(idx1, idx2, typeWT) = intData{1}(idx1, idx2, typeWT) + 1;
        
        for c = 1:length(mrmCuts)
            cut = mrmCuts{c};
            mrm_m1_mono = monoData{strcmp(monoData.Model, sp1), [pref '_MRM_' cut]};
            mrm_m2_mono = monoData{strcmp(monoData.Model, sp2), [pref '_MRM_' cut]};
            mrm_m1_comm = commData{i, ['MRM_' cut '_M1_Max_' diet]};
            mrm_m2_comm = commData{i, ['MRM_' cut '_M2_Max_' diet]};
            
            alpha1MRM = calcAlpha(mrm_m1_comm, mrm_m1_mono);
            alpha2MRM = calcAlpha(mrm_m2_comm, mrm_m2_mono);
            typeMRM = classifyInteraction(alpha1MRM, alpha2MRM, threshold);
            intData{c+1}(idx1, idx2, typeMRM) = intData{c+1}(idx1, idx2, typeMRM) + 1;
        end
    end
    
    figure('Name', sprintf('Phylum Interactions - %s', dietTitles{d}), 'Position', [50, 50, 1400, 400]);
    titles = [{'WT'}, mrmTitles];
    
    for k = 1:4
        subplot(1, 4, k);
        drawPhylumGrid(intData{k}, phyla, colorMap);
        title(sprintf('%s - %s', dietTitles{d}, titles{k}), 'Interpreter', 'none');
    end
end

function drawPhylumGrid(dataMatrix, labels, colors)
    n = length(labels);
    hold on;
    maxTotal = max(sum(dataMatrix, 3), [], 'all');
    if maxTotal == 0, maxTotal = 1; end
    
    for i = 1:n
        for j = i:n
            counts = squeeze(dataMatrix(i, j, :));
            total = sum(counts);
            if total > 0
                x = j; y = n - i + 1;
                drawScaledPie(x, y, counts, colors, total, maxTotal, 0.45);
            end
        end
    end
    
    xlim([0.5 n+0.5]); ylim([0.5 n+0.5]);
    xticks(1:n); xticklabels(labels); xtickangle(90);
    yticks(1:n); yticklabels(flipud(labels'));
    set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 7);
    axis square; box on; hold off;
end

function drawScaledPie(x, y, counts, colors, total, maxTotal, maxRadius)
    radius = maxRadius * sqrt(total / maxTotal);
    if total == 0, return; end
    
    angles = 2 * pi * (counts / total);
    startAngle = 0;
    
    for k = 1:length(counts)
        if counts(k) > 0
            endAngle = startAngle + angles(k);
            theta = linspace(startAngle, endAngle, 50);
            xf = x + radius * [0, cos(theta), 0];
            yf = y + radius * [0, sin(theta), 0];
            fill(xf, yf, colors(k,:), 'EdgeColor', 'none');
            startAngle = endAngle;
        end
    end
    theta = linspace(0, 2*pi, 100);
    plot(x + radius*cos(theta), y + radius*sin(theta), 'k', 'LineWidth', 0.5);
end

function taxMap = getPhylumMapping()
    taxMap = containers.Map();
    taxMap('Acidaminococcus_fermentans_DSM_20731') = 'Firmicutes';
    taxMap('Akkermansia_muciniphila_ATCC_BAA_835') = 'Verrucomicrobia';
    taxMap('Alistipes_finegoldii_DSM_17242') = 'Bacteroidetes';
    taxMap('Alistipes_putredinis_DSM_17216') = 'Bacteroidetes';
    taxMap('Alistipes_shahii_WAL_8301') = 'Bacteroidetes';
    taxMap('Bacteroides_caccae_ATCC_43185') = 'Bacteroidetes';
    taxMap('Bacteroides_coprocola_M16_DSM_17136') = 'Bacteroidetes';
    taxMap('Bacteroides_eggerthii_DSM_20697') = 'Bacteroidetes';
    taxMap('Bacteroides_fragilis_NCTC_9343') = 'Bacteroidetes';
    taxMap('Bacteroides_intestinalis_341_DSM_17393') = 'Bacteroidetes';
    taxMap('Bacteroides_stercoris_ATCC_43183') = 'Bacteroidetes';
    taxMap('Bacteroides_thetaiotaomicron_VPI_5482') = 'Bacteroidetes';
    taxMap('Bacteroides_uniformis_ATCC_8492') = 'Bacteroidetes';
    taxMap('Bifidobacterium_adolescentis_ATCC_15703') = 'Actinobacteria';
    taxMap('Bifidobacterium_bifidum_PRL2010') = 'Actinobacteria';
    taxMap('Bifidobacterium_longum_NCC2705') = 'Actinobacteria';
    taxMap('Bilophila_wadsworthia_3_1_6') = 'Proteobacteria';
    taxMap('Blautia_wexlerae_DSM_19850') = 'Firmicutes';
    taxMap('Citrobacter_amalonaticus_Y19') = 'Proteobacteria';
    taxMap('Clostridium_clostridioforme_CM201') = 'Firmicutes';
    taxMap('Collinsella_tanakaei_YIT_12063') = 'Actinobacteria';
    taxMap('Coprococcus_catus_GD_7') = 'Firmicutes';
    taxMap('Desulfovibrio_desulfuricans_subsp_desulfuricans_DSM_642') = 'Proteobacteria';
    taxMap('Desulfovibrio_piger_ATCC_29098') = 'Proteobacteria';
    taxMap('Dialister_invisus_DSM_15470') = 'Firmicutes';
    taxMap('Dorea_longicatena_DSM_13814') = 'Firmicutes';
    taxMap('Enterobacter_cloacae_EcWSU1') = 'Proteobacteria';
    taxMap('Enterococcus_faecium_TX1330') = 'Firmicutes';
    taxMap('Escherichia_coli_str_K_12_substr_MG1655') = 'Proteobacteria';
    taxMap('Eubacterium_rectale_M104_1') = 'Firmicutes';
    taxMap('Faecalibacterium_prausnitzii_L2_6') = 'Firmicutes';
    taxMap('Fusobacterium_varium_ATCC_27725') = 'Fusobacteria';
    taxMap('Haemophilus_parainfluenzae_T3T1') = 'Proteobacteria';
    taxMap('Helicobacter_pylori_26695') = 'Proteobacteria';
    taxMap('Klebsiella_pneumoniae_pneumoniae_MGH78578') = 'Proteobacteria';
    taxMap('Lactobacillus_mucosae_LM1') = 'Firmicutes';
    taxMap('Lactobacillus_reuteri_SD2112_ATCC_55730') = 'Firmicutes';
    taxMap('Megasphaera_elsdenii_DSM_20460') = 'Firmicutes';
    taxMap('Parabacteroides_distasonis_ATCC_8503') = 'Bacteroidetes';
    taxMap('Parabacteroides_johnsonii_DSM_18315') = 'Bacteroidetes';
    taxMap('Phascolarctobacterium_succinatutens_YIT_12067') = 'Firmicutes';
    taxMap('Prevotella_ruminicola_23') = 'Bacteroidetes';
    taxMap('Pseudoflavonifractor_capillosus_strain_ATCC_29799') = 'Firmicutes';
    taxMap('Roseburia_intestinalis_L1_82') = 'Firmicutes';
    taxMap('Roseburia_inulinivorans_DSM_16841') = 'Firmicutes';
    taxMap('Ruminococcus_bromii_L2_63') = 'Firmicutes';
    taxMap('Ruminococcus_callidus_ATCC_2776001') = 'Firmicutes';
    taxMap('Ruminococcus_torques_ATCC_27756') = 'Firmicutes';
    taxMap('Shigella_flexneri_2002017') = 'Proteobacteria';
    taxMap('Streptococcus_salivarius_JIM8777') = 'Firmicutes';
    taxMap('Subdoligranulum_variabile_DSM_15176') = 'Firmicutes';
    taxMap('Veillonella_atypica_ACS_049_V_Sch6') = 'Firmicutes';
end
