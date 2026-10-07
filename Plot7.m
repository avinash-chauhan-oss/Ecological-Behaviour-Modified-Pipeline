% Plot7.m - Jaccard Similarity vs Interaction Type (Modified Pipeline)
% Shows interaction composition by Jaccard similarity tertiles for WT and MRMs across diets.
% Uses WT-based Jaccard similarity for categorization.

%% Configuration
commFile = 'Community_Results_Wide_52Models.csv';
monoFile = 'minReactModels/minReactModels_WideSummary.csv';
base_dir = pwd; % Load original .mat models from base directory

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

threshold = 0.1;

% 1: Positive (Mutualism, Commensalism), 2: Negative (Amensalism, Parasitism, Competition), 3: Neutral
colorMap = [0.2 0.8 0.2; 0.8 0.2 0.2; 0.8 0.8 0.8]; 

%% Calculate Jaccard Similarities
jaccardSim = zeros(height(commData), 1);

disp('Calculating Jaccard similarities...');
for i = 1:height(commData)
    names = strsplit(commData.Community{i}, '_AND_');
    sp1 = names{1};
    sp2 = names{2};
    
    S1 = load(fullfile(base_dir, [sp1 '.mat']));
    fn1 = fieldnames(S1); m1 = S1.(fn1{1}); if iscell(m1), m1 = m1{1}; end
    
    S2 = load(fullfile(base_dir, [sp2 '.mat']));
    fn2 = fieldnames(S2); m2 = S2.(fn2{1}); if iscell(m2), m2 = m2{1}; end
    
    rxns1 = strtrim(cellstr(m1.rxns));
    rxns2 = strtrim(cellstr(m2.rxns));
    
    intersectRxns = intersect(rxns1, rxns2);
    unionRxns = union(rxns1, rxns2);
    
    jaccardSim(i) = length(intersectRxns) / length(unionRxns);
end

% Tertiles
p33 = prctile(jaccardSim, 33.33);
p66 = prctile(jaccardSim, 66.67);
tertIdx = ones(height(commData), 1);
tertIdx(jaccardSim > p33 & jaccardSim <= p66) = 2;
tertIdx(jaccardSim > p66) = 3;

%% Process & Plot
for d = 1:length(dietLabels)
    diet = dietLabels{d};
    pref = dietPrefixes{d};
    
    compData = zeros(3, 3, 4);
    
    for i = 1:height(commData)
        names = strsplit(commData.Community{i}, '_AND_');
        sp1 = names{1};
        sp2 = names{2};
        t = tertIdx(i);
        
        wt_m1_mono = monoData{strcmp(monoData.Model, sp1), [pref '_WT_Growth']};
        wt_m2_mono = monoData{strcmp(monoData.Model, sp2), [pref '_WT_Growth']};
        wt_m1_comm = commData{i, ['WT_M1_Max_' diet]};
        wt_m2_comm = commData{i, ['WT_M2_Max_' diet]};
        
        alpha1WT = calcAlpha(wt_m1_comm, wt_m1_mono);
        alpha2WT = calcAlpha(wt_m2_comm, wt_m2_mono);
        typeWT = classifyMeta(alpha1WT, alpha2WT, threshold);
        compData(t, typeWT, 1) = compData(t, typeWT, 1) + 1;
        
        for c = 1:length(mrmCuts)
            cut = mrmCuts{c};
            mrm_m1_mono = monoData{strcmp(monoData.Model, sp1), [pref '_MRM_' cut]};
            mrm_m2_mono = monoData{strcmp(monoData.Model, sp2), [pref '_MRM_' cut]};
            mrm_m1_comm = commData{i, ['MRM_' cut '_M1_Max_' diet]};
            mrm_m2_comm = commData{i, ['MRM_' cut '_M2_Max_' diet]};
            
            alpha1MRM = calcAlpha(mrm_m1_comm, mrm_m1_mono);
            alpha2MRM = calcAlpha(mrm_m2_comm, mrm_m2_mono);
            typeMRM = classifyMeta(alpha1MRM, alpha2MRM, threshold);
            compData(t, typeMRM, c+1) = compData(t, typeMRM, c+1) + 1;
        end
    end
    
    figure('Name', sprintf('Jaccard vs Interaction - %s', dietTitles{d}), 'Position', [100, 100, 800, 600]);
    titles = [{'WT'}, mrmTitles];
    yLabels = {'Low', 'Medium', 'High'};
    
    maxTotal = max(sum(compData, 2), [], 'all');
    
    hold on;
    for sc = 1:4
        for t = 1:3
            counts = squeeze(compData(t, :, sc));
            total = sum(counts);
            if total > 0
                drawScaledPie(sc, 4-t, counts, colorMap, total, maxTotal, 0.4);
            end
        end
    end
    
    xlim([0.5 4.5]); ylim([0.5 3.5]);
    xticks(1:4); xticklabels(titles);
    yticks(1:3); yticklabels(flipud(yLabels')); ylabel('Jaccard Similarity Tertile');
    set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 10);
    title(sprintf('Interaction Composition by Similarity - %s', dietTitles{d}));
    
    legH = gobjects(3,1);
    for k=1:3
        legH(k) = fill(NaN, NaN, colorMap(k,:));
    end
    legend(legH, {'Positive', 'Negative', 'Neutral'}, 'Location', 'northeastoutside');
    
    axis normal; box on; hold off;
end

function metaType = classifyMeta(alpha1, alpha2, threshold)
    type = classifyInteraction(alpha1, alpha2, threshold);
    if ismember(type, [1, 2])
        metaType = 1;
    elseif ismember(type, [3, 4, 5])
        metaType = 2;
    else
        metaType = 3;
    end
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
            fill(xf, yf, colors(k,:), 'EdgeColor', 'k', 'LineWidth', 0.5);
            startAngle = endAngle;
        end
    end
end
