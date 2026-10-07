% Plot1.m - Interaction Types Stacked Bar

commFile = 'Community_Results_Wide_52Models.csv';
monoFile = 'minReactModels/minReactModels_WideSummary.csv';

T_comm = readtable(commFile, 'PreserveVariableNames', true);
T_mono = readtable(monoFile, 'PreserveVariableNames', true);

diets = {'High_Fiber', 'Western'};
monoPrefix = {'HF', 'WES'};
scenarios = {'WT', 'MRM_cut1', 'MRM_cut05', 'MRM_ext'};
monoScenarios = {'WT_Growth', 'MRM_cut1', 'MRM_cut05', 'MRM_ext'};

numPairs = height(T_comm);
th = 0.1;

counts = zeros(length(diets) * length(scenarios), 6); 

rowIdx = 1;
for d = 1:length(diets)
    diet = diets{d};
    pref = monoPrefix{d};
    
    for s = 1:length(scenarios)
        scen = scenarios{s};
        mScen = monoScenarios{s};
        
        for p = 1:numPairs
            names = strsplit(T_comm.Community{p}, '_AND_');
            m1_name = names{1};
            m2_name = names{2};
            
            idx1 = find(strcmp(T_mono.Model, m1_name));
            idx2 = find(strcmp(T_mono.Model, m2_name));
            
            if isempty(idx1) || isempty(idx2)
                continue;
            end
            
            gm1 = T_mono{idx1, [pref '_' mScen]};
            gm2 = T_mono{idx2, [pref '_' mScen]};
            
            if strcmp(scen, 'WT')
                colM1 = ['WT_M1_Max_' diet];
                colM2 = ['WT_M2_Max_' diet];
            else
                colM1 = [scen '_M1_Max_' diet];
                colM2 = [scen '_M2_Max_' diet];
            end
            
            gc1 = T_comm{p, colM1};
            gc2 = T_comm{p, colM2};
            
            a1 = calcAlpha(gc1, gm1);
            a2 = calcAlpha(gc2, gm2);
            
            type = classifyInteraction(a1, a2, th);
            counts(rowIdx, type) = counts(rowIdx, type) + 1;
        end
        rowIdx = rowIdx + 1;
    end
end

figure('Position', [100, 100, 800, 600]);
bar(counts, 'stacked');
set(gca, 'XTick', 1:8);
labels = {'HF WT', 'HF cut1', 'HF cut05', 'HF ext', 'WES WT', 'WES cut1', 'WES cut05', 'WES ext'};
set(gca, 'XTickLabel', labels);
xtickangle(45);
ylabel('Number of Pairs');
title('Interaction Types Across Diets and Scenarios');
legend({'Mutualism (+/+)', 'Commensalism (+/0)', 'Neutralism (0/0)', 'Amensalism (-/0)', 'Parasitism (+/-)', 'Competition (-/-)'}, 'Location', 'eastoutside');
