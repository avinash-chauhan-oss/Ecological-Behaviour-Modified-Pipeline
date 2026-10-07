% Plot2.m - Monoculture vs Community Growth Scatter

commFile = 'Community_Results_Wide_52Models.csv';
monoFile = 'minReactModels/minReactModels_WideSummary.csv';

T_comm = readtable(commFile, 'PreserveVariableNames', true);
T_mono = readtable(monoFile, 'PreserveVariableNames', true);

diets = {'High_Fiber', 'Western'};
monoPrefix = {'HF', 'WES'};
scenarios = {'WT', 'MRM_cut1', 'MRM_cut05', 'MRM_ext'};
monoScenarios = {'WT_Growth', 'MRM_cut1', 'MRM_cut05', 'MRM_ext'};

numPairs = height(T_comm);

figure('Position', [100, 100, 800, 1000]);
plotIdx = 1;

for s = 1:length(scenarios)
    scen = scenarios{s};
    mScen = monoScenarios{s};
    
    for d = 1:length(diets)
        diet = diets{d};
        pref = monoPrefix{d};
        
        subplot(4, 2, plotIdx);
        hold on;
        
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
                gc1 = T_comm{p, ['WT_M1_Max_' diet]};
                gc2 = T_comm{p, ['WT_M2_Max_' diet]};
            else
                gc1 = T_comm{p, [scen '_M1_Max_' diet]};
                gc2 = T_comm{p, [scen '_M2_Max_' diet]};
            end
            
            scatter(gm1, gc1, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
            scatter(gm2, gc2, 30, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
        end
        
        plot([0 0.5], [0 0.5], 'k--');
        xlabel('Monoculture Growth');
        ylabel('Community Growth');
        title(sprintf('%s - %s', scen, diet), 'Interpreter', 'none');
        
        plotIdx = plotIdx + 1;
    end
end
sgtitle('Monoculture vs Community Growth');
