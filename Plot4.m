% Plot4.m - 52x52 Heatmap of Interactions (Modified Pipeline)
% Shows Wild-Type (WT) interactions in the upper triangle and MRM interactions in the lower triangle.
% Generated for each diet and each MRM cutoff (6 heatmaps total).

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

colorMap = [
    0 1 0;      % 1: Mutualism (Green)
    0.5 1 0.5;  % 2: Commensalism (Light Green)
    1 0.5 0.5;  % 3: Amensalism (Light Red)
    1 0.5 0;    % 4: Parasitism (Orange)
    1 0 0;      % 5: Competition (Red)
    0.8 0.8 0.8 % 6: Neutral (Gray)
];

speciesSet = {};
for i = 1:height(commData)
    names = strsplit(commData.Community{i}, '_AND_');
    speciesSet = [speciesSet, names{1}, names{2}];
end
speciesList = unique(speciesSet);
numSpecies = length(speciesList);
speciesMap = containers.Map(speciesList, 1:numSpecies);

threshold = 0.1;

%% Plotting Loop
for d = 1:length(dietLabels)
    diet = dietLabels{d};
    pref = dietPrefixes{d};
    
    for c = 1:length(mrmCuts)
        cut = mrmCuts{c};
        
        matrixWT = 6 * ones(numSpecies, numSpecies); % Default Neutral
        matrixMRM = 6 * ones(numSpecies, numSpecies);
        
        % Populate matrices
        for i = 1:height(commData)
            names = strsplit(commData.Community{i}, '_AND_');
            sp1 = names{1};
            sp2 = names{2};
            
            idx1 = speciesMap(sp1);
            idx2 = speciesMap(sp2);
            
            wt_m1_mono = monoData{strcmp(monoData.Model, sp1), [pref '_WT_Growth']};
            wt_m2_mono = monoData{strcmp(monoData.Model, sp2), [pref '_WT_Growth']};
            wt_m1_comm = commData{i, ['WT_M1_Max_' diet]};
            wt_m2_comm = commData{i, ['WT_M2_Max_' diet]};
            
            mrm_m1_mono = monoData{strcmp(monoData.Model, sp1), [pref '_MRM_' cut]};
            mrm_m2_mono = monoData{strcmp(monoData.Model, sp2), [pref '_MRM_' cut]};
            mrm_m1_comm = commData{i, ['MRM_' cut '_M1_Max_' diet]};
            mrm_m2_comm = commData{i, ['MRM_' cut '_M2_Max_' diet]};
            
            alpha1WT = calcAlpha(wt_m1_comm, wt_m1_mono);
            alpha2WT = calcAlpha(wt_m2_comm, wt_m2_mono);
            typeWT = classifyInteraction(alpha1WT, alpha2WT, threshold);
            matrixWT(idx1, idx2) = typeWT;
            
            alpha1MRM = calcAlpha(mrm_m1_comm, mrm_m1_mono);
            alpha2MRM = calcAlpha(mrm_m2_comm, mrm_m2_mono);
            typeMRM = classifyInteraction(alpha1MRM, alpha2MRM, threshold);
            matrixMRM(idx2, idx1) = typeMRM; 
        end
        
        combinedMatrix = triu(matrixWT, 1) + tril(matrixMRM, -1);
        for i = 1:numSpecies
            combinedMatrix(i,i) = NaN; 
        end
        
        figure('Name', sprintf('%s - WT vs %s', dietTitles{d}, mrmTitles{c}), 'Position', [100, 100, 800, 800]);
        imagesc(combinedMatrix);
        colormap([colorMap; 1 1 1]); 
        clim([1 7]);
        
        title(sprintf('%s: WT (Upper) vs %s (Lower)', dietTitles{d}, mrmTitles{c}), 'Interpreter', 'none');
        
        xticks(1:numSpecies);
        yticks(1:numSpecies);
        xticklabels(strrep(speciesList, '_', ' '));
        yticklabels(strrep(speciesList, '_', ' '));
        xtickangle(90);
        set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 6);
        
        hold on;
        legNames = {'Mutualism', 'Commensalism', 'Amensalism', 'Parasitism', 'Competition', 'Neutral'};
        for k = 1:6
            plot(NaN, NaN, 's', 'MarkerSize', 10, 'MarkerFaceColor', colorMap(k,:), 'MarkerEdgeColor', 'k');
        end
        legend(legNames, 'Location', 'northeastoutside');
        hold off;
        
        axis square;
    end
end
