% Plot5.m - Sankey Transition Diagram (Modified Pipeline)
% Shows interaction transitions from WT to MRM models across diets.

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

intCategories = {'Mutualism', 'Commensalism', 'Amensalism', 'Parasitism', 'Competition', 'Neutral'};
threshold = 0.1;

%% Plotting
for d = 1:length(dietLabels)
    diet = dietLabels{d};
    pref = dietPrefixes{d};
    
    for c = 1:length(mrmCuts)
        cut = mrmCuts{c};
        
        transitionMatrix = zeros(6, 6);
        
        for i = 1:height(commData)
            names = strsplit(commData.Community{i}, '_AND_');
            sp1 = names{1};
            sp2 = names{2};
            
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
            idxWT = classifyInteraction(alpha1WT, alpha2WT, threshold);
            
            alpha1MRM = calcAlpha(mrm_m1_comm, mrm_m1_mono);
            alpha2MRM = calcAlpha(mrm_m2_comm, mrm_m2_mono);
            idxMRM = classifyInteraction(alpha1MRM, alpha2MRM, threshold);
            
            transitionMatrix(idxWT, idxMRM) = transitionMatrix(idxWT, idxMRM) + 1;
        end
        
        figure('Name', sprintf('Transitions - %s - WT to %s', dietTitles{d}, mrmTitles{c}), 'Position', [100, 100, 800, 600]);
        drawSankey(transitionMatrix, intCategories, intCategories);
        title(sprintf('%s: WT \x2192 %s', dietTitles{d}, mrmTitles{c}), 'Interpreter', 'none');
    end
end

function drawSankey(matrix, labelsLeft, labelsRight)
    numLeft = length(labelsLeft);
    numRight = length(labelsRight);
    
    leftTotals = sum(matrix, 2);
    rightTotals = sum(matrix, 1);
    
    maxTotal = max(sum(leftTotals), sum(rightTotals));
    if maxTotal == 0, return; end
    
    nodeWidth = 0.1;
    ySpacing = maxTotal * 0.05;
    
    colors = lines(numLeft);
    
    yLeftStart = zeros(numLeft, 1);
    yRightStart = zeros(numRight, 1);
    
    currY = 0;
    for i = 1:numLeft
        yLeftStart(i) = currY;
        currY = currY + leftTotals(i) + ySpacing;
    end
    
    currY = 0;
    for i = 1:numRight
        yRightStart(i) = currY;
        currY = currY + rightTotals(i) + ySpacing;
    end
    
    hold on;
    
    for i = 1:numLeft
        if leftTotals(i) > 0
            rectangle('Position', [0, yLeftStart(i), nodeWidth, leftTotals(i)], 'FaceColor', colors(i,:), 'EdgeColor', 'k');
            text(-0.05, yLeftStart(i) + leftTotals(i)/2, sprintf('%s (%d)', labelsLeft{i}, leftTotals(i)), 'HorizontalAlignment', 'right');
        end
    end
    
    for i = 1:numRight
        if rightTotals(i) > 0
            rectangle('Position', [1-nodeWidth, yRightStart(i), nodeWidth, rightTotals(i)], 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
            text(1.05, yRightStart(i) + rightTotals(i)/2, sprintf('%s (%d)', labelsRight{i}, rightTotals(i)), 'HorizontalAlignment', 'left');
        end
    end
    
    yLeftFlow = yLeftStart;
    yRightFlow = yRightStart;
    
    for i = 1:numLeft
        for j = 1:numRight
            val = matrix(i, j);
            if val > 0
                xFlow = linspace(nodeWidth, 1-nodeWidth, 100);
                curve = 1 ./ (1 + exp(-10 * (xFlow - 0.5) / (1-2*nodeWidth)));
                curve = curve - min(curve);
                curve = curve / max(curve);
                
                yFlowTop = yLeftFlow(i) + val + (yRightFlow(j) + val - (yLeftFlow(i) + val)) * curve;
                yFlowBottom = yLeftFlow(i) + (yRightFlow(j) - yLeftFlow(i)) * curve;
                
                fill([xFlow, fliplr(xFlow)], [yFlowTop, fliplr(yFlowBottom)], colors(i,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                
                yLeftFlow(i) = yLeftFlow(i) + val;
                yRightFlow(j) = yRightFlow(j) + val;
            end
        end
    end
    
    xlim([-0.3, 1.3]);
    ylim([-ySpacing, max(currY, max(yLeftStart(end)+leftTotals(end), yRightStart(end)+rightTotals(end))) + ySpacing]);
    axis off;
    hold off;
end
