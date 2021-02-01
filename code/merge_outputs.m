% VERTEX-WISE STRUCTURAL COVARIANCE: merge outputs
%   This script combines outputs (network/smoothing/study) from main
%   script into a behavioural dataset for further analysis. Assumes behavioural
%   data has same number of rows/same order as outputs, behavioural data is
%   located in present working directory, and outputs are in a containing
%   folder called "netfeats".

% Options
sm = {'10';'20';'40'};
nets = {'VIS';'SOM'; 'DAN'; 'VAN'; 'LIM'; 'FPN'; 'DMN'};
study = 'Insight';
behfile = 'InsightBehData.mat';

% Code
load(behfile)

for i = 1:size(sm,1)
    for j = 1:size(nets,1)
        load(fullfile(pwd, 'netfeats', ['vertexConnectivity_' study '_' sm{i} '_' nets{j} '_features.mat']));
        beh.([nets{j} sm{i} '_modularity']) = gr.Q;
        beh.([nets{j} sm{i} '_meanPC']) = sum(gr.PC)';
        beh.([nets{j} sm{i} '_meanStrength']) = sum(gr.S)';
        beh.([nets{j} sm{i} '_eigenspect']) = gr.ES;
    end
end

writetable(beh,fullfile(pwd, 'vertexConnectivity_results.csv'))
save('vertexConnectivity_results.mat')
