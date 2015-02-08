%% Jean: this script will test the PAA package by sequentially running 
% the synthetic data, the mouse and the cancer, both with PAA and PAA_light.

%% Synthetic data
clear all
close all
pause(2);
disp('*** Synthetic data');
global ForceNArchetypes; ForceNArchetypes = 3;
global lowIterations; lowIterations = 1;
exampleSynthetic

%% also try in PAAM_lite
disp('Done with PAA. Now trying PAA_lite with all algorithms and too many dimensions.');
for algNum = 1:5
    fprintf('Algorithm %d\n', algNum);
    for nA = 1:4
        ForceNArchetypes = nA;
        fprintf('... with %d archetypes\n', nA);
        close all
        pause(2);
        try
            [arc, arcFinal] = PAAM_lite(data, algNum, 5, discrAttrNames, ...
                discrAttr, 0, contAttrNames, contAttr, 0.2, 'Synthetic_enrichmentAnalysis');
        catch
            disp('Caught error message!');
        end
    end
end

%% Mouse data
clear all
close all
pause(2);
disp('*** Mouse data');
global ForceNArchetypes; ForceNArchetypes = 4;
global lowIterations; lowIterations = 1;
exampleMouse

%% also try in PAAM_lite
close all
pause(2);
[arc, arcFinal] = PAAM_lite(geneExpression, 2, 8, discrAttrNames, ...
    discrAttr, 0, contAttrNames, contAttr, 0.2, 'Mouse_enrichmentAnalysis');

%% Cancer data
clear all
close all
pause(2);
disp('*** Cancer data');
global ForceNArchetypes; ForceNArchetypes = 4;
global lowIterations; lowIterations = 1;
exampleCancer

%% also try in PAAM_lite
close all
pause(2);
[arc, arcFinal] = PAAM_lite(geneExpression, 1, 8, discrAttrNames, ...
    discrAttr, 0, contAttrNames, contAttr, 0.05, 'Cancer_enrichmentAnalysis');

%% Cancer RNAseq data
clear all
close all
pause(2);
disp('*** CancerRNAseq data');
global ForceNArchetypes; ForceNArchetypes = 4;
global lowIterations; lowIterations = 1;
exampleCancerRNAseq

%% also try in PAAM_lite
close all
pause(2);
[arc, arcFinal] = PAAM_lite(geneExpression, 1, 8, discrAttrNames, ...
    discrAttr, 0, contAttrNames, contAttr, 0.05, 'CancerRNAseq_enrichmentAnalysis');