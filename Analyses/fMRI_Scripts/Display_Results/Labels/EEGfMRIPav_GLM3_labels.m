%% Names and labels of contrasts

GLMlabels = { % sID, [filter to select]
                1, 'Cue valence (Win - Avoid)', '\beta_{Win} - \beta_{Avoid} (a.u.)', 'pos'; 
                2, 'Executed action (Go - NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                3, 'Correct action (Go - NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos';
                4, 'Incorrect action (Go - NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                5, 'Accuracy (Correct - Incorrect)', '\beta_{Correct} - \beta_{Incorrect} (a.u.)', 'pos'; 
                6, 'Conflict (Incongruent - Congruent)', '\beta_{Incongruent} - \beta_{Congruent} (a.u.)', 'pos'; 
                7, 'Correct & Conflict (Incongruent - Congruent)', '\beta_{Correct Incongruent} - \beta_{Correct Congruent} (a.u.)', 'pos'; 
                8, 'InCorrect & Conflict (Incongruent - Congruent)', '\beta_{Incorrect Incongruent} - \beta_{Incorrect Congruent} (a.u.)', 'pos'; 
                9, 'Conflict & Accuracy (Correct - Incorrect)', '\beta_{Conflict Correct} - \beta_{Conflict Incorrect} (a.u.)', 'pos'; 
                10, 'EEG Alpha Synchronization (150 - 350 ms)', '\beta_{EEG Alpha Synchronization 150 -350 ms} (a.u.)', 'pos'; 
                11, 'EEG Theta Synchronization (500 - 1300 ms)', '\beta_{EEG Theta Synchronization 500 -1300 ms} (a.u.)', 'pos'; 
                };

rootdir  = '/project/3017042.02'; # root directory--needs to be adapted to users' folder structure
filename = fullfile(rootdir,'Analyses/fMRI_Scripts/Display_Results/Labels/GLM3labels.mat';
save(filename,'GLMlabels')
