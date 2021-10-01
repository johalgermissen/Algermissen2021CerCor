%% Names and labels of contrasts

GLMlabels = { % sID, [filter to select]
                1, 'Cue valence all trials (Win - Avoid)', '\beta_{Win} - \beta_{Avoid} (a.u.)', 'pos'; 
                2, 'Executed action all trials (Go - NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                3, 'Conflict all trials (Incongruent - Congruent)', '\beta_{Incongruent} - \beta_{Congruent} (a.u.)', 'pos'; 
                4, 'Cue valence correct trials (Win - Avoid)', '\beta_{Win} - \beta_{Avoid} (a.u.)', 'pos'; 
                5, 'Executed action correct trials (Go - NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                6, 'Conflict correct trials (Incongruent - Congruent)', '\beta_{Incongruent} - \beta_{Congruent} (a.u.)', 'pos'; 
                7, 'Cue valence incorrect trials (Win - Avoid)', '\beta_{Win} - \beta_{Avoid} (a.u.)', 'pos'; 
                8, 'Executed action incorrect trials (Go - NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                9, 'Conflict incorrect trials (Incongruent - Congruent)', '\beta_{Incongruent} - \beta_{Congruent} (a.u.)', 'pos'; 
                10, 'Cue valence correct minus incorrect (Win - Avoid)', '\beta_{Win} - \beta_{Avoid} (a.u.)', 'pos'; 
                11, 'Executed action correct minus incorrect (Go - NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                12, 'Conflict correct minus incorrect (Incongruent - Congruent)', '\beta_{Incongruent} - \beta_{Congruent} (a.u.)', 'pos'; 
                13, 'Accuracy (Correct - incorrect)', '\beta_{Correct} - \beta_{Incorrect} (a.u.)', 'pos'; 
                14, 'Hand (Left - Right)', '\beta_{Left Hand} - \beta_{Right Hand} (a.u.)', 'pos'; 
                15, 'Outcome Valence', '\beta_{Positive} - \beta_{Negative} (a.u.)', 'pos'; 
                };

rootdir  = '/project/3017042.02'; # root directory--needs to be adapted to users' folder structure
filename = fullfile(rootdir,'Analyses/fMRI_Scripts/Display_Results/Labels/GLM2labels.mat';
save(filename,'GLMlabels')
