%% Names and labels of contrasts

GLMlabels = { % sID, [filter to select]
                1, 'Cue valence (Win - Avoid)', '\beta_{Win} - \beta_{Avoid} (a.u.)', 'pos'; 
                2, 'Executed action (Go - NoGo)', '\beta_{Go} - \beta_{NoGo} (a.u.)', 'pos'; 
                3, 'Conflict (Incongruent - Congruent)', '\beta_{Incongruent} - \beta_{Congruent} (a.u.)', 'pos'; 
                4, 'Hand (Left - Right)', '\beta_{Left Hand} - \beta_{Right Hand} (a.u.)', 'pos'; 
                5, 'Valence Go actions', '\beta_{Go2Win} - \beta_{Go2Avoid} (a.u.)', 'pos'; 
                };

rootdir  = '/project/3017042.02'; # root directory--needs to be adapted to users' folder structure
filename = fullfile(rootdir,'Analyses/fMRI_Scripts/Display_Results/Labels/GLM1labels.mat';
save(filename,'GLMlabels')
