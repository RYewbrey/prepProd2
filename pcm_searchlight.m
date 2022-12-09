function pcm_searchlight(Y,M,run,c)

[Tgroup,theta] = pcm_fitModelGroup(Y,M,run,c,'runEffect',runEffect,'fitScale',1);

% Fit the models through cross-subject crossvalidation
[Tcross,thetaCr] = pcm_fitModelGroupCrossval(Y,M,run,c,'runEffect',runEffect,'groupFit',theta,'fitScale',1);