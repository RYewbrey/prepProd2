function prepProdSimu_parametricCurve(prepProd,vord,vtemp,vinter,vnoise,n)

% prepProdSimu_parametricCurve()
%
% plots PCM performance and likelihood based on simulated data with various parameters
%
%
%

%%% Paths
% addpath(genpath('G:\projectsBackup\rhys\prepProd2\matlab')); %Adjust! genpath command also loads subfolders
% addpath(genpath('D:\projects\toolboxes\tools')); %joern's extensions for spm
% addpath(genpath('D:\projects\toolboxes\userfun')); %joern's util tools (open source)
% addpath(genpath('D:\projects\toolboxes\rsatoolbox_matlab')); %RSA toolbox
% addpath(genpath('D:\projects\toolboxes\pcm_toolbox')); %PCM toolbox

%adjust directory to where simulated data is stored
simuDir = 'G:\projectsBackup\rhys\prepProd2\data\imaging\simulations';
modelDir=[simuDir '\models'];

sn = 3; %number of simulated subjects

%% load representational models            %% and format ready for PCM
models = prepProdSimu_loadModels(modelDir); M = prepProdSimu_prepModels(models);

%% label conditions
run = repmat(1:6,8,1); %labelling runs 1:6 (imaging runs)
run = reshape(run,1,[])';
c = repmat([1;2;3;4;5;6;7;8], 6,1); %labelling conditions 1:4 (sequences)

%% Generate data and run PCM in loop

%%%generate data
for i=1:sn %for 'subject' n
    for j=1:n %for iterations
        
        %make the data based on input specifications
        Y=prepProdSimu_makedataPP('prepProd',prepProd,'vord',vord,'vtemp',vtemp,'vinter',vinter,'vnoise',vnoise);
        
        D.data(:,:,j) = Y;
        clear('Y')
    end%for iterations
    
    A.data{i} = mean(D.data,3);
end%for subject n

%%% PCM Functions

% Treat the run effect as random or fixed?
% We are using a fixed run effect here, as we are not interested in the
% activity relative the the baseline (rest) - so as in RSA, we simply
% subtract out the mean patttern across all conditions.
runEffect  = 'fixed';

% Fit the models on the group level
[Tgroup,theta] = pcm_fitModelGroup(Y,M,run,c,'runEffect',runEffect,'fitScale',1);

% Fit the models through cross-subject crossvalidation
[Tcross,thetaCr] = pcm_fitModelGroupCrossval(Y,M,run,c,'runEffect',runEffect,'groupFit',theta,'fitScale',1);


