function varargout = prepProdSimu_runPCMPP(prepProd,vord,vtemp,vinter,vnoise,n)

% prepProdSimu_runPCM(prepProd,vtemp,vord,vinter,vnoise,n)
%
% compares models generated from simulated (no noise) data to new simulated data
%
% inputs (parameters of simulated data to test models against):
%   prepProd - prep correlates with prod, set to zero
%              prep doesn't correlate with prod, set to one
%
%   vord     - variance across order conditions
%
%   vtemp    - variance across timing conditions
%
%   vinter   - variance across four sequences
%
%   vnoise   - amount of noise in data set
%
%   n        - number of iterations for simulation
%
% models defined by prepProdSimu_genmodels, make sure you generate these first.
% specification for models is found below
% RY 10/2022

% addpath(genpath('G:\projectsBackup\rhys\prepProd2\matlab')); %Adjust! loaded with subdirectories (genpath command)
% addpath(genpath('D:\projects\toolboxes\tools')); %joern's extensions for spm
% addpath(genpath('D:\projects\toolboxes\userfun')); %joern's util tools (open source)
% addpath(genpath('D:\projects\toolboxes\rsatoolbox_matlab')); %RSA toolbox
% addpath(genpath('D:\projects\toolboxes\pcm_toolbox')); %PCM toolbox

%adjust directory to where simulated data is stored
simuDir = 'Z:\rhys\prepProd2\data\imaging\simulations';
modelDir=[simuDir '\models'];

sn = 5; %number of simulated subjects

D=[];
for i=1:sn %for 'subject' n
    for j=1:n %for iterations
        
        %make the data based on input specifications
        Y=prepProdSimu_makedataPP('prepProd',prepProd,'vord',vord,'vtemp',vtemp,'vinter',vinter,'vnoise',vnoise);
        
        D.data(:,:,j) = Y;
        clear('Y')
    end%for iterations
    
    A.data{i} = mean(D.data,3);
end%for subject n

run = repmat(1:6,8,1); %labelling runs 1:6 (imaging runs)
run = reshape(run,1,[])';

c = repmat([1;2;3;4;5;6;7;8], 6,1); %labelling conditions 1:4 (sequences)

models = prepProdSimu_loadModels(modelDir); [M, R] = prepProdSimu_prepModels(models,run,c);

for s=1:length(A.data) %estimate second moment matrix for each 'subject'
    G_hat(:,:,s)=pcm_estGCrossval(A.data{s},run,c);
end;
Gm = mean(G_hat,3); % Mean estimate


figure %plot empirical data
subplot(3,7,[1, 8]);
H = eye(8)-ones(8)/8;
imagesc(H*Gm*H');
title('Empirical Data')

C= pcm_indicatorMatrix('allpairs',(1:8)'); %multi-dimensional scaling
[COORD,l]=pcm_classicalMDS(Gm,'contrast',C);
subplot(3,7,15);
plot(COORD(:,1),COORD(:,2),'o');
axis equal;

%%% visualise models
%order
subplot(3,7,2);
imagesc(R.orderPrepModel.G);
title('Order Prep')

subplot(3,7,9);
imagesc(R.orderProdModel.G);
title('Order Prod')

subplot(3,7,16);
imagesc(R.orderMaintModel.G);
title('Order Maint')

%timing
subplot(3,7,3);
imagesc(R.timingPrepModel.G);
title('Timing Prep')

subplot(3,7,10);
imagesc(R.timingProdModel.G);
title('Timing Prod')

subplot(3,7,17);
imagesc(R.timingMaintModel.G);
title('Timing Maint')

%integrated
subplot(3,7,4);
imagesc(R.integratedPrepModel.G);
title('Integrated Prep')

subplot(3,7,11);
imagesc(R.integratedProdModel.G);
title('Integrated Prod')

subplot(3,7,18);
imagesc(R.integratedMaintModel.G);
title('Integrated Maint')

% Treat the run effect as random or fixed?
% We are using a fixed run effect here, as we are not interested in the
% activity relative the the baseline (rest) - so as in RSA, we simply
% subtract out the mean patttern across all conditions.
runEffect  = 'fixed';

Y = A.data; % function requires cell input, struct -> cell

% Fit the models on the group level
[Tgroup,theta] = pcm_fitModelGroup(Y,M,run,c,'runEffect',runEffect,'fitScale',1);

% Fit the models through cross-subject crossvalidation
[Tcross,thetaCr] = pcm_fitModelGroupCrossval(Y,M,run,c,'runEffect',runEffect,'groupFit',theta,'fitScale',1);

% Provide a plot of the crossvalidated likelihoods
subplot(3,7,[5 6 7 12 13 14 19 20 21]);
T = pcm_plotModelLikelihood_RY(Tcross,M,'upperceil',Tgroup.likelihood(:,length(M)));

title(['prepProd: ' num2str(prepProd) ', Ord: ' num2str(vord) ', Temp: ' num2str(vtemp) ', Int: ' num2str(vinter) ', Noise: ' num2str(vnoise)])

varargout={T,M};
