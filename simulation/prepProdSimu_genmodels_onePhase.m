function prepProdSimu_genmodels_onePhase(model)

%==============================================================================
% Simulates data across 1000 iterations to generate PCM models for prepProd
%
% prepProdSimu_genmodels('model')
% 'model' can be any of the following:
%     ord   - order control
%     temp  - timing control
%     int   - integrated control
%
% RY 10/2022

% addpath(genpath('G:\projectsBackup\rhys\prepProd2\matlab')); %Adjust! loaded with subdirectories (genpath command)
% addpath(genpath('D:\projects\toolboxes\tools')); %joern's extensions for spm
% addpath(genpath('D:\projects\toolboxes\userfun')); %joern's util tools (open source)
% addpath(genpath('D:\projects\toolboxes\rsatoolbox_matlab')); %RSA toolbox
% addpath(genpath('D:\projects\toolboxes\pcm_toolbox')); %PCM toolbox

%generate all models:
% prepProdSimu_genmodels('ordPrep'); prepProdSimu_genmodels('ordProd'); prepProdSimu_genmodels('ordMaint'); prepProdSimu_genmodels('ordSwitch')
% prepProdSimu_genmodels('tempPrep'); prepProdSimu_genmodels('tempProd'); prepProdSimu_genmodels('tempMaint'); prepProdSimu_genmodels('tempSwitch')
% prepProdSimu_genmodels('intPrep'); prepProdSimu_genmodels('intProd'); prepProdSimu_genmodels('intMaint'); prepProdSimu_genmodels('intSwitch')

saveDir = 'G:\projectsBackup\rhys\prepProd2\data\imaging\simulations\models\onePhase';

D=[];

classes=4;
nr   = 6;%8;   % trial type/class repetitions
ns   = classes*nr;  % trials/samples <<
nv   = 160;%160; % voxels
np   = 8;  % parameters or treatment effects (ncolumn of design matrix)
vnoise = 0; %no noise for models
n    = 500; %n iterations

%Generate design matrix: (9 x 8) trials x 15 effects
%  order   timing  sequence    prep/prod
Z =[1 0    1 0     1 0 0 0;
    1 0    0 1     0 1 0 0;
    0 1    1 0     0 0 1 0;
    0 1    0 1     0 0 0 1];

X = Z; %matrix for just one run
Z = repmat(Z,[nr,1]); %matrix for number of runs
Y = [];

switch(model)
    
    case 'ord'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            H = [0 1 0]; %Theta (parameter variance)
            %[temp ord int]^ - order control
            T = blockdiag(zeros(2),eye(2),zeros(4)).*H(1);
            O = blockdiag(eye(2),zeros(2),zeros(4)).*H(2);
            I = blockdiag(zeros(2),zeros(2),eye(4)).*H(3);
            
            G = T+O+I;

            U = mvnrnd(zeros(np,1),G,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            Y = Z*U + normrnd(0,sqrt(vnoise),ns,nv);
            
            D.data{i,1} = Y;
            
        end%for iterations

    case 'temp'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            H = [1 0 0]; %Theta (parameter variance)
            %[temp ord int]^ - order control
            T = blockdiag(zeros(2),eye(2),zeros(4)).*H(1);
            O = blockdiag(eye(2),zeros(2),zeros(4)).*H(2);
            I = blockdiag(zeros(2),zeros(2),eye(4)).*H(3);
            
            G = T+O+I;

            U = mvnrnd(zeros(np,1),G,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            Y = Z*U + normrnd(0,sqrt(vnoise),ns,nv);
            
            D.data{i,1} = Y;
            
        end%for iterations

    case 'int'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            H = [0 0 1]; %Theta (parameter variance)
            %[temp ord int]^ - order control
            T = blockdiag(zeros(2),eye(2),zeros(4)).*H(1);
            O = blockdiag(eye(2),zeros(2),zeros(4)).*H(2);
            I = blockdiag(zeros(2),zeros(2),eye(4)).*H(3);
            
            G = T+O+I;

            U = mvnrnd(zeros(np,1),G,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            Y = Z*U + normrnd(0,sqrt(vnoise),ns,nv);
            
            D.data{i,1} = Y;
            
        end%for iterations

end%switch

D.modelName = model;
%average data across iterations for RSA
D.data = cat(3, D.data{:});
D.meanData = mean(D.data,3);

dataDistance(:,:,1) = squareform(pdist(squeeze(D.meanData(1:4,:))));
dataDistance(:,:,2) = squareform(pdist(squeeze(D.meanData(5:8,:))));
dataDistance(:,:,3) = squareform(pdist(squeeze(D.meanData(9:12,:))));
dataDistance(:,:,4) = squareform(pdist(squeeze(D.meanData(13:16,:))));
dataDistance(:,:,5) = squareform(pdist(squeeze(D.meanData(17:20,:))));
dataDistance(:,:,6) = squareform(pdist(squeeze(D.meanData(21:24,:))));

meanDataDistance = mean(dataDistance,3);
D.distance = meanDataDistance;

figure
imagesc(meanDataDistance)
colorbar

if ~isdir(saveDir)
    mkdir(saveDir)
end
save([saveDir '\' model '_model_prepProdSimu.mat'], 'D')