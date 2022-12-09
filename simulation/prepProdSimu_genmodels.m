function prepProdSimu_genmodels(model)

%==============================================================================
% Simulates data across 1000 iterations to generate PCM models for prepProd
%
% prepProdSimu_genmodels('model')
% 'model' can be any of the following:
%     ordPrep   - order during prep only
%     ordProd   - order during prod only
%     ordMaint  - order stays online with same distribution
%     ordShift  - order stays online but distribution changes (a la Kaufman)
%     tempPrep  - timing during prep
%     tempProd  - timing during prod
%     tempMaint - timing stays online with same distribution
%     tempShift - timing stays online but distribution changes
%     intPrep   - integration during prep
%     intProd   - integration during prod
%     intMaint  - integration in both with same distribution
%     intShift  - integration in both but distribution changes
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

saveDir = 'G:\projectsBackup\rhys\prepProd2\data\imaging\simulations\models';

D=[];

%Change the below parameters to suit desired simulation
classes=4;
nr   = 6;%8;   % trial type/class repetitions
ns   = classes*nr;  % trials/samples <<
nv   = 160;%160; % voxels
np   = 8;  % parameters or treatment effects (ncolumn of design matrix)
nseq = 4; %number of sequences collapsed across phase
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
    
    case 'ordPrep'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [0 1 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [0 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - no control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'ordProd'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [0 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - no control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [0 1 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'ordMaint'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [0 1 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [0 1 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = prepU;
            %             prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'ordSwitch'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [0 1 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [0 1 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'tempPrep'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [1 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [0 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - no control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'tempProd'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [0 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - no control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [1 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'tempMaint'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [1 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [1 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = prepU;
            %             prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'tempSwitch'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [1 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [1 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'intPrep'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [0 0 1]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [0 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - no control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'intProd'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [0 0 0]; %Theta (parameter variance)
            %   [temp ord int]^ - no control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [0 0 1]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'intMaint'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [0 0 1]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [0 0 1]; %Theta (parameter variance)
            %   [temp ord int]^ - order control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = prepU;
            %             prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations
        
    case 'intSwitch'
        
        for i=1:n %for iterations
            % Generate variance-covariance matrix G:
            prepH = [0 0 1]; %Theta (parameter variance)
            %   [temp ord int]^ - int control
            prepT = blockdiag(zeros(2),eye(2),zeros(4)).*prepH(1);
            prepO = blockdiag(eye(2),zeros(2),zeros(4)).*prepH(2);
            prepI = blockdiag(zeros(2),zeros(2),eye(4)).*prepH(3);
            
            prepG = prepT+prepO+prepI;
            
            % Generate variance-covariance matrix G:
            prodH = [0 0 1]; %Theta (parameter variance)
            %   [temp ord int]^ - int control
            prodT = blockdiag(zeros(2),eye(2),zeros(4)).*prodH(1);
            prodO = blockdiag(eye(2),zeros(2),zeros(4)).*prodH(2);
            prodI = blockdiag(zeros(2),zeros(2),eye(4)).*prodH(3);
            
            prodG = prodT+prodO+prodI;
            
            prepU = mvnrnd(zeros(np,1),prepG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            prodU = mvnrnd(zeros(np,1),prodG,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            
            prepP = X*prepU + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
            prodP = X*prodU + normrnd(0,sqrt(vnoise),classes,nv); %generate prod
            PP = [prepP;prodP];
            Y = repmat(PP,nr,1);
            
            D.data{i,1} = Y;
            
        end%for iterations

end%switch

D.modelName = model;
%average data across iterations for RSA
D.data = cat(3, D.data{:});
D.meanData = mean(D.data,3);

dataDistance(:,:,1) = squareform(pdist(squeeze(D.meanData(1:8,:))));
dataDistance(:,:,2) = squareform(pdist(squeeze(D.meanData(9:16,:))));
dataDistance(:,:,3) = squareform(pdist(squeeze(D.meanData(17:24,:))));
dataDistance(:,:,4) = squareform(pdist(squeeze(D.meanData(25:32,:))));
dataDistance(:,:,5) = squareform(pdist(squeeze(D.meanData(33:40,:))));
dataDistance(:,:,6) = squareform(pdist(squeeze(D.meanData(41:48,:))));

meanDataDistance = mean(dataDistance,3);
D.distance = meanDataDistance;

figure
imagesc(meanDataDistance)
colorbar

if ~isdir(saveDir)
    mkdir(saveDir)
end
save([saveDir '\' model '_model_prepProdSimu.mat'], 'D')


