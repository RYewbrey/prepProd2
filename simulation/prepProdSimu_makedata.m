function [Y]=prepProdSimu_makedata(varargin)

%==========================================================================
% Modelling data for multivariate classification for prepProd
% Used by function prepProdSimu_runsimu
%
vtemp=2;     % variance across different timings
vord=1;      % variance across different orders
vinter=0.01; % variance across different sequences
vnoise=1;    % variance of the noise
classes=4;
nr   = 6;%8;   % trial type/class repetitions
ns   = classes*nr;  % trials/samples <<
nv   = 160;%160; % voxels
np   = 8;  % parameters or treatment effects (ncolumn of design matrix)



vararginoptions(varargin,{'nr','np','vtemp','vord','vinter','vnoise'});

%Generate design matrix: (9 x 8) trials x 15 effects
%  order   timing  sequence
Z =[1 0    1 0     1 0 0 0;
    1 0    0 1     0 1 0 0;
    0 1    1 0     0 0 1 0;
    0 1    0 1     0 0 0 1];

X =[1 0    1 0     1 0 0 0;
    1 0    0 1     0 1 0 0;
    0 1    1 0     0 0 1 0;
    0 1    0 1     0 0 0 1];




% repeat nr of times:
for i=1:(nr-1)
    Z = [Z;X];
end;


% Generate variance-covariance matrix G:
h = [vtemp vord vinter]; %Theta (parameter variance)
T = blockdiag(zeros(2),eye(2),zeros(4)).*h(1);
O = blockdiag(eye(2),zeros(2),zeros(4)).*h(2);
I = blockdiag(zeros(2),zeros(2),eye(4)).*h(3);


G = T+O+I;
%imagesc(G);

U = mvnrnd(zeros(np,1),G,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
Y = Z*U + normrnd(0,sqrt(vnoise),ns,nv);

%%% If the 9 patterns scale 
% pattern=Z*U;
% pattern=pattern(1,:);
% minVal=abs(min(pattern)); minVal=repmat(minVal,1,length(pattern));
% pattern=pattern+(1+minVal); %shift pattern to > 0
% pattern=repmat(pattern,54,1);
%  scaleFactor=randperm(9)'*0.05;%0.07;
% %   scaleFactor=(1:9)'*0.005;
% scaleFactor=repmat(scaleFactor,nr,nv);
% Y=pattern.*scaleFactor + normrnd(0,sqrt(vnoise),ns,nv);


% %plot 3d:
% cond=(1:9)';
% D.cond=repmat(cond,nr,1);
% D.Y1=Y(:,1); D.Y2=Y(:,2);
% scatterplot(D.Y1,D.Y2,...
% 'split',D.cond);

%%%%%Mean T1, T2, T3 averaged over O1, O2, O3 Additional classification analysis
k=0;
for i=1:3:ns
    k=k+1;
    YtempMean(k,:)=mean(Y(i:(i+2),:),1);
end;

%%Subtract mean T1, T2, T2, from Y, respectively
YtempMean2Subtract=kron(YtempMean,[1;1;1]); %triple mean of voxels
YcombiCorrected4Temp=Y-YtempMean2Subtract;
%figure;
