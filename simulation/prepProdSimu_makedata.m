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
end


% Generate variance-covariance matrix G:
h = [vtemp vord vinter]; %Theta (parameter variance)
T = blockdiag(zeros(2),eye(2),zeros(4)).*h(1);
O = blockdiag(eye(2),zeros(2),zeros(4)).*h(2);
I = blockdiag(zeros(2),zeros(2),eye(4)).*h(3);


G = T+O+I;
%imagesc(G);

U = mvnrnd(zeros(np,1),G,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
Y = Z*U + normrnd(0,sqrt(vnoise),ns,nv);

