function [Y]=prepProdSimu_makedataPP(varargin)

%=======================================================================================
% [Y]=prepProdSimu_makedataPP(varargin)
%
% Modelling data for multivariate classification for prepProd
% Used by function prepProdSimu_runsimu
%
% The following inputs can be provided:
% prepProd; % preparation and production phases correlate? 0 yes 1 no
%             if 1, provide 1x2 vectors for vord, vtemp, and vinter
% vord;      % variance across different orders
% vtemp;     % variance across different timings
% vinter; % variance across different sequences
% vnoise=;    % variance of the noise
%
% RY 10/2022

%Change the below parameters to suit desired simulation
classes=4;
nr   = 6;%8;   % trial type/class repetitions
ns   = classes*nr;  % trials/samples <<
nv   = 1000; %160 voxels
% nv   = 4000;
np   = 8;  %parameters or treatment effects (ncolumn of design matrix)

vararginoptions(varargin,{'prepProd','vord','vtemp','vinter','vnoise'});

%Generate design matrix: (9 x 8) trials x 15 effects
%  order   timing  sequence    prep/prod
Z =[1 0    1 0     1 0 0 0;
    1 0    0 1     0 1 0 0;
    0 1    1 0     0 0 1 0;
    0 1    0 1     0 0 0 1];

X = Z;

% repeat nr of times:
for i=1:(nr-1)
    Z = [Z;X];
end

if ~prepProd %if prep and prod are chosen to come from the same distribution
    % Generate variance-covariance matrix G:
    h = [vord(1,1) vtemp(1,1) vinter(1,1)]; %Theta (parameter variance)
    O = blockdiag(eye(2),zeros(2),zeros(4)).*h(1);
    T = blockdiag(zeros(2),eye(2),zeros(4)).*h(2);
    I = blockdiag(zeros(2),zeros(2),eye(4)).*h(3);
    
    G = O+T+I;
    %imagesc(G);
    Y = [];
    
    for i=1:nr%for each run
        if i==1
            U = mvnrnd(zeros(np,1),G,nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
        end
        P = X*U + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
        R = X*U + normrnd(0,sqrt(vnoise),classes,nv); %generate prod from same multivar rand distribution
        PR = [P;R];
        Y = [Y ; PR]; %concatenate
    end%for each run
    
elseif prepProd == 1 %if prep and prod are to come from different distributions
    for i=1:2%for prep & prod
        
        % Generate variance-covariance matrix G:
        h = [vtemp(i) vord(i) vinter(i)]; %Theta (parameter variance)
        T = blockdiag(zeros(2),eye(2),zeros(4)).*h(1);
        O = blockdiag(eye(2),zeros(2),zeros(4)).*h(2);
        I = blockdiag(zeros(2),zeros(2),eye(4)).*h(3);
        
        G(:,:,i) = T+O+I;
    end%for prep & prod
    
    Y = [];
    
    for i=1:nr%for each run
        if i==1
            Uprep = mvnrnd(zeros(np,1),G(:,:,1),nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
            Uprod = mvnrnd(zeros(np,1),G(:,:,2),nv)'; % -> Diedrichsen et al. 2011: distribute variance across voxels
        end
        P = X*Uprep + normrnd(0,sqrt(vnoise),classes,nv); %generate prep
        R = X*Uprod + normrnd(0,sqrt(vnoise),classes,nv); %generate prod from different multivar rand distribution
        PR = [P;R];
        Y = [Y ; PR]; %concatenate
    end%for each run
    
else
    error('Provide 1 or 0 to prepProd input')
end

% figure
% imagesc(Y)