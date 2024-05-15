function [Y, Y_hat, model, snr, tgtSnr, vnoise]=prepProdSimu_makedataPP_region(nVox, tgtSnr, varargin)

%=======================================================================================
% [Y]=prepProdSimu_makedataPP_region(R, varargin)
%
% Modelling data for multivariate classification for prepProd
% Models regions based on the number of voxels and signal to noise ratio
% (snr ratio calculated by noiseNormalizeBeta_RY)
%
% Inputs:
%  nVox;    % number of voxels in target region (can be avg acros subj)
%  snr;     % signal to noise ratio of region (calculated by
%               noiseNormalizeBeta_RY)
%
% Output:
%  Y          % Resulting data including noise
%  model      % 'Perfect' model, excludes noise
%  residual   % data - model
%  Yhat       % Multivariately prewhitened Y (see Walther et al 2016)
%
% RY 10/2022

%%%Calculate required input noise for each region to achieve target snr.
%%%Signal is consistent across simulations. To adjust snr, we adjust the
%%%residuals (calculated as Y - X*Beta). Residuals scale linearly with the
%%%amount of noise input into the simulations (by a factor of ~120, checked
%%%by seeing squared residuals at different noise levels: mean(sum(res.^2))
simuSignal = 0.8556; %consistent signal across simulations - calculated as spm.xX.trRV*mean(diag(spm.xX.Bcov))
resScaling = 48;     %consistent scaling from residual to noise input (residual / vnoise)
%found in noiseNormalizeBeta_RY

resTarget = simuSignal / tgtSnr;
vnoise    = resTarget / resScaling; %input = tgtres / resScaling (see above)
nv        = nVox; %n voxels according to input

%Change the below parameters to suit desired simulation
classes     = 4;
nr          = 6;  % trial type/class repetitions
np          = 8;  %parameters or treatment effects (ncolumn of design matrix)

prepProd    = 1;         %0: prep = prod, 1: prep ~= prod
prodscaling = 1;         %scaling of prod relative to prep (multiplicative) - prod = prep * prodscaling
vord        = [0.3 0.1]; %distance between two orders
vtemp       = [0.4 0.7]; %distance between two timings
vinter      = [0.6 0.8]; %distance between four sequences

vararginoptions(varargin,{'prepProd','prodscaling','vord','vtemp','vinter'});

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

model    = [];

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
        
        R = X*U;                                    %generate prod from same multivar rand distribution, first without noise
        R = R + prodscaling;                        %then scale by input variable...
        R = R + normrnd(0,sqrt(vnoise),classes,nv); %then add the noise (so we don't also scale the noise, we add it after)
        model    = [model; X*U; (X*U)+prodscaling]; %add to model run-by-run
        
        
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
        R = R + prodscaling;
        PR = [P;R];
        Y = [Y ; PR]; %concatenate
        model    = [model; X*Uprep; ((X*Uprod) + prodscaling)]; %add to model run-by-run
    end%for each run
else
    error('Provide 1 or 0 to prepProd input')
end

% SPM.xY.VY  - nScan x 1 struct array of file handles
% SPM.xX     - structure containing design matrix information
% SPM.xX.W   - optional whitening/weighting matrix
% SPM.xVi    - structure describing intrinsic non-sphericity
% SPM.xM     - structure containing masking information

%%%Run prewhitening, based on rsa.spm.noiseNormalizeBeta.mat (rsa toolbox)
%%%first we have to build a version of some SPM variables:
spm.xX.X    = [Z;Z];    [nScan, ~] = size(spm.xX.X);

spm.xX.trRV = spm_SpUtil('trRV',spm.xX.X);
spm.xX.K    = 1;
spm.xX.W    = speye(nScan,nScan);
spm.xVi     = struct('form', 'i.i.d.',...
                     'V',    speye(nScan,nScan));
spm.xX.xKXs = spm_sp('Set',spm_filter(spm.xX.K,spm.xX.W*spm.xX.X));    % KWX
spm.xX.pKX  = spm_sp('x-',spm.xX.xKXs);
spm.xX.V    = spm_filter(spm.xX.K,spm_filter(spm.xX.K,spm.xX.W*spm.xVi.V*spm.xX.W')'); % KWVW'K'
spm.xX.Bcov = spm.xX.pKX*spm.xX.V*spm.xX.pKX';
%for info on spm struct see people.duke.edu/~njs28/spmdatastructure.htm

beta_hat = Y;
KWY      = spm_filter(spm.xX.K,spm.xX.W*Y);                               %%% filter out low-frequence trends in Y
%res      = spm_sp('r',spm.xX.xKXs,KWY); %%% residuals: res  = Y - X*beta
res      = Y - model; %%% residuals: res  = Y - X*beta

%%%Then run the prewhitening as in the RSA toolbox
% in the scaling of the noise, take into account mean beta-variance
Opt.shrinkage = [];
[Sw_reg,~,~]=rsa.stat.covdiag(res,spm.xX.trRV/mean(diag(spm.xX.Bcov)),'shrinkage',Opt.shrinkage);   %%% regularize Sw_hat through optimal shrinkage

% Postmultiply by the inverse square root of the estimated matrix
[V,L]=eig(Sw_reg);  % in the scaling of the noise, take into account mean beta-variance
l=diag(L);
sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
Y_hat=beta_hat*sq;

if (nargout>1)
    snr = spm.xX.trRV*mean(diag(spm.xX.Bcov)) / mean(sum(res.^2));
end

% figure
% imagesc(Y)