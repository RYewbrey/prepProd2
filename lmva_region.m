function lmva_region(region, betas, mva_func,varargin)

%MVA Function contained to region provided as a .nii file
%%% INPUT:
% region: .nii file of region mask
% betas: beta filenames of interest
% mva_func: multivariate analysis function, should take X (data, P x N)
%           as the first argument 
% varargin: c,run,train,test
%
%RY July 2022

params = {}; %Extra parameters pased to mva_func

vararginoptions(varargin,{'params'});

% Load region .nii file and get voxel index to feed to MVA function
V = spm_vol(region);
X=spm_read_vols(V);

[i,j,k]=ind2sub(size(X),find(X~=0)); %find all of the voxels that belong to the specified mask
vox=[i j k];

LI = sub2ind(V)

feval(mva_func,full(X(T.LI{i},:)'),params{:});