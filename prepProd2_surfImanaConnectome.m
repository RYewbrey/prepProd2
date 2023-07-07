function prepProd2_surfImanaConnectome(what,varargin)

% addpath(genpath('E:\projects\rhys\prepProd2\matlab')); %Ajust! loaded with subdirectories (genpath command)
% addpath(genpath('E:\projects\toolboxes\spm12'));
% addpath(genpath('E:\projects\toolboxes\tools')); %joern's extensions for spm
% addpath(genpath('E:\projects\toolboxes\userfun')); %joern's util tools (open source)
% addpath(genpath('E:\projects\toolboxes\region')); %joern's region toolbox for spm
% addpath(genpath('E:\projects\toolboxes\surfAnalysis')); %joern's toolbox for connectome workbench

%%%definition and variables

%% Data paths:
baseDir= 'E:\projects\rhys\prepProd2\data';
anatDir=[baseDir filesep 'imaging' filesep 'anatomicals'];
epiDir=[baseDir filesep 'imaging' filesep 'epi'];
glmDir=[baseDir filesep 'imaging' filesep 'GLM_firstlevel'];
behDir=[baseDir filesep 'behavioural'];
groupDir=[baseDir filesep 'imaging' filesep 'GLM_secondlevel'];
scndDir=[baseDir filesep 'imaging' filesep 'GLM_secondlevel' filesep 'data'];
regDir=[baseDir filesep 'imaging' filesep 'RegionOfInterest'];
freesurferDir=[baseDir filesep 'imaging' filesep 'surfaceFreesurfer'];
connectomeDir=[baseDir filesep 'imaging' filesep 'surfaceConnectome'];


subj = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 ...
    37 38 39 40 41 42];
subj_name={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16',...
    's17','s18','19','s20','s21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34',...
    's35','s36','s37','s38','s39','s40','s41','s42'}; %isotropic new data
anaSubj=[3 5 6 7 9 10 13 16 17 18 20 21 22 25 26 31 32 34 36 38 39 40 41 42]; %participants who met behavioural performance threshold



switch(what)
    
    
    case 'surf_reslice'
        
        %Reslices participant data into fs_LR_164 space (https://wiki.humanconnectome.org/download/attachments/63078513/Resampling-FreeSurfer-HCP_5_8.pdf)
        %Inputs surfaceFreesurfer folder for individual participant (after freesurfer recon) and outputs surfaceConnectome folder for individual.
        
        sn=varargin{1};
        
        surf_resliceFS2WB(subj_name{sn}, freesurferDir, [connectomeDir '\' subj_name{sn}])
    
    
    case 'surf_define_search'
        
        %Defines a surface searchlight on participant data which has been aligned to FS atlas using 'surf_reslice'
        %For use in later MVA_do functions
        
        
    
    
    
    
    
    
    
end
