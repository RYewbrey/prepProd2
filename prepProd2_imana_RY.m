function prepProd2_imana_RY(what,varargin)


%% Load relevant Matlab directories
% addpath(genpath('E:\projects\rhys\prepProd2\matlab')); %Adjust! loaded with subdirectories (genpath command)
% addpath(genpath('G:\projectsBackup\rhys\prepProd2\matlab')); %Adjust! loaded with subdirectories (genpath command)
% addpath('D:\projects\toolboxes\spm12');
% addpath(genpath('D:\projects\toolboxes\tools')); %joern's extensions for spm
% addpath(genpath('D:\projects\toolboxes\userfun')); %joern's util tools (open source)
% addpath(genpath('D:\projects\toolboxes\region-master')); %joern's region toolbox for spm
% addpath(genpath('D:\projects\toolboxes\spm12\marsbar-0.44')) %marsbar ROI toolbox
% addpath(genpath('D:\projects\toolboxes\spm12\toolbox\suit')); %SUIT toolbox for cerebellar analysis
% addpath(genpath('D:\projects\toolboxes\spm12\toolbox\DARTEL')); %DARTEL toolbox for SUIT reslice function
% addpath(genpath('D:\projects\toolboxes\permutest')); %permutest for crossSection analysis
% addpath(genpath('D:\projects\toolboxes\rsatoolbox_matlab')); %RSA toolbox
% addpath(genpath('D:\projects\toolboxes\pcm_toolbox')); %PCM toolbox

%%%definition and variables

%% Data paths:
% baseDir= 'E:\projects\rhys\prepProd2\data';
% baseDirR= 'E:\projects\rhys\prepProd2\data';
baseDir= 'G:\projectsBackup\rhys\prepProd2\data';
baseDirR= 'G:\projectsBackup\rhys\prepProd2\data';
rawDir=[baseDirR filesep 'imaging' filesep 'raw']; %original files before conversion
anatDir=[baseDir filesep 'imaging' filesep 'anatomicals'];
epiDir=[baseDir filesep 'imaging' filesep 'epi'];
glmDir=[baseDir filesep 'imaging' filesep 'GLM_firstlevel'];
behDir=[baseDir filesep 'behavioural'];
groupDir=[baseDir filesep 'imaging' filesep 'GLM_secondlevel'];
scndDir=[baseDir filesep 'imaging' filesep 'GLM_secondlevel' filesep 'data'];
roiDir=[baseDir filesep 'imaging' filesep 'ROI'];
regDir=[baseDir filesep 'imaging' filesep 'RegionOfInterest'];
suitDir=[baseDir filesep 'imaging' filesep 'suit'];
suitGroupDir=[baseDir filesep 'imaging' filesep 'suit_secondlevel'];
suitScndDir=[baseDir filesep 'imaging' filesep 'suit_secondlevel' filesep 'data'];
subcorticalSearchDir=[baseDir filesep 'imaging' filesep 'subcortical' filesep 'search'];
subcorticalSearchGroupDir=[baseDir filesep 'imaging' filesep 'subcortical_secondlevel' filesep 'search'];
subcorticalAreaDir=[baseDir filesep 'imaging' filesep 'subcortical' filesep 'area'];
crossSectionDir=[baseDir filesep 'imaging' filesep 'crossSection'];
pcmDir=[baseDir filesep 'imaging' filesep 'pcm'];
pcmGroupDir=[baseDir filesep 'imaging' filesep 'GLM_secondlevel' filesep 'PCM'];
modelDir=[baseDir filesep 'imaging' filesep 'simulations\models'];

%% Names and IDs
subj=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42];
anaSubj = [3,5,6,7,9,10,13,16,17,18,20,21,22,25,26,31,32,34,36,38,39,40,41,42];

subj_name={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24','s25','s26',...
    's27','s28','s29','s30','s31','s32','s33','s34','s35','s36','s37','s38','s39','s40','s41','s42'};

subjID={'placeholder','SS0094_09','SS0094_10','SS0094_11','SS0094_14','SS0094_13','SS0094_12','SS0094_15','SS0094_16','SS0094_17','SS0094_18','SS0094_19','SS0094_20',...
    'placeholder','SS0094_21','SS0094_22','SS0094_23','SS0094_24','placeholder','SS0094_25_2','SS0094_26','SS0094_27','placeholder','placeholder','SS0094_28','SS0094_29',...
    'placeholder','placeholder','placeholder','placeholder','SS0094_30','SS0094_31','SS0094_32_2','SS0094_33','placeholder','SS0094_34','SS0094_35','SS0094_36','SS0094_37',...
    'SS0094_38','SS0094_39','SS0094_40'}; %folder

subjIDsession={'placeholder','SS0094_09_20201113_631887115','SS0094_10_20201119_632400213','SS0094_11_20201119_632410210','SS0094_14_20201126_633022974',...
    'SS0094_13_20201126_633008222','SS0094_12_20201126_632998125','SS0094_15_20201203_631017618','PrepProd2_s09','PrepProd2_s10','PrepProd2_s11','SS0094_19_s12',...
    'SS0094_20_s13','placeholder','SS0094_21_s15','SS0094_22_s16','SS0094_23','SS0094_24','placeholder','SS0094_25_2','SS0094_26','SS0094_27','placeholder',...
    'placeholder','SS0094_28','SS0094_29','placeholder','placeholder','placeholder','placeholder','SS0094_30','SS0094_31','SS0094_32_2','SS0094_33','placeholder',...
    'SS0094_34','SS0094_35','SS0094_36','SS0094_37','SS0094_38','SS0094_39','SS0094_40'}; %beginning of filename raw files

scanID={'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s01
    '301','401','501','601','701','801'; ... %s02
    '401','501','601','701','801','901'; ... %s03
    '301','401','501','601','701','901'; ... %s04
    '301','401','501','601','701','801'; ... %s05
    '301','401','501','601','701','801'; ... %s06
    '301','501','701','801','901','1001'; ... %s07
    '301','401','501','601','701','801'; ... %s08
    '301','401','501','601','701','801'; ... %s09
    '301','401','501','601','701','801'; ... %s10
    '301','401','501','601','701','801'; ... %s11
    '301','401','501','601','701','801'; ... %s12
    '301','401','501','601','701','801'; ... %s13
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s14
    '301','401','501','601','701','801'; ... %s15
    '301','401','501','601','701','801'; ... %s16
    '301','501','601','701','801','901'; ... %s17
    '301','401','501','601','701','901'; ... %s18
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s19
    '401','501','601','701','801','901'; ... %s20
    '301','401','501','601','701','801'; ... %s21
    '301','401','501','601','701','801'; ... %s22
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s23
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s24
    '301','401','501','601','701','801'; ... %s25
    '301','401','501','601','701','801'; ... %s26
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s27
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s28
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s29
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s30
    '301','401','501','601','701','801'; ... %s31
    '301','401','501','601','701','801'; ... %s32
    '301','401','501','601','701','801'; ... %s33
    '301','401','501','601','701','801'; ... %s34
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s35
    '401','501','601','701','801','901'; ... %s36
    '301','401','501','601','701','801'; ... %s37
    '301','501','601','701','801','901'; ... %s38
    '301','401','501','601','701','801'; ... %s39
    '401','601','701','801','901','1001'; ... %s40
    '301','401','501','601','701','801'; ... %s41
    '301','401','501','601','701','801'; ... %s42
    };

fieldMapID={};
fieldMapID2={};
fieldMapName={'magnitude','phase'};
% runID={'R1','R2','R3','R4','R5','R6'};
run={'1','2','3','4','5','6'};
% dataext='run';
% hemName={'LeftHem','RightHem'};
% hem={'lh','rh'};
% atlasA={'i','x'};
% atlas= 2;
% atlasname={'fsaverage','fsaverage_sym'};


%% Parameters
delay=1; %SPM starts with TR 0; designfiles start counting with 1; -> has to be substracted out CHECK!!!
TR=2; %sec
nrVolumes=230; %number of volumes
sn8NrVolumes=218;
nrSlices=60; %different across subjects
radius=16; %mm
numVox=160;

%Slice Acquisition
MBsliceAcquisition = [1:nrSlices]';        %multiband acquisition
multiBand=2;       %N slices captured at one time
interSlice = TR/(nrSlices/multiBand);  %for multiband


%% Generate vector for multiband acquisition. Comment if using standard acquisition
%Default (interleaved) design
interSliceLoop = 0; %start at 0 for first slice
MBCol = MBsliceAcquisition(1:nrSlices/multiBand); %Create single column of slices
for i= 1:2:length(MBCol)    %loop through odd slices first
    MBCol(i) = interSliceLoop;
    interSliceLoop = interSliceLoop + interSlice;
end
for i=2:2:length(MBCol) %loop through even slices
    MBCol(i) = interSliceLoop;
    interSliceLoop = interSliceLoop + interSlice;
end
MBsliceAcquisition = repmat (MBCol, multiBand,1);    %replicate for all slice columns

sliceAcquisition = MBsliceAcquisition;


%Determine AC origin by hand Invert y and z coordinate
% originAC=[124, -131, -72]; %s01 (RY)
originAC=[1, 1, 1, 1; ... %s01 placeholder
    127,-131,-72, 1; ... %s02
    123,-129,-76, 1; ... %s03
    119,-131,-71, 1; ... %s04
    123,-135,-75, 1; ... %s05
    122,-131,-76, 1; ... %s06
    121,-129,-72, 1; ... %s07
    122,-138,-83, 1; ... %s08
    125,-141,-72, 1; ... %s09
    123,-130,-75, 1; ... %s10
    130,-138,-70, 1; ... %s11
    126,-130,-82, 1; ... %s12
    119,-132,-79, 1; ... %s13
    1, 1, 1, 1; ... %s14 placeholder
    124,-131,-83, 1; ... %s15
    121,-139,-76, 1; ... %s16
    120,-128,-74, 1; ... %s17
    127,-136,-74, 1; ... %s18
    1, 1, 1, 1; ... %s19 placeholder
    124,-135,-75, 1; ... %s20
    128,-131,-78, 1; ... %s21
    120,-130,-77, 1; ... %s22
    1, 1, 1, 1; ... %s23 placeholder
    1, 1, 1, 1; ... %s24 placeholder
    120,-132,-73, 1; ... %s25
    127,-130,-83, 1; ... %s26
    1, 1, 1, 1; ... %s27 placeholder
    1, 1, 1, 1; ... %s28 placeholder
    1, 1, 1, 1; ... %s29 placeholder
    1, 1, 1, 1; ... %s30 placeholder
    126,-132,-78, 1; ... %s31
    124,-131,-78, 1; ... %s32
    124,-135,-82, 1; ... %s33
    126,-135,-80, 1; ... %s34
    1, 1, 1, 1; ... %s35 placeholder
    121,-128,-81, 1; ... %s36
    122,-129,-78, 1; ... %s37
    127,-128,-81, 1; ... %s38
    125,-132,-73, 1; ... %s39
    124,-134,-70, 1; ... %s40
    124,-131,-77, 1; ... %s41
    125,-141,-84, 1; ... %s42
    ];

subcortStructs = ... %subcortical structures of interest for later analysis, defined from freesurfer aseg
    {'left_thalamus','left_caudate','left_putamen','left_pallidum','brain_stem','left_hippocampus',...
    'left_amygdala','left_accumbens_area','left_ventral_dc','left_vessel','left_choroid_plexus',...
    'right_thalamus','right_caudate','right_putamen','right_pallidum','right_hippocampus',...
    'right_amygdala','right_accumbens_area','right_ventral_dc','right_choroid_plexus'};
% 'right_vessel,' %reinsert before right_choroid_plexus

%%%
switch(what)
    
    case 'make_nii'
        sn=varargin{1};
        cd(fullfile(rawDir,subjID{sn}));
        
        %%Fieldmaps:
        try
            for i=1:2
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn}  '_' fieldMapID{sn,i} '_Fieldmap-'  num2str(i) '_te' fieldMapID2{sn,i} '.nii.gz']);
                dest = fullfile(rawDir, subjID{sn});
                gunzip(source,dest); %unzip
                
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn}  '_' fieldMapID{sn,i} '_Fieldmap-'  num2str(i) '_te' fieldMapID2{sn,i} '.nii']);
                dest   = fullfile(epiDir, subj_name{sn}, [subj_name{sn} '_' fieldMapName{sn,i} '.nii']); %EPI directory
                
                destFolder   = fullfile(epiDir, subj_name{sn}); %EPI directory
                k=exist(destFolder);
                if k==0
                    mkdir(destFolder)  ;
                end;
                
                movefile(source,dest); %move
                
            end;
            disp('fieldmap done')
        catch
            disp('No field map. Most likely a typo in file name or missing destination directory!');
        end;
        
        %%Anatomical:
        try
            %             source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} '_202_scombi.nii.gz']); %naming convention for s1-s16
            source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_202.nii.gz']); %naming convention changed on Odin from S17 onwards
            
            if sn == 20 %SS0094_25_2 had a different number assigned to scombi scan
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_302.nii.gz']); %naming convention changed on Odin from S17 onwards
            end
            
            if sn == 40 %SS0094_25_2 had a different number assigned to scombi scan
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_302.nii.gz']); %naming convention changed on Odin from S17 onwards
            end
            
            dest = fullfile(rawDir, subjID{sn});
            gunzip(source,dest); %unzip
            
            %             source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} '_202_scombi.nii']); %naming convention for s1-s16
            source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_202.nii']); %naming convention changed on Odin from S17 onwards
            
            if sn == 20 %SS0094_25_2 and SS0094_38 had a different number assigned to scombi scan
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_302.nii']); %naming convention changed on Odin from S17 onwards
            end
            
            if sn == 40 %SS0094_25_2 and SS0094_38 had a different number assigned to scombi scan
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_302.nii']); %naming convention changed on Odin from S17 onwards
            end
            
            dest   = fullfile(anatDir, subj_name{sn}, [subj_name{sn} '_anatomical_orig.nii']);
            
            destFolder   = fullfile(anatDir, subj_name{sn}); %EPI directory
            k=exist(destFolder);
            if k==0
                mkdir(destFolder)  ;
            end;
            
            movefile(source,dest); %move
            disp('anatomical done')
        catch
            disp('No anatomicals. Most likely a typo in file name or missing destination directory!');
        end;
        
        %%EPI:
        try
            for i=1:length(run) %for each of the 6 runs
                %                 source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} '_' scanID{sn,i} '_FE_EPI_2_iso_MB2' '.nii.gz']); %naming convention for s1-s16
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'FE_EPI_2_iso_MB21_' scanID{sn,i} '.nii.gz']); %naming convention changed on Odin from S17 onwards
                dest = fullfile(rawDir, subjID{sn});
                gunzip(source,dest); %unzip
                
                %                 source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} '_' scanID{sn,i} '_FE_EPI_2_iso_MB2' '.nii']); %naming convention for s1-s16
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'FE_EPI_2_iso_MB21_' scanID{sn,i} '.nii']); %naming convention changed on Odin from S17 onwards
                dest = fullfile(epiDir, subj_name{sn}, [subj_name{sn} '_run' run{i} '.nii']);
                
                destFolder   = fullfile(epiDir, subj_name{sn}); %EPI directory
                k=exist(destFolder);
                if k==0
                    mkdir(destFolder)  ;
                end;
                
                movefile(source,dest); %move
            end;
            disp('EPI done')
        catch
            disp('No EPI. Most likely a typo in file name or missing destination directory!');
        end;
        
    case 'set_AC' %reslice into LPI and set anterior commisure
        sn=varargin{1}; disp(subj_name{sn});
        cd(fullfile(anatDir,subj_name{sn}));
        V=spm_vol([subj_name{sn} '_anatomical_orig.nii']);
        X=spm_read_vols(V);
        V.mat(:,4)= originAC(sn,:); %put in AC coordinates (reverse the sign for each coordinate!)
        V.fname= [subj_name{sn} '_anatomical.nii'];
        spm_write_vol(V,X);
        
    case 'slice_timing'
        sn=varargin{1};
        cd(fullfile(epiDir, subj_name{sn}));
        for r= 1:length(run); %for each volume and each run, generate a path to the file and store it in N
            
            if sn == 8 && r == 4
                nrVolumes = sn8NrVolumes;
            else
                nrVolumes = nrVolumes;
            end
            
            for i=1:nrVolumes
                N{i,1} = [fullfile(epiDir, subj_name{sn}, [ subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
            end;
            J.scans{r} = N;
        end
        matlabbatch{1}.spm.temporal.st.scans = J.scans; %provide all the details required to the 'matlabbatch' variable ready for the function
        matlabbatch{1}.spm.temporal.st.nslices = nrSlices;
        matlabbatch{1}.spm.temporal.st.tr = TR;
        matlabbatch{1}.spm.temporal.st.ta = TR-(TR/nrSlices);
        matlabbatch{1}.spm.temporal.st.so = (sliceAcquisition);
        matlabbatch{1}.spm.temporal.st.refslice = 1;
        matlabbatch{1}.spm.temporal.st.prefix = 'a';
        spm_jobman('run',matlabbatch);
        
    case 'realign_unwarp' %extracts motion regressors for glm; http://www.diedrichsenlab.org/imaging/robustWLS.html
        sn=varargin{1};
        cd(fullfile(epiDir, subj_name{sn}));
        prefix='a';
        spmj_realign_unwarp(epiDir, subj_name{sn}, run, 1, nrVolumes,'prefix',prefix)
        
    case 'meanepi'
        sn=varargin{1};
        
        cd(fullfile(epiDir, subj_name{sn}));
        
        matlabbatch{1}.spm.util.imcalc.input = {
            ['ua', subj_name{sn}, '_run1.nii,1']
            ['ua', subj_name{sn}, '_run2.nii,1']
            ['ua', subj_name{sn}, '_run3.nii,1']
            ['ua', subj_name{sn}, '_run4.nii,1']
            ['ua', subj_name{sn}, '_run5.nii,1']
            ['ua', subj_name{sn}, '_run6.nii,1']
            };
        matlabbatch{1}.spm.util.imcalc.output = ['meanepi_',subj_name{sn},'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {''};
        matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6)/6';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        
        spm_jobman('run',matlabbatch);
        
    case 'coreg'
        
        sn=varargin{1}; disp(subj_name{sn});
        
        %% *1. Coregister meanEPI to anatomical (no reslicing into anat space, dimension of EPI preserved!): %%%
        %         NOTE: If original image completely off in terms of alignment to anat - First coregister per hand via 'coregtool' command
        %         then run the spm algorithm:
        %         ref=fullfile(anatDir, subj_name{sn},[subj_name{sn}, '_anatomical.nii']);
        %         source=fullfile(epiDir, subj_name{sn},['meanepi_' subj_name{sn} '.nii']);
        %         matlabbatch{1}.spm.spatial.coreg.estimate.ref =  {ref};
        %         matlabbatch{1}.spm.spatial.coreg.estimate.source = {source};
        %         matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        %         matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        %         matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        %         matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        %         matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        %         spm_jobman('run',matlabbatch); %meanepi will be overwritten by a version that is corgegistered
        
        %%% *2. Check with coregtool (first input anat, second meanepi), adjust and/or repeat as necessary until sulci aligned
        %%make sure you overwrite and work with the correct meanepi* files
        %                 coregtool;
        
        %%% *3. Coregister EPI runs to meanEPI: BACK UP ua*_run* epi data
        %%% FIRST for each subject.
        
        %         source=fullfile(epiDir, subj_name{sn},['meanepi_' subj_name{sn} '.nii']);
        %         for r=1:length(run)
        %             Epi = fullfile(epiDir,subj_name{sn}, ['ua' subj_name{sn},'_run',run{r},'.nii']);
        %             spmj_makesamealign_nifti(source,Epi);
        %         end;
        %
        %%% *4. Check whether runs are really aligned to the rmeanepi on
        %%% MRIcron: using anatomical and ua* data overlap
        
    case 'glm_set'
        sn=varargin{1};
        prefix='ua'; %% if bias correction was performed previously
        delay=0;
        dur=nrVolumes;  % TRs per Block
        fMRIblockNr=(47:52); %BN behavioural file
        T=[];
        
        for r=1:numel(run)
            
            subjName=subj_name;
            
            %Load behavioural file (structure B)
            fname=fullfile(behDir,subj_name{sn},['exp_BN' num2str(fMRIblockNr(r)) '.mat']);
            load(fname);
            
            subj_name=subjName;
            
            %Compute condition index
            B.cond=nan(size(B.seqID,1),1);
            B.error=(B.points==0);
            B.goCueDur=nan(size(B.trialType,1),1);
            B.goCueDur(B.trialType==1,1)=4000; %CHECK trial type!
            B.goCueDur(B.trialType==2,1)=1000;
            
            c=0;
            for seqID=1:4
                for trialType=1:2 %1: seq1,prod;  2: seq1,prep; 3: seq2,prod; 4: seq2,prep etc.
                    c=c+1;
                    B.cond(B.seqID==seqID & B.trialType==trialType,1)=c;
                end;
            end;
            
            %General info GLM
            J.dir = {fullfile(glmDir, subj_name{sn})}; % Working directory of the GLM
            J.timing.units = 'secs';                   % Units (scans or sec for ons)
            J.timing.RT = TR;                     % TR
            J.timing.fmri_t = nrSlices;                   % Mircotime resolution
            J.timing.fmri_t0 = 1; % reference slice
            
            %Gather info for each run
            %%% Load all images/TRs (N)
            for i=1:nrVolumes                 % Loop over all images and put into the data structure
                N{i,1} = [fullfile(epiDir, subj_name{sn}, [prefix subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
                %                 N{i,1} = [fullfile(epiDir, subj_name{sn}, [subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
            end;
            
            J.sess(r).scans= N;         % Number of images
            
            c=[];
            
            
            
            %%% Assign Conditions
            %%% Sequence cue
            for c=2:2:8                  % Loop over the Prep conditions (catch trials only) cond: 2     4     6     8
                idx=find(B.cond==c);
                
                % Make new structure(S) that reflects what each of the
                % regressors in the design matrix -> goes into T
                S.RN(c,1)=r; %run n
                S.SN(c,1)=sn; %subject nr
                S.condition(c,1)=c;
                S.numtrials(c,1)=length(idx);
                
                % Info on condition name, onset and duration
                J.sess(r).cond(c).name = sprintf('%d_%d',r,c);
                J.sess(r).cond(c).onset = [B.tZero(idx)]/1000;
                %                 J.sess(r).cond(c).duration = 0; %does not set a duration for preparation
                J.sess(r).cond(c).duration = max(B.tZero- B.tZeroCue)/1000; %max prep time in sec; %[(B.tZero(idx) - B.tZeroCue(idx))/1000];
                J.sess(r).cond(c).tmod = 0;
                J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
            end;
            
            
            
            %%% Sequence initiation/1st press
            for c=1:2:8                  % Loop over the conditions (production trials only!!) cond: 1     3     5     7
                idx=find(B.cond==c);
                
                nanidx = ~isnan(B.timing(:,1)); %ignore any production trials where there were no presses (all timings were NaN)
                nanidx = nanidx(idx);
                idx = idx(nanidx);
                
                % Make new structure(S) that reflects what each of the
                % regressors in the design matrix -> goes into T
                S.RN(c,1)=r; %run n
                S.SN(c,1)=sn; %subject nr
                S.condition(c,1)=c;
                S.numtrials(c,1)=length(idx);
                
                % Info on condition name, onset and duration
                J.sess(r).cond(c).name = sprintf('%d_%d',r,c);
                J.sess(r).cond(c).onset = [(B.tZero(idx)+ B.timing(idx,1))/1000];
                J.sess(r).cond(c).duration = 0;
                %                 J.sess(r).cond(c).duration = [(max(B.timing(idx),[],2)-B.timing(idx,1))/1000];
                %                 J.sess(r).cond(c).duration = [(B.timing(idx,5)-B.timing(idx,1))/1000];
                J.sess(r).cond(c).tmod = 0;
                J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
            end;
            
            
            
            
            %%% Error (feedback)
            c=9;
            idx=[];
            idx=find(B.error==1);
            
            % Make new structure(S) that reflects what each of the
            % regressors in the design matrix -> goes into T
            S.RN(c,1)=r; %run n
            S.SN(c,1)=sn; %subject nr
            S.condition(c,1)=c;
            S.numtrials(c,1)=length(idx);
            
            % Info on condition name, onset and duration
            J.sess(r).cond(c).name = sprintf('%d_%d',r,c);
            
            k=isempty(idx);
            if k==0
                J.sess(r).cond(c).onset = [B.tZeroCue(idx)/1000];
                J.sess(r).cond(c).duration = [( (B.tZero(idx)- B.tZeroCue(idx))+B.goCueDur(idx)+B.crossDur(idx)+B.feedbackCalcDur(idx) + 1000 )/1000 ]; %%%TODO
                
            else
                J.sess(r).cond(c).onset = -99;
                J.sess(r).cond(c).duration = 1;
                
            end;
            
            J.sess(r).cond(c).tmod = 0;
            J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
            
            c=[];
            
            
            
            %%% Separate Prep regressor during production to explain the
            %%% prep variance in production trials
            vectorProdTrials=(1:2:8);
            % for c=1:2:8                  % Loop over the conditions (production trials only!!) conf: 1     3     5     7
            idx=find(ismember(B.cond,vectorProdTrials));
            cond=10;
            % Make new structure(S) that reflects what each of the
            % regressors in the design matrix -> goes into T
            S.RN(c,1)=r; %run n
            S.SN(c,1)=sn; %subject nr
            S.condition(c,1)=cond;
            S.numtrials(c,1)=length(idx);
            
            % Info on condition name, onset and duration
            J.sess(r).cond(cond).name = sprintf('%d_%d',r,cond);
            J.sess(r).cond(cond).onset = [B.tZeroCue(idx)/1000];
            J.sess(r).cond(cond).duration = [(B.tZero(idx) - B.tZeroCue(idx))/1000]; %Time between Sequence and Go cues
            J.sess(r).cond(cond).tmod = 0;
            J.sess(r).cond(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
            %end;
            c=cond;
            
            
            
            %             %%% Visual regressor (fractal image)
            %             c=11;
            %             idx= 1:length(B.cond);
            %
            %             % Make new structure(S) that reflects what each of the
            %             % regressors in the design matrix -> goes into T
            %             S.RN(c,1)=r; %run n
            %             S.SN(c,1)=sn; %subject nr
            %             S.condition(c,1)=c;
            %             S.numtrials(c,1)=length(idx);
            %
            %             % Info on condition name, onset and duration
            %             J.sess(r).cond(c).name = sprintf('%d_%d',r,c);
            %             J.sess(r).cond(c).onset = [B.tZeroCue(idx)/1000];
            %             J.sess(r).cond(c).duration = [0.4];
            %             J.sess(r).cond(c).tmod = 0;
            %             J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
            
            
            
            %             %%% Points-based rumination regressor
            %             c=11;
            %             idx= 1:length(B.cond);
            %
            %             % Make new structure(S) that reflects what each of the
            %             % regressors in the design matrix -> goes into T
            %             S.RN(c,1)=r; %run n
            %             S.SN(c,1)=sn; %subject nr
            %             S.condition(c,1)=c;
            %             S.numtrials(c,1)=length(idx);
            %
            %             % Info on condition name, onset and duration
            %             J.sess(r).cond(c).name = sprintf('%d_%d',r,c);
            %             J.sess(r).cond(c).onset = [(B.tZeroCue(idx) + B.goCueDur(idx) + B.crossDur(idx) + B.feedbackCalcDur(idx))/1000];
            %             J.sess(r).cond(c).duration = [(1000 + B.iti(idx))/1000];
            %             J.sess(r).cond(c).tmod = 0;
            %
            %             %Parametric modulation for points-based regressor
            %             currentBlock = B.points;
            %
            %             fname=fullfile(behDir,subj_name{sn},['exp_BN' num2str(fMRIblockNr(r)-1) '.mat']); %load block previous to current block, to allow for z score transformation based on previous scores
            %             prevBlock = load(fname); prevBlock = prevBlock.B.points;
            % %             currentBlock = rmfield (currentBlock,{'cond', 'error', 'goCueDur'});
            % %             currentBlock.error = [];
            % %             currentBlock.goCueDur = [];
            %
            %             rawPoints = [prevBlock; currentBlock]; %concatenate previous with current data
            %
            % %             rawPoints = addstruct(prevBlock, currentBlock); %concatenate previous with current data
            %             zScoreConcat = NaN(1,48); %variable to store zScores and be wiped before each pass %%%FUTRE-PROOF%%%
            %
            %             for i=length(prevBlock)+1 : length(rawPoints) %produce z-scored points for current block, using previous blocks to begin value calculations
            %
            %                 window = i-10 : i-1; %adjust window for STDev calcs as needed
            %
            %                 meanPoints = mean(rawPoints(window)); %calc the mean and STDev of the window around point of interest
            %                 stDevPoints = std(rawPoints(window));
            %
            %
            %                 zScore = (rawPoints(i) - meanPoints) / stDevPoints;
            %
            %                 if stDevPoints == 0
            %                     zScore = abs(rawPoints(i) - meanPoints);
            %                 elseif isinf(zScore)
            %                     zScore = abs(rawPoints(i) - meanPoints);
            %                 end
            %
            %
            %                 zScoreConcat(i-length(prevBlock)) = zScore;
            %             end
            %
            %             J.sess(r).cond(c).pmod = struct('name', {'feedback'}, 'param', {zScoreConcat}, 'poly', {});
            % %             J.sess(r).cond(c).pmod.name = 'feedback';
            % %             J.sess(r).cond(c).pmod.param = zScoreConcat;
            % %             J.sess(r).cond(c).pmod.poly = 1;
            
            
            
            
            
            % Make new structure(S) that reflects what each of the
            % regressors in the design matrix -> goes into T
            idx=sum(S.RN); %%%
            S.RN(c,1)=r; %run n
            S.SN(c,1)=sn; %subject nr
            S.condition(c,1)=c;
            S.numtrials(c,1)=idx; %Nr of trials in a run
            
            J.sess(r).multi = {''};
            J.sess(r).regress = struct('name', {}, 'val', {});
            J.sess(r).multi_reg = {''};
            J.sess(r).hpf = 128;
        end;
        
        
        T=addstruct(T,S);
        %J.fact = struct('name', {}, 'levels', {});
        
        %%%Derivatives:
        J.bases.hrf.derivs = [1 0];
        J.volt = 1;
        J.global = 'None';
        
        %%J.mask = {[fullfile(anatDir,  subj_name{sn},[subj_name{sn},'_mask4glm.img'])],1};
        J.mthresh = 0.8;
        J.mask = {''};
        
        %             J.cvi =  'AR(1)';
        %             matlabbatch{1}.spm.stats.fmri_spec=J; %SPM standard
        J.cvi =  'wls';
        matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec=J; %joern's toolbox http://www.diedrichsenlab.org/imaging/robustWLS_spm12.html; code updates as documented inline
        
        spm_jobman('run',matlabbatch);
        dsave(fullfile(glmDir, subj_name{sn},'SPM_info.ana'),T);  % save the information file for the SPM
        disp('GLM set DONE.');
        
        
        %%%CHECK: load('SPM.mat'); imagesc(SPM.xX.X)    ; use plot as well
        %http://people.duke.edu/~njs28/spmdatastructure.htm
        
    case 'glm_estimate'
        %%% Estimate GLM
        sn=varargin{1};
        fname=fullfile(glmDir,subj_name{sn},['SPM.mat']);
        
        
        matlabbatch{1}.spm.tools.rwls.fmri_rwls_est.spmmat = { fname}; %joern's toolbox; TRs with motion artefact weighted down in the GLM; code update as documented inline
        matlabbatch{1}.spm.tools.rwls.fmri_rwls_est.method.Classical = 1;
        
        
        spm_jobman('run',matlabbatch);
        disp('GLM estimate DONE.');
        
    case 'glm_contrast'
        
        sn=varargin{1};
        cd(fullfile(glmDir, subj_name{sn}));
        load SPM;
        SPM=rmfield(SPM,'xCon'); %remove old contrasts
        
        %1: seq1,prod; 2: deriv
        %3: seq1,prep; 4; deriv
        %5: seq2,prod; 6; deriv
        %7: seq2,prep  8; deriv
        %9: seq3,prod; 10; deriv
        %11: seq3,prep; 12; deriv
        %13: seq4,prod; 14; deriv
        %15: seq4,prep;  16; deriv
        %17: error/feedback  18; deriv
        %19: prep in prod trials 20; deriv
        
        % % % % % % % % % % % % % % % % % % % % % % %         %21: points ; 22% deriv %not currently included.
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        
        %% Define contrasts:
        %%% Production
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,numel(run)) runBSL];
        prod=prod/sum(prod); %make the sum of the vector 1
        SPM.xCon(1)=spm_FcUtil('Set','Move', 'T', 'c',prod',SPM.xX.xKXs);
        
        %%% Preparation
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,numel(run)) runBSL];
        prep=prep/sum(prep);
        SPM.xCon(2)=spm_FcUtil('Set','Prep', 'T', 'c',prep',SPM.xX.xKXs);
        
        %%% Error
        feedbackE =[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]; %error feedback
        feedbackE=[repmat(feedbackE,1,numel(run)) runBSL];
        feedbackE=feedbackE/sum(feedbackE);
        SPM.xCon(3)=spm_FcUtil('Set','FeedbackE', 'T', 'c',feedbackE',SPM.xX.xKXs);
        
        %%% Prep vs Prod
        prepVSprod =[-1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0 0 0 0 0]; %PrepVProd
        prepVSprod=[repmat(prepVSprod,1,numel(run)) runBSL];
        %prepVSprod=prepVSprod/sum(prepVSprod);
        SPM.xCon(4)=spm_FcUtil('Set','PrepVProd', 'T', 'c',prepVSprod',SPM.xX.xKXs);
        
        %%% Prod vs Prep
        prodVSprep =[1 0 -1 0 1 0 -1 0 1 0 -1 0 1 0 -1 0 0 0 0 0]; %ProdVPrep
        prodVSprep=[repmat(prodVSprep,1,numel(run)) runBSL];
        %prepVSprod=prepVSprod/sum(prepVSprod);
        SPM.xCon(5)=spm_FcUtil('Set','ProdVPrep', 'T', 'c',prodVSprep',SPM.xX.xKXs);
        
        %%% Rest
        rest =[-1 0 -1 0 -1 0 -1 0 -1 0 -1 0 -1 0 -1 0 -1 0 -1 0]; %rest
        rest=[repmat(rest,1,numel(run)) runBSL];
        rest=rest/sum(abs(rest));
        SPM.xCon(6)=spm_FcUtil('Set','Rest', 'T', 'c',rest',SPM.xX.xKXs);
        
        %% Do the contrasts
        SPM=spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save SPM SPM;
        
        %    SPM.xCon(1)=spm_FcUtil('Set','R_move', 'T', 'c',con,SPM.xX.xKXs);
        %
        %             %-----Error
        %             con=zeros(length(D.condition),1); %define new contrast vector
        %             con(D.condition==10,1)=1; %choose the correct conditions
        %             con=con/sum(con); %so that it adds up to 1
        %             con=[con;rest]; %add rest regressors
        %             SPM.xCon(2)=spm_FcUtil('Set','Error', 'T', 'c',con,SPM.xX.xKXs);
        %
        %             %____do the constrasts
        %             SPM=spm_contrasts(SPM,[1:length(SPM.xCon)]);
        %             save SPM SPM;
        
    case 'con_smooth'
        
        
        s=varargin{1};
        
        con=fullfile(glmDir, subj_name{s},'con_0001.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0001.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(glmDir, subj_name{s},'con_0002.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0002.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(glmDir, subj_name{s},'con_0003.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0003.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(glmDir, subj_name{s},'con_0004.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0004.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(glmDir, subj_name{s},'con_0005.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0005.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(glmDir, subj_name{s},'con_0006.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0006.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
    case 'glm_contrastGroup'
        
        %Produces second-level contrasts. Edit contrast folders, image names, and subject nifti files.
        
        dataDir = {'Mov', 'Prep', 'Error', 'PrepProd', 'ProdPrep', 'Rest'}; %%Save folders for each contrast
        images = {'scon_0001';'scon_0002';'scon_0003';'scon_0004';'scon_0005';'scon_0006'};
        
        %         dataDir = {'Mov','Prep','Points'}; %%Save folders for each contrast
        %         images = {'scon_0001';'scon_0002';'scon_0006'};
        
        %         dataDir = {'Prep'}; %%Save folders for each contrast
        %         images = {'scon_0002'};
        
        subNii = {'_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s09.nii','_s10.nii','_s13.nii','_s16.nii','_s17.nii','_s18.nii','_s20.nii','_s21.nii','_s22.nii','_s25.nii','_s26.nii'...
            '_s31.nii','_s32.nii','_s34.nii','_s36.nii','_s38.nii','_s39.nii','_s40.nii','_s41.nii','_s42.nii'};
        contrastN = length(dataDir);
        images = repmat(images,1,length(subNii));
        subNii = repmat (subNii,length(dataDir),1);
        fileName = strcat (images,subNii);  %%Concatenate contrast files and subject names
        
        for i=1:contrastN  %%Loop across contrasts, plugging parameters into SPM.
            glmscndDir = fullfile(scndDir, dataDir(i));
            matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (groupDir, fileName(i,:))';  %%Select files from vectors above
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            spm_jobman('run',matlabbatch);  %Run SPM
        end
        
    case 'glm_contrastEstimate'
        
        dataDir = {'Mov', 'Prep', 'Error', 'PrepProd', 'ProdPrep', 'Rest'}; %%Save folders for each contrast
        contrastN = length(dataDir);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat = fullfile (scndDir, dataDir(i), 'SPM.mat');
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
    case 'segment' %Segment grey matter + white matter + csf
        sn=varargin{1};
        
        J.data = {fullfile(anatDir, subj_name{sn},[subj_name{sn} '_anatomical.nii,1'])};
        J.output.GM = [0 0 1];
        J.output.WM = [0 0 1];
        J.output.CSF = [0 0 0];
        J.output.biascor = 1;
        J.output.cleanup = 0;
        J.opts.tpm = {
            'D:\projects\toolboxes\spm12\tpm\OldSeg\grey.nii'
            'D:\projects\toolboxes\spm12\tpm\OldSeg\white.nii'
            'D:\projects\toolboxes\spm12\tpm\OldSeg\csf.nii'
            };
        J.opts.ngaus = [2
            2
            2
            4];
        J.opts.regtype = 'mni';
        J.opts.warpreg = 1;
        J.opts.warpco = 25;
        J.opts.biasreg = 0.0001;
        J.opts.biasfwhm = 60;
        J.opts.samp = 3;
        J.opts.msk = {''};
        matlabbatch{1}.spm.spatial.preproc=J;
        spm_jobman('run',matlabbatch);
        
    case 'make_mask'  % Makes restricted analysis mask for MVA
        s=varargin{1};
        mask=fullfile(glmDir, subj_name{s},'mask.nii');
        omask=fullfile(glmDir, subj_name{s},'maskbrain.nii'); %output mask to be used in the future
        P1=fullfile(anatDir,subj_name{s},sprintf('c1%s_anatomical.nii',subj_name{s}));
        P2=fullfile(anatDir,subj_name{s},sprintf('c2%s_anatomical.nii',subj_name{s}));
        sP1=fullfile(anatDir,subj_name{s},sprintf('sc1%s_anatomical.nii',subj_name{s}));
        sP2=fullfile(anatDir,subj_name{s},sprintf('sc2%s_anatomical.nii',subj_name{s}));
        spm_smooth(P1,sP1,[4 4 4]); %smooth with 4mm kernel
        spm_smooth(P2,sP2,[4 4 4]);
        spm_imcalc_ui({mask,sP1,sP2},omask,'i1 & (i2 +i3)>0.01',{}); %recorded activity in brain (grey + white matter)
        
    case 'MVA_search' % Define the search lights for the MVA analysis
        
        s=varargin{1};
        radius=16;
        numVox=160;
        cd(fullfile(glmDir, subj_name{s}));
        V=spm_vol('maskbrain.nii'); %if preceded by case MVA_mask
        X=spm_read_vols(V);
        [i,j,k]=ind2sub(size(X),find(X~=0));
        vox=[i j k];
        [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[radius numVox],V.mat,V.dim,[],'mva160_numvox.nii');
        save volsearch160.mat vox LI voxmin voxmax n
        
    case 'MVA_do_overallMov'                 % Conduct the classification analysis 4 sequences
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Comb_160_Mov.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prod>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_do_overallPrep'                 % Conduct the classification analysis 4 sequences
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Comb_160_Prep.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_do_spatOneout_Mov'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Spat_160_Mov.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 0; 0 1];
            oneout=repmat(oneout,1,12);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_spatOneout_Prep'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Spat_160_Prep.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 0; 0 1];
            oneout=repmat(oneout,1,12);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_tempOneout_Mov'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Temp_160_Mov.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 1 0 0; 0 0 1 1];
            oneout=repmat(oneout,1,6);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_tempOneout_Prep'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Temp_160_Prep.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 1 0 0; 0 0 1 1];
            oneout=repmat(oneout,1,6);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_Int_Mov'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Int_160_Mov.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prod>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_do_Int_Prep'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Int_160_Prep.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'MVA_do_overallMov_points'                 % Conduct the classification analysis 4 sequences
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Comb_160_Mov.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prod>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_do_overallPrep_points'                 % Conduct the classification analysis 4 sequences
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Comb_160_Prep.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_do_spatOneout_Mov_points'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Spat_160_Mov.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 0; 0 1];
            oneout=repmat(oneout,1,12);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_spatOneout_Prep_points'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Spat_160_Prep.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 0; 0 1];
            oneout=repmat(oneout,1,12);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_tempOneout_Mov_points'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Temp_160_Mov.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 1 0 0; 0 0 1 1];
            oneout=repmat(oneout,1,6);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_tempOneout_Prep_points'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Temp_160_Prep.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 1 0 0; 0 0 1 1];
            oneout=repmat(oneout,1,6);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_Int_Mov_points'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Int_160_Mov.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prod>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_do_Int_Prep_points'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Int_160_Prep.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_group'
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
        images = {'szacc_Comb_160_Mov';'szacc_Comb_160_Prep';'szacc_Int_160_Mov';'szacc_Int_160_Prep';'szacc_Spat_160_Mov';'szacc_Spat_160_Prep';'szacc_Temp_160_Mov';'szacc_Temp_160_Prep'};
        
        %         dataDir = {'MVA_comb_mov', 'MVA_comb_prep'}; %%Save folders for each contrast
        %         images = {'szacc_Comb_160_Mov';'szacc_Comb_160_Prep'};
        
        %ONLY PARTICIPANTS WHO MODULATED TIMING
        subNii = {'_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s09.nii','_s10.nii','_s13.nii','_s16.nii','_s17.nii','_s18.nii','_s20.nii','_s21.nii','_s22.nii','_s25.nii','_s26.nii'...
            '_s31.nii','_s32.nii','_s34.nii','_s36.nii','_s38.nii','_s39.nii','_s40.nii','_s41.nii','_s42.nii'}; %3 5 6 7 9 10 13 16 17 18 20 21 22 25 26 31 32 34 36 38 39 40 41 42
        
        %ALL PARTICIPANTS regardless of timing modulation ***TO DO***
        %         subNii = {'_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s09.nii','_s10.nii','_s13.nii','_s16.nii','_s17.nii','_s18.nii','_s20.nii','_s21.nii','_s22.nii','_s25.nii','_s26.nii'...
        %             '_s31.nii','_s32.nii','_s34.nii','_s36.nii','_s38.nii','_s39.nii','_s40.nii','_s41.nii'}; %3 4 5 6 7 9 10 11 12 13 15 16 17 18 20 21 22 25 26 31 32 33 34 36 37 38 39 40 41 42
        
        contrastN = length(dataDir);
        images = repmat(images,1,length(subNii));
        subNii = repmat (subNii,length(dataDir),1);
        fileName = strcat (images,subNii);  %%Concatenate contrast files and subject names
        
        for i=1:contrastN  %%Loop across contrasts
            glmscndDir = fullfile(groupDir, dataDir(i));
            matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;  %Adjust directory
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (groupDir, fileName(i,:))';  %%Select files from matrix
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
        %AVERAGE GROUP DATA
        %open spm fmri and select 'imcalc'
        %choose all images from participants of interest
        %(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+16+i17+i18+i19+i20+i21+i22+i23+i24)/24
        
    case 'MVA_estimate'
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
        %         dataDir = {'MVA_comb_mov', 'MVA_comb_prep'}; %%Save folders for each contrast
        contrastN = length(dataDir);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat = fullfile (groupDir, dataDir(i), 'SPM.mat');  %Adjust directory
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
    case 'MVA_zValue'
        
        s=varargin{1};
        cd(fullfile(glmDir,subj_name{s}));
        
        numTests=6;
        numCat=4;
        mu=1/numCat; %mu=0.25;
        N=numTests*numCat;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        images= {'_accuracy_Comb_160_Mov','_accuracy_Comb_160_Prep','_accuracy_Int_160_Mov','_accuracy_Int_160_Prep'};
        %         images= {'_accuracy_Comb_160_Mov','_accuracy_Comb_160_Prep'};
        
        outimages={'_zacc_Comb_160_Mov','_zacc_Comb_160_Prep','_zacc_Int_160_Mov','_zacc_Int_160_Prep'};
        %         outimages={'_zacc_Comb_160_Mov','_zacc_Comb_160_Prep'};
        
        
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},[subj_name{s} images{j} '.nii']);
            output_image= fullfile(glmDir,subj_name{s},[subj_name{s} outimages{j} '.nii']);
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
        end;
        
    case 'MVA_zValue_oneOut'
        
        s=varargin{1};
        cd(fullfile(glmDir,subj_name{s}));
        
        takeOneOutIter=2;
        numTests=6;
        numCat=2;
        mu=1/numCat; %mu=0.5;
        N=numTests*numCat*takeOneOutIter;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        images= {'_accuracy_Spat_160_Mov','_accuracy_Spat_160_Prep','_accuracy_Temp_160_Mov','_accuracy_Temp_160_Prep'};
        
        outimages={'_zacc_Spat_160_Mov','_zacc_Spat_160_Prep','_zacc_Temp_160_Mov','_zacc_Temp_160_Prep'};
        
        
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},[subj_name{s} images{j} '.nii']);
            output_image= fullfile(glmDir,subj_name{s},[subj_name{s} outimages{j} '.nii']);
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
        end;
        
    case 'MVA_smooth' %%%Smoothing in subject space
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Comb_160_Mov.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Comb_160_Mov.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Comb_160_Prep.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Comb_160_Prep.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Spat_160_Mov.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Spat_160_Mov.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Spat_160_Prep.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Spat_160_Prep.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Temp_160_Mov.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Temp_160_Mov.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Temp_160_Prep.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Temp_160_Prep.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Int_160_Mov.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Int_160_Mov.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Int_160_Prep.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Int_160_Prep.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        %smooth other images here as required
        
    case 'MNI_normalization'
        mkdir(fullfile(groupDir,'data')); % folder for each contrast
        
        
        %MVPA accuracy maps
        images= {'szacc_Int_160_Mov.nii','szacc_Int_160_Prep.nii'...
            'szacc_Comb_160_Mov.nii','szacc_Comb_160_Prep.nii',...
            'szacc_Spat_160_Mov.nii','szacc_Spat_160_Prep.nii','szacc_Temp_160_Mov.nii','szacc_Temp_160_Prep.nii'}; % please add other images as required, e.g. spmT_...
        %
        %
        % 'scon_0001.nii','scon_0002.nii','scon_0003.nii','scon_0004.nii','scon_0005.nii','scon_0006.nii'
        
        
        
        %         images= {'szacc_Comb_160_Mov.nii','szacc_Comb_160_Prep.nii',...
        %             'scon_0001.nii','scon_0002.nii','scon_0003.nii','scon_0004.nii','scon_0005.nii','scon_0006.nii'}; % please add other images as required, e.g. spmT_...
        
        for i=1:length(images)
            anaImages{1} = images{i};
            s=varargin{1};
            defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
            for j=1:numel(anaImages)
                [~,name,ext]=spm_fileparts(anaImages{j});
                sn_images{j}= fullfile(glmDir,subj_name{s},[subj_name{s} '_' anaImages{j}]);
                out_images{j}= fullfile(groupDir,[name '_' subj_name{s} '.nii']);
                spmj_normalization_write(defor, sn_images,'outimages',out_images); %Trilinear interpolation
            end
        end
        
        
        %         images= {'szacc_Int_160_Prep.nii'}; % please add other images as required, e.g. spmT_...
        %
        %         s=varargin{1};
        %         defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
        %         for j=1:numel(images)
        %             [dir,name,ext]=spm_fileparts(images{j});
        %             sn_images{j}= fullfile(glmDir,subj_name{s},[subj_name{s} '_' images{j}]);
        %             out_images{j}= fullfile(groupDir,[name '_' subj_name{s} '.nii']);
        %             spmj_normalization_write(defor, sn_images,'outimages',out_images); %Trilinear interpolation
        %         end
        
        %Contrasts
        %                 images= {'con_0001.nii','con_0002.nii'}; % please add other images as required, e.g. spmT_...
        %                 images= {'scon_0001.nii','scon_0002.nii','scon_0003.nii','scon_0004.nii','scon_0005.nii','scon_0006.nii'}; % please add other images as required, e.g. spmT_...
        
        %          images= {'szacc_Comb_160_Prep.nii'}; % please add other images as required, e.g. spmT_...
        
        %                 s=varargin{1};
        %                 defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
        %                 for j=1:numel(images)
        %                     [dir,name,ext]=spm_fileparts(images{j});
        %                     sn_images{j}= fullfile(glmDir,subj_name{s},[images{j}]);
        %                     out_images{j}= fullfile(groupDir,[name '_' subj_name{s} '.nii']);
        %                 end
        %                 spmj_normalization_write(defor, sn_images,'outimages',out_images); %Trilinear interpolation
        
    case 'ROI_define' %%% 1. Define ROI via GUI:
        R=region_getregions; %opens gui
        
        %%% from Elife 2014:
        %Overall-cM1; %Overall-iPMd ( x coordinate flipped to left);
        %%Overall_cSMA; %Overall_cSupParietal
        %         roiOUT='D:\projects\rhys\prepProd\data\imaging\ROI\roi_Elife2014.mat';
        %             %%% Select/Copy paste parameters as follows:
        % % %         mni= [-36 -22 53; -21 -12 60; -8 -12 57; -32 -54 56];
        % % %         radius=6; %mm
        % % %         name={'L_M1','L_PMd','L_SMA','L_SPC'};
        
        %%% from Elife 2014 bilateral:
        roiOUT='E:\projects\rhys\prepProd2\data\imaging\ROI\roi_Elife2014_bilateral.mat';
        
        
        %%% Select/Copy paste parameters as follows:
        %         point
        mni=[-36 -22 53 ; 36 -22 53 ; -21 -12 60  ; 21 -12 60; -8 12 57;  8 18 49; -32 -54 56; 30 -59 46]; %Overall
        
        % mni=[-33 -23 59 ; 33 -23 59 ; -31 -13 53  ; 31 -13 53; -9 1 54; -9 1 54; -32 -54 56; 32 -54 56]; %from Integrated & Spatial
        %
        %         ROIname={'L_M1','R_M1','L_PMd','R_PMd','L_SMA','R_SMA','L_SPC','R_SPC'};
        %         click *done
        
        %%% Then save:
        save(roiOUT,'R');
        R=region_saveasimg(R,'avg152T1.nii','name','ROI');
        
    case 'ROI_run'        %%% 2. Extract values via GUI
        %         roiIN='D:\projects\rhys\prepProd\data\imaging\ROI\roi_Elife2014.mat';
        roiIN='E:\projects\rhys\prepProd2\data\imaging\ROI\roi_Elife2014_bilateral.mat';
        
        load(roiIN);
        P={};
        %%%% Comb
        cd('E:\projects\rhys\prepProd2\data\imaging\ROI');
        MVA_ID={'Comb','Temp','Spat','Int'};%%%TODO
        Phase_ID={'Prep','Mov'};
        ROI=[];
        K=[];
        
        subjects=anaSubj;
        loopCounter=1;
        subj_name = subj_name;
        
        for phase=1:2 %prep and production
            for m=1:numel(MVA_ID) %loop for MVA type
                %                 for s=1:numel(subj_name)
                for s=subjects
                    %                P{s}=fullfile(groupDir, ['szacc_' MVA_ID{m} '_160_Mov_' subj_name{s}, '.nii']);
                    
                    P{loopCounter}=fullfile(groupDir, ['szacc_' MVA_ID{m} '_160_' Phase_ID{phase} '_' subj_name{s}, '.nii']);
                    
                    loopCounter=loopCounter+1;
                end;
                loopCounter=1;
                V=spm_vol(char(P));
                stats = region_getdata(V,R);
                D=[];
                S=[];
                
                for i=1:size(stats,2) %loop over ROIs
                    meanROI(:,i)=nanmean(stats{i},2); %average across voxels
                    D.meanROI=meanROI(:,i);
                    D.subj=subj(anaSubj)';
                    D.roiID=repmat(i,numel(anaSubj),1);
                    D.mvaID=repmat(1,numel(anaSubj),1)*m;
                    D.phase=repmat(1,numel(anaSubj),1)*phase;
                    S=addstruct(S,D);%add ROI by ROI
                    D=[];
                end;
                ROI=addstruct(ROI,S); %add MVA type by type
                S=[];
                %            fileOut=[MVA_ID{m} '_Mov_ROI_spss.mat'];
                fileOut=[MVA_ID{m} '_' Phase_ID{phase} '_ROI_spss.mat'];
                save(fileOut,'meanROI');
            end;
            K=addstruct(K,ROI); %add phase by phase
            ROI=[];
            
            %            if phase == 1
            %                ROIname={'Prep L-M1','Prep R-M1','Prep L-PMd','Prep R-PMd','Prep L-SPC','Prep R-SPC'};
            %            else
            %                ROIname={'Mov L-M1','Mov R-M1','Mov L-PMd','Mov R-PMd','Mov L-SPC','Mov R-SPC'};
            %            end;
            
            
            
        end;
        %%% Figure
        %black: combined; red: temporal; blue: spatial; green: integrated
        color={[0 0 0],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
        %color={[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
        figure('Renderer', 'painters', 'Position', [10 10 400 1000])
        ROIname={'L-M1','R-M1','L-PMd','R-PMd','L-SMA','R-SMA','L-SPC','R-SPC'};
        for i=1:size(stats,2) %loop over ROIs
            %                subplot(3,(size(stats,2))/3,i)
            subplot(4,2,i)
            %                myboxplot([K.phase K.mvaID],K.meanROI,'split',K.mvaID,'subset',K.roiID==i&K.mvaID>1,'fillcolor',color,...
            %                     'leglocation','northeast') ;
            barplot([K.phase K.mvaID],K.meanROI,'split',K.mvaID,'subset',K.roiID==i&K.mvaID>0,'facecolor',color,...
                'leglocation','northeast') ;
            
            title(ROIname{i})
            drawline(0,'dir','horz');
            ylim([-1 2]);
            %axis square
        end;
        
        save('PrepProd_ROI.mat','-struct','K');
        
        ROI;
        
    case 'ROI_runComb'        %%% 2. Extract values via GUI
        %         roiIN='D:\projects\rhys\prepProd\data\imaging\ROI\roi_Elife2014.mat';
        roiIN='E:\projects\rhys\prepProd2\data\imaging\ROI\roi_Elife2014_bilateral.mat';
        
        load(roiIN);
        P={};
        %%%% Comb
        cd('E:\projects\rhys\prepProd2\data\imaging\ROI');
        MVA_ID={'Comb'};%%%TODO
        Phase_ID={'Prep','Mov'};
        ROI=[];
        K=[];
        
        subjects=[3 4 5 6 7 9 10 11 12 13 15 16 17 18 20 21 22 25 26 31 32 33 34 36 37 38 39 40 41];
        loopCounter=1;
        %         subj_name = subj_name(subjects);
        
        for phase=1:2 %prep and production
            for m=1:numel(MVA_ID) %loop for MVA type
                %                 for s=1:numel(subj_name)
                for s=subjects
                    %                P{s}=fullfile(groupDir, ['szacc_' MVA_ID{m} '_160_Mov_' subj_name{s}, '.nii']);
                    
                    P{loopCounter}=fullfile(groupDir, ['szacc_' MVA_ID{m} '_160_' Phase_ID{phase} '_' subj_name{s}, '.nii']);
                    
                    loopCounter=loopCounter+1;
                end;
                V=spm_vol(char(P));
                stats = region_getdata(V,R);
                D=[];
                S=[];
                
                for i=1:size(stats,2) %loop over ROIs
                    meanROI(:,i)=nanmean(stats{i},2); %average across voxels
                    D.meanROI=meanROI(:,i);
                    D.subj=subj';
                    D.roiID=repmat(i,numel(subj),1);
                    D.mvaID=repmat(1,numel(subj),1)*m;
                    D.phase=repmat(1,numel(subj),1)*phase;
                    S=addstruct(S,D);%add ROI by ROI
                    D=[];
                end;
                ROI=addstruct(ROI,S); %add MVA type by type
                S=[];
                %            fileOut=[MVA_ID{m} '_Mov_ROI_spss.mat'];
                fileOut=[MVA_ID{m} '_' Phase_ID{phase} '_ROI_spss.mat'];
                save(fileOut,'meanROI');
            end;
            K=addstruct(K,ROI); %add phase by phase
            ROI=[];
            
            %            if phase == 1
            %                ROIname={'Prep L-M1','Prep R-M1','Prep L-PMd','Prep R-PMd','Prep L-SPC','Prep R-SPC'};
            %            else
            %                ROIname={'Mov L-M1','Mov R-M1','Mov L-PMd','Mov R-PMd','Mov L-SPC','Mov R-SPC'};
            %            end;
            
            
            
        end;
        %%% Figure
        %black: combined; red: temporal; blue: spatial; green: integrated
        color={[0 0 0],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
        %color={[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
        figure('Renderer', 'painters', 'Position', [10 10 400 1000])
        ROIname={'L-M1','R-M1','L-PMd','R-PMd','L-SMA','R-SMA','L-SPC','R-SPC'};
        for i=1:size(stats,2) %loop over ROIs
            %                subplot(3,(size(stats,2))/3,i)
            subplot(4,2,i)
            %                myboxplot([K.phase K.mvaID],K.meanROI,'split',K.mvaID,'subset',K.roiID==i&K.mvaID>1,'fillcolor',color,...
            %                     'leglocation','northeast') ;
            barplot([K.phase K.mvaID],K.meanROI,'split',K.mvaID,'subset',K.roiID==i&K.mvaID>0,'facecolor',color,...
                'leglocation','northeast') ;
            
            title(ROIname{i})
            drawline(0,'dir','horz');
            ylim([-1 3]);
            %axis square
        end;
        
        save('PrepProd_ROI.mat','-struct','K');
        
        ROI;
        
    case 'ROI_runCon' %for contrasts only
        roiIN='D:\projects\rhys\prepProd\data\imaging\ROI\roi_Elife2014_bilateral.mat';
        
        load(roiIN);
        P={};
        %%%% Comb
        cd('D:\projects\rhys\prepProd\data\imaging\ROI');
        MVA_ID={'0002','0001','0004','0005'};%%%Prep>Rest, Mov>Rest, Prep>Mov, Mov>Prep
        
        ROI=[];
        K=[];
        
        for m=1:numel(MVA_ID) %loop for MVA type
            for s=1:numel(subj_name)
                %                P{s}=fullfile(groupDir, ['szacc_' MVA_ID{m} '_160_Mov_' subj_name{s}, '.nii']);
                P{s}=fullfile(groupDir, ['scon_' MVA_ID{m} '_' subj_name{s}, '.nii']);
            end;
            V=spm_vol(char(P));
            stats = region_getdata(V,R);
            D=[];
            S=[];
            
            for i=1:size(stats,2) %loop over ROIs
                meanROI(:,i)=nanmean(stats{i},2); %average across voxels
                D.meanROI=meanROI(:,i);
                D.subj=subj';
                D.roiID=repmat(i,numel(subj),1);
                D.mvaID=repmat(1,numel(subj),1)*m;
                S=addstruct(S,D);%add ROI by ROI
                D=[];
            end;
            ROI=addstruct(ROI,S); %add MVA type by type
            S=[];
            %            fileOut=[MVA_ID{m} '_Mov_ROI_spss.mat'];
            fileOut=[MVA_ID{m} '_ROI_spss.mat'];
            save(fileOut,'meanROI');
        end;
        
        K=ROI;
        
        %%% Figure
        %black: combined; red: temporal; blue: spatial; green: integrated
        color={[0 0 0],[1 1 1],[.2 .2 .2],[.9 .9 .9]};
        %color={[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
        figure('Renderer', 'painters', 'Position', [10 10 400 1000])
        ROIname={'L-M1','R-M1','L-PMd','R-PMd','L-SMA','R-SMA','L-SPC','R-SPC'};
        for i=1:size(stats,2) %loop over ROIs
            %                subplot(3,(size(stats,2))/3,i)
            subplot(4,2,i)
            %                myboxplot([K.phase K.mvaID],K.meanROI,'split',K.mvaID,'subset',K.roiID==i&K.mvaID>1,'fillcolor',color,...
            %                     'leglocation','northeast') ;
            barplot([K.mvaID],K.meanROI,'split',K.mvaID,'subset',K.roiID==i&K.mvaID<3,'facecolor',color,...
                'leglocation','northeast') ;
            
            title(ROIname{i})
            drawline(0,'dir','horz');
            ylim([-1 15]);
            %axis square
        end;
        
        save('Con_ROI.mat','-struct','K');
        
    case 'ROI_crossSectionPermutest'
        
        cd(crossSectionDir)
        load('PrepProd_crossSection.mat')
        
        null = zeros(24,120); %null results to compare to, acts as one-sample T test
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(spatPrep',null','true',.05,10^4,'false');
        spatPrepClust = {clusters; p_values; t_sums; permutation_distribution};
        save('spatPrepClusters','spatPrepClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(spatMov',null','true',.05,10^4,'false');
        spatMovClust = {clusters; p_values; t_sums; permutation_distribution};
        save('spatMovClusters','spatMovClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(tempPrep',null','true',.05,10^4,'false');
        tempPrepClust = {clusters; p_values; t_sums; permutation_distribution};
        save('tempPrepClusters','tempPrepClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(tempMov',null','true',.05,10^4,'false');
        tempMovClust = {clusters; p_values; t_sums; permutation_distribution};
        save('tempMovClusters','tempMovClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(intPrep',null','true',.05,10^4,'false');
        intPrepClust = {clusters; p_values; t_sums; permutation_distribution};
        save('intPrepClusters','intPrepClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(intMov',null','true',.05,10^4,'false');
        intMovClust = {clusters; p_values; t_sums; permutation_distribution};
        save('intMovClusters','intMovClust')
        
    case 'ROI_crossSectionPermutest_0.001'
        
        cd(crossSectionDir)
        load('PrepProd_crossSection.mat')
        
        null = zeros(24,120); %null results to compare to, acts as one-sample T test
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(spatPrep,null,'true',.001,10^4,'false');
        spatPrepClust = {clusters; p_values; t_sums; permutation_distribution};
        save('spatPrepClusters','spatPrepClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(spatMov,null,'true',.001,10^4,'false');
        spatMovClust = {clusters; p_values; t_sums; permutation_distribution};
        save('spatMovClusters','spatMovClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(tempPrep,null,'true',.001,10^4,'false');
        tempPrepClust = {clusters; p_values; t_sums; permutation_distribution};
        save('tempPrepClusters','tempPrepClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(tempMov,null,'true',.001,10^4,'false');
        tempMovClust = {clusters; p_values; t_sums; permutation_distribution};
        save('tempMovClusters','tempMovClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(intPrep,null,'true',.001,10^4,'false');
        intPrepClust = {clusters; p_values; t_sums; permutation_distribution};
        save('intPrepClusters','intPrepClust')
        
        [clusters, p_values, t_sums, permutation_distribution] = permutest(intMov,null,'true',.001,10^4,'false');
        intMovClust = {clusters; p_values; t_sums; permutation_distribution};
        save('intMovClusters','intMovClust')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGINNING OF SUIT ANALYSIS (CEREBELLUM)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        
    case 'suit_segment' %%%%%%%%%%%%%%%%%%% BEGINNING OF SUIT ANALYSIS (CEREBELLUM) %%%%%%%%%%%%%%%%%%%%%%%
        %Isolate the cerebellum of each participant - produces 'c_<source>_pcereb' (cerebellar mask) and '<source>_seg1/2' (grey and white matter respectively) images
        
        sn=varargin{1};
        cd([baseDir '\imaging\anatomicals\' subj_name{sn}]);
        disp(['suit_segmenting ' subj_name{sn}])
        
        anatomical = {[subj_name{sn} '_anatomical.nii']};
        
        suit_isolate_seg(anatomical, 'maskp', 0.2) %change maskp for probability value. Higher = tighter mask. Hand-correct using MRIcron if necessary
        
    case 'suit_make_mask' %restrict area of analysis to grey matter - produces 'maskbrainSUIT.nii'
        
        s=varargin{1};
        
        if isdir([baseDir '\imaging\suit\' subj_name{s}]) == 0
            mkdir([baseDir '\imaging\suit\' subj_name{s}])
        end
        
        mask=fullfile(glmDir, subj_name{s},'mask.nii');
        %         suit=fullfile(anatDir, subj_name{s},['c_', subj_name{s},'_anatomical_pcereb.nii']); %pcereb holds all cerebellum-related regions to a value of 1...
        suit=fullfile(anatDir, subj_name{s},[subj_name{s}, '_anatomical_seg1.nii']); %whereas _seg1 is only grey matter and sets extra-cerebellar regions (e.g. pons) to values other than 1...
        omask=fullfile(suitDir, subj_name{s},'maskbrainSUIT.nii');
        
        spm_imcalc_ui({mask,suit},omask,'i1>0 & i2>0.999',{}); %so including a mask of 0.999 makes sure we only include cerebellar regions.
        
    case 'MVA_searchSUIT' % Define the search lights for the MVA analysis
        
        s=varargin{1};
        
        radius=16;
        numVox=160;
        
        cd(fullfile(suitDir, subj_name{s}));
        V=spm_vol('maskbrainSUIT.nii'); %preceded by case suit_make_mask
        X=spm_read_vols(V);
        [i,j,k]=ind2sub(size(X),find(X~=0));
        vox=[i j k];
        
        [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[radius numVox],V.mat,V.dim,[],'mva160_numvoxSUIT.nii');
        save volsearch160SUIT.mat vox LI voxmin voxmax n
        
    case 'MVA_do_overallMov_suit'
        
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(suitDir, subj_name{s}, [subj_name{s}, '_accuracy_Comb_160_MovSUIT.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prod>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_do_overallPrep_suit'                 % Conduct the classification analysis 4 sequences
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(suitDir, subj_name{s}, [subj_name{s}, '_accuracy_Comb_160_Prep.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_do_spatOneout_Mov_suit'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(suitDir, subj_name{s}, [subj_name{s}, '_accuracy_Spat_160_Mov.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 0; 0 1];
            oneout=repmat(oneout,1,12);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_spatOneout_Prep_suit'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(suitDir, subj_name{s}, [subj_name{s}, '_accuracy_Spat_160_Prep.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 0; 0 1];
            oneout=repmat(oneout,1,12);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_tempOneout_Mov_suit'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(suitDir, subj_name{s}, [subj_name{s}, '_accuracy_Temp_160_Mov.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 1 0 0; 0 0 1 1];
            oneout=repmat(oneout,1,6);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_tempOneout_Prep_suit'
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(suitDir, subj_name{s}, [subj_name{s}, '_accuracy_Temp_160_Prep.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 1 0 0; 0 0 1 1];
            oneout=repmat(oneout,1,6);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_Int_Mov_suit'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(suitDir, subj_name{s}, [subj_name{s}, '_accuracy_Int_160_Mov.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prod>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'),Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_do_Int_Prep_suit'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(suitDir, subj_name{s}, [subj_name{s}, '_accuracy_Int_160_Prep.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        tstart = tic
        lmva_spm(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'),Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_zValue_suit'
        
        s=varargin{1};
        cd(fullfile(suitDir,subj_name{s}));
        
        numTests=6;
        numCat=4;
        mu=1/numCat; %mu=0.25;
        N=numTests*numCat;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        images= {'_accuracy_Comb_160_Mov','_accuracy_Comb_160_Prep','_accuracy_Int_160_Mov','_accuracy_Int_160_Prep'};
        %         images= {'_accuracy_Comb_160_Mov','_accuracy_Comb_160_Prep'};
        
        outimages={'_zacc_Comb_160_Mov','_zacc_Comb_160_Prep','_zacc_Int_160_Mov','_zacc_Int_160_Prep'};
        %         outimages={'_zacc_Comb_160_Mov','_zacc_Comb_160_Prep'};
        
        
        for j=1:numel(images)
            input_image= fullfile(suitDir,subj_name{s},[subj_name{s} images{j} '.nii']);
            output_image= fullfile(suitDir,subj_name{s},[subj_name{s} outimages{j} '.nii']);
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
        end;
        
    case 'MVA_zValue_oneOut_suit'
        
        s=varargin{1};
        cd(fullfile(suitDir,subj_name{s}));
        
        takeOneOutIter=2;
        numTests=6;
        numCat=2;
        mu=1/numCat; %mu=0.5;
        N=numTests*numCat*takeOneOutIter;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        images= {'_accuracy_Spat_160_Mov','_accuracy_Spat_160_Prep','_accuracy_Temp_160_Mov','_accuracy_Temp_160_Prep'};
        
        outimages={'_zacc_Spat_160_Mov','_zacc_Spat_160_Prep','_zacc_Temp_160_Mov','_zacc_Temp_160_Prep'};
        
        
        for j=1:numel(images)
            input_image= fullfile(suitDir,subj_name{s},[subj_name{s} images{j} '.nii']);
            output_image= fullfile(suitDir,subj_name{s},[subj_name{s} outimages{j} '.nii']);
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
        end;
        
    case 'MVA_smooth_suit'
        
        s=varargin{1};
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} '_zacc_Comb_160_Mov.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_szacc_Comb_160_Mov.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} '_zacc_Comb_160_Mov.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_szacc_Comb_160_Mov.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} '_zacc_Comb_160_Prep.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_szacc_Comb_160_Prep.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} '_zacc_Spat_160_Mov.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_szacc_Spat_160_Mov.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} '_zacc_Spat_160_Prep.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_szacc_Spat_160_Prep.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} '_zacc_Temp_160_Mov.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_szacc_Temp_160_Mov.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} '_zacc_Temp_160_Prep.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_szacc_Temp_160_Prep.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} '_zacc_Int_160_Mov.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_szacc_Int_160_Mov.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} '_zacc_Int_160_Prep.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_szacc_Int_160_Prep.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
    case 'con_smooth_suit'
        
        s = varargin{1};
        
        con=fullfile(suitDir, subj_name{s},'con_0001.nii'); %%contrast smoother
        scon=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0001.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(suitDir, subj_name{s},'con_0002.nii'); %%contrast smoother
        scon=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0002.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(suitDir, subj_name{s},'con_0003.nii'); %%contrast smoother
        scon=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0003.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(suitDir, subj_name{s},'con_0004.nii'); %%contrast smoother
        scon=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0004.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(suitDir, subj_name{s},'con_0005.nii'); %%contrast smoother
        scon=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0005.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(suitDir, subj_name{s},'con_0006.nii'); %%contrast smoother
        scon=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0006.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
    case 'suit_normalize' %normalise the isolated cerebellum to the suit atlas - produces 'affine_<source>.mat' and 'u_a_<name>.nii'
        
        sn=varargin{1};
        cd([baseDir '\imaging\anatomicals\' subj_name{sn}]);
        disp(['suit_normalizing ' subj_name{sn}])
        
        gray = {[subj_name{sn} '_anatomical_seg1.nii']}; %grey and white matter images from previous suit stage
        white = {[subj_name{sn} '_anatomical_seg2.nii']};
        isoMask = {['c_' subj_name{sn} '_anatomical_pcereb.nii']}; %isolated cerebellum
        
        job.subjND.gray = gray; %put them all into a struct...
        job.subjND.white = white;
        job.subjND.isolation = isoMask;
        
        suit_normalize_dartel(job) %run the function with the struct as the input
        
    case 'suit_normalise' %normalisation into suit space
        
        sn=varargin{1};
        cd([baseDir '\imaging\suit\']);
        
        if isdir(suitGroupDir) == 0
            mkdir(suitGroupDir)
        end
        disp(['suit_reslicing ' subj_name{sn}])
        
        inDir = [suitDir '\' subj_name{sn} '\']; %path to where data is stored (to be normalised)
        outDir = suitGroupDir;
        filenames = {'szacc_Comb_160_Prep', 'szacc_Comb_160_Mov', 'szacc_Spat_160_Prep', 'szacc_Spat_160_Mov', 'szacc_Temp_160_Prep', 'szacc_Temp_160_Mov', 'szacc_Int_160_Prep', 'szacc_Int_160_Mov'};
        
        % prepare files for input
        affine = {[anatDir '\' subj_name{sn} '\' 'Affine_' subj_name{sn} '_anatomical_seg1.mat']};
        flowfield = {[anatDir '\' subj_name{sn} '\' 'u_a_' subj_name{sn} '_anatomical_seg1.nii']};
        
        dataFiles = cell(length(filenames),1); %loop to put all full input file directories into a cell
        for i=1:length(filenames)
            dataFiles{i} = [inDir, subj_name{sn}, '_', filenames{i}, '.nii'];
        end
        
        mask = {[anatDir '\' subj_name{sn} '\' 'c_' subj_name{sn} '_anatomical_pcereb.nii']};
        
        outFiles = cell(length(filenames),1);
        for i=1:length(filenames)
            outFiles{i} = [outDir, '\', filenames{i}, '_', subj_name{sn}, '.nii'];
        end
        
        %% prepare struct for function
        job.subj.affineTr = affine; %fill job.subj. struct with respective items
        job.subj.flowfield = flowfield;
        job.subj.resample = dataFiles;
        job.subj.mask = mask;
        job.subj.outname = outFiles;
        
        %function
        suit_reslice_dartel(job)
        
    case 'suit_reslice_contrast' %reslice smoothed, normalised, individual contrast maps into SUIT space
        
        sn=varargin{1};
        cd([baseDir '\imaging\suit\' subj_name{sn}]);
        
        disp(['suit_reslicing_contrast ' subj_name{sn}])
        
        inDir = groupDir; %path to where data is stored (to be normalised)
        
        contrasts = {'scon_0001', 'scon_0002', 'scon_0003', 'scon_0004', 'scon_0005', 'scon_0006'};
        
        outDir = suitGroupDir;
        
        % prepare files for input
        affine = {[anatDir '\' subj_name{sn} '\' 'Affine_' subj_name{sn} '_anatomical_seg1.mat']};
        flowfield = {[anatDir '\' subj_name{sn} '\' 'u_a_' subj_name{sn} '_anatomical_seg1.nii']};
        
        dataFiles = cell(length(contrasts),1);
        for i=1:length(contrasts)
            dataFiles{i} = [inDir, '\', contrasts{i}, '_', subj_name{sn} '.nii'];
        end
        
        mask = {[anatDir '\' subj_name{sn} '\' 'c_' subj_name{sn} '_anatomical_pcereb.nii']};
        
        outFiles = cell(length(contrasts),1);
        for i=1:length(contrasts)
            outFiles{i} = [outDir, '\', contrasts{i}, '_', subj_name{sn}, '.nii'];
        end
        
        %%% prepare struct for function
        job.subj.affineTr = affine; %fill job.subj. struct with respective items
        job.subj.flowfield = flowfield;
        job.subj.resample = dataFiles;
        job.subj.mask = mask;
        job.subj.outname = outFiles;
        
        %function
        suit_reslice_dartel(job)
        
    case 'glm_contrastGroup_suit'
        
        %Produces second-level contrasts. Edit contrast folders, image names, and subject nifti files.
        
        dataDir = {'Mov', 'Prep', 'Error', 'PrepProd', 'ProdPrep', 'Rest'}; %%Save folders for each contrast
        images = {'scon_0001';'scon_0002';'scon_0003';'scon_0004';'scon_0005';'scon_0006'};
        
        subNii = {'_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s09.nii','_s10.nii','_s13.nii','_s16.nii','_s17.nii','_s18.nii','_s20.nii','_s21.nii','_s22.nii','_s25.nii','_s26.nii'...
            '_s31.nii','_s32.nii','_s34.nii','_s36.nii','_s38.nii','_s39.nii','_s40.nii','_s41.nii','_s42.nii'};
        contrastN = length(dataDir);
        images = repmat(images,1,length(subNii));
        subNii = repmat (subNii,length(dataDir),1);
        fileName = strcat (images,subNii);  %%Concatenate contrast files and subject names
        
        for i=1:contrastN  %%Loop across contrasts, plugging parameters into SPM.
            glmscndDir = fullfile(suitScndDir, dataDir(i));
            matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (suitGroupDir, fileName(i,:))';  %%Select files from vectors above
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            spm_jobman('run',matlabbatch);  %Run SPM
        end
        
    case 'glm_contrastEstimate_suit'
        
        dataDir = {'Mov', 'Prep', 'Error', 'PrepProd', 'ProdPrep', 'Rest'}; %%Save folders for each contrast
        contrastN = length(dataDir);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat = fullfile (suitScndDir, dataDir(i), 'SPM.mat');
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
    case 'MVA_group_suit'
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
        images = {'szacc_Comb_160_Mov';'szacc_Comb_160_Prep';'szacc_Int_160_Mov';'szacc_Int_160_Prep';'szacc_Spat_160_Mov';'szacc_Spat_160_Prep';'szacc_Temp_160_Mov';'szacc_Temp_160_Prep'};
        
        %         dataDir = {'MVA_comb_mov', 'MVA_comb_prep'}; %%just overall decoding
        %         images = {'szacc_Comb_160_Mov';'szacc_Comb_160_Prep'};
        
        %ONLY PARTICIPANTS WHO MODULATED TIMING
        subNii = {'_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s09.nii','_s10.nii','_s13.nii','_s16.nii','_s17.nii','_s18.nii','_s20.nii','_s21.nii','_s22.nii','_s25.nii','_s26.nii'...
            '_s31.nii','_s32.nii','_s34.nii','_s36.nii','_s38.nii','_s39.nii','_s40.nii','_s41.nii','_s42.nii'}; %3 5 6 7 9 10 13 16 17 18 20 21 22 25 26 31 32 34 36 38 39 40 41 42
        
        %ALL PARTICIPANTS regardless of timing modulation ***TO DO***
        %         subNii = {'s03.nii','s05.nii','s06.nii','s07.nii','s09.nii','s10.nii','s13.nii','s16.nii','s17.nii','s18.nii','s20.nii','s21.nii','s22.nii','s25.nii','s26.nii'...
        %             's31.nii','s32.nii','s34.nii','s36.nii','s38.nii','s39.nii','s40.nii','s41.nii'}; %3 4 5 6 7 9 10 11 12 13 15 16 17 18 20 21 22 25 26 31 32 33 34 36 37 38 39 40 41 42
        
        contrastN = length(dataDir);
        images = repmat(images,1,length(subNii));
        subNii = repmat (subNii,length(dataDir),1);
        fileName = strcat (images,subNii);  %%Concatenate contrast files and subject names
        
        for i=1:contrastN  %%Loop across contrasts
            glmscndDir = fullfile(suitGroupDir, dataDir(i));
            matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;  %Adjust directory
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (suitGroupDir, fileName(i,:))';  %%Select files from matrix
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
        %AVERAGE GROUP DATA
        %open spm fmri and select 'imcalc'
        %choose all images from participants of interest
        %(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+i18+i19+i20+i21+i22+i23+i24)/24
        %(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+i18+i19+i20+i21+i22+i23+i24)/24;
        
    case 'MVA_estimate_suit'
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
        %         dataDir = {'MVA_comb_mov', 'MVA_comb_prep'}; %%Save folders for each contrast
        contrastN = length(dataDir);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat = fullfile (suitGroupDir, dataDir(i), 'SPM.mat');  %Adjust directory
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Whole brain image reslicing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'suit_reslice_whole'
        
        sn=varargin{1};
        cd([baseDir '\imaging\suit\']);
        
        %         if isdir([baseDir '\imaging\suit\' subj_name{sn}]) == 0
        %             mkdir([baseDir '\imaging\suit\' subj_name{sn}])
        %         end
        disp(['suit_reslicing ' subj_name{sn}])
        
        inDir = [glmDir '\' subj_name{sn} '\']; %path to where data is stored (to be normalised)
        outDir = [baseDir '\imaging\suit\'];
        filenames = {'_szacc_Comb_160_Prep.nii', '_szacc_Comb_160_Mov.nii', '_szacc_Spat_160_Prep.nii', '_szacc_Spat_160_Mov.nii', '_szacc_Temp_160_Prep.nii', '_szacc_Temp_160_Mov.nii', '_szacc_Int_160_Prep.nii', '_szacc_Int_160_Mov.nii'};
        filenamesT = {'spmT_0001.nii', 'spmT_0002.nii', 'spmT_0003.nii', 'spmT_0004.nii', 'spmT_0005.nii', 'spmT_0006.nii'};
        
        %% prepare files for input
        affine = {[anatDir '\' subj_name{sn} '\' 'Affine_' subj_name{sn} '_anatomical_seg1.mat']};
        flowfield = {[anatDir '\' subj_name{sn} '\' 'u_a_' subj_name{sn} '_anatomical_seg1.nii']};
        
        dataFiles = cell(length(filenames),1); %separate loops to fill filenames for classifiers and T contrast maps
        for i=1:length(filenames)
            dataFiles{i} = [inDir, subj_name{sn}, filenames{i}];
        end
        dataFilesT = cell(length(filenamesT),1);
        for i=1:length(filenamesT)
            dataFilesT{i} = [inDir, filenamesT{i}];
        end
        dataFiles = [dataFiles; dataFilesT];
        
        mask = {[anatDir '\' subj_name{sn} '\' 'c_' subj_name{sn} '_anatomical_pcereb.nii']};
        
        outFiles = cell(length(filenames),1);
        for i=1:length(filenames)
            outFiles{i} = [outDir, '\', subj_name{sn}, filenames{i}];
        end
        outFilesT = cell(length(filenamesT),1);
        for i=1:length(filenamesT)
            outFilesT{i} = [outDir, '\', subj_name{sn}, '_' filenamesT{i}];
        end
        outFiles = [outFiles; outFilesT];
        
        %% prepare struct for function
        job.subj.affineTr = affine; %fill job.subj. struct with respective items
        job.subj.flowfield = flowfield;
        job.subj.resample = dataFiles;
        job.subj.mask = mask;
        job.subj.outname = outFiles;
        
        %function
        suit_reslice_dartel(job)
        
    case 'suit_group_whole'
        
        %% Produces second-level contrasts. Edit contrast folders, image names, and subject nifti files.
        
        dataDir = {'Mov', 'Prep', 'Error', 'PrepProd', 'ProdPrep', 'Rest'}; %%Save folders for each contrast
        images = {'_spmT_0001.nii';'_spmT_0002.nii';'_spmT_0003.nii';'_spmT_0004.nii';'_spmT_0005.nii';'_spmT_0006.nii'};
        
        %         dataDir = {'Mov','Prep','Points'}; %%Save folders for each contrast
        %         images = {'scon_0001';'scon_0002';'scon_0006'};
        
        %         dataDir = {'Prep'}; %%Save folders for each contrast
        %         images = {'scon_0002'};
        
        subNii = {'s03','s04','s05','s06','s07','s09','s10','s11'};
        contrastN = length(dataDir);
        images = repmat(images,1,length(subNii));
        subNii = repmat (subNii,length(dataDir),1);
        fileName = strcat (subNii, images);  %%Concatenate contrast files and subject names
        
        for i=1:contrastN  %%Loop across contrasts, plugging parameters into SPM.
            glmscndDir = fullfile(suitDir, 'data', dataDir(i));
            matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (suitDir, fileName(i,:))';  %%Select files from vectors above
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            spm_jobman('run',matlabbatch);  %Run SPM
        end
        
        %% Produces second level MVPA accuracy maps. Edit folders, image names, and subject nifti files.
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
        images = {'_szacc_Comb_160_Mov.nii';'_szacc_Comb_160_Prep.nii';'_szacc_Int_160_Mov.nii';'_szacc_Int_160_Prep.nii';'_szacc_Spat_160_Mov.nii';'_szacc_Spat_160_Prep.nii';'_szacc_Temp_160_Mov.nii';'_szacc_Temp_160_Prep.nii'};
        %         subNii = {'_s01.nii','_s02.nii','_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s08.nii','_s09.nii','_s10.nii'};
        subNii = {'s03','s04','s05','s06','s07','s09','s10','s11'};
        contrastN = length(dataDir);
        images = repmat(images,1,length(subNii));
        subNii = repmat (subNii,length(dataDir),1);
        fileName = strcat (subNii, images);  %%Concatenate contrast files and subject names
        
        for i=1:contrastN  %%Loop across contrasts
            glmscndDir = fullfile(suitDir, dataDir(i));
            matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;  %Adjust directory
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (suitDir, fileName(i,:))';  %%Select files from matrix
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
    case 'suit_estimate_whole'
        
        %%Contrast
        dataDir = {'Mov', 'Prep', 'Error', 'PrepProd', 'ProdPrep', 'Rest'}; %%Save folders for each contrast
        contrastN = length(dataDir);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat = fullfile (suitDir, 'data', dataDir(i), 'SPM.mat');
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
        %% MVPA
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
        contrastN = length(dataDir);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat = fullfile (suitDir, dataDir(i), 'SPM.mat');  %Adjust directory
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
    case 'suit_surface_map_whole'
        
        %MANUALLY OPEN RESULTS IN SPM FIRST
        % % SPM -> Results -> suit -> second level folders (contrasts under 'data', MVPA in individual folders.
        dataFolder = {'data\Mov','data\Prep','data\PrepProd','data\ProdPrep','data\Error','data\Rest',... %folders containing data to be projected
            'MVA_comb_mov','MVA_comb_prep','MVA_spat_mov','MVA_spat_prep','MVA_temp_mov','MVA_temp_prep','MVA_int_mov','MVA_int_prep'};
        
        for i =1:length(dataFolder)
            cd([suitDir, '\', dataFolder{i}]) %go into folder location for data and to save into.
            
            flatmapVector = suit_map2surf('spmT_0001.nii'); %produces vector information regarding flatmap
            flatmapGifti = gifti(flatmapVector); %convert to gifti format (compatibility with connectome workbench)
            save(flatmapGifti,'spmT_0001.func.gii')
            save('spmT_0001_flatmap.mat','flatmapVector')
        end
        
        %% for within-matlab visualisation of flatmap, load respective 'spmT_0001_flatmap.mat'
        % then use suit_plotflatmap(flatmapVector, 'cmap', hot, 'cscale', [0 7.00], 'threshold', 3.48) replacing numbers with desired values.
        
    case 'suit_roi_whole'
        
        dataFolder = {'data\Mov','data\Prep','data\PrepProd','data\ProdPrep','data\Error','data\Rest',... %folders containing data to be projected
            'MVA_comb_mov','MVA_comb_prep','MVA_spat_mov','MVA_spat_prep','MVA_temp_mov','MVA_temp_prep','MVA_int_mov','MVA_int_prep'};
        
        dataNames = {'Mov','Prep','PrepProd','ProdPrep','Error','Rest',... %names for saving
            'MVA_comb_mov','MVA_comb_prep','MVA_spat_mov','MVA_spat_prep','MVA_temp_mov','MVA_temp_prep','MVA_int_mov','MVA_int_prep'};
        
        for i =1:length(dataFolder)
            
            cd([suitDir, '\', dataFolder{i}]) %go into folder location for data.
            saveDir = [suitDir, '\', 'RoI', '\', dataNames{i}, '.txt'];
            filename = {'spmT_0001.nii'};
            
            suit_ROI_summarize(filename, 'outfilename', saveDir) %For definitions of RoIs using standard atlas visit http://www.diedrichsenlab.org/imaging/mdtb.htm
            
        end
        
        
        
        
        
        
        
    case 'subcortical_make_nii' %%%%%%%%%%%%%%%%%%% BEGINNING OF SUBCORTICAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%
        
        sn = varargin{1}; %requires recon-all to be completed through freesurfer for each participant, and for 'aseg' file to be moved to relevant directory
        
        %%Anatomical:
        try
            source = fullfile(subcorticalAreaDir, subj_name{sn}, [subj_name{sn} '_' 'aseg.nii.gz']);
            dest = fullfile(subcorticalAreaDir, subj_name{sn});
            
            gunzip(source,dest); %unzip
            
            disp('segmentation unzip done')
        catch
            disp('No segmentation files. Most likely a typo in file name or missing destination directory!');
        end;
        
    case 'subcortical_make_structs' %extracts each subcortical structure from aseg file, and creates respective nii files
        
        sn = varargin{1};
        
        subcortValues = [10, 11, 12, 13, 16, 17, 18, 26, 28, 30, 31, 49, 50, 51, 52, 53, 54, 58, 60, 87, 63];
        
        for i=1:length(subcortStructs)
            
            cd([subcorticalAreaDir, '\', subj_name{sn}])
            
            matlabbatch{1}.spm.util.imcalc.input = {[subcorticalAreaDir, '\' subj_name{sn} '\' subj_name{sn}, '_aseg.nii,1']};
            matlabbatch{1}.spm.util.imcalc.output = [subj_name{sn}, '_', subcortStructs{i}];
            matlabbatch{1}.spm.util.imcalc.expression = ['i1 == ', num2str(subcortValues(i))];
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            
            spm_jobman('run',matlabbatch);
        end
        
    case 'subcortical_reslice_structs' %reslices subcortical niftis into functional scan resolution
        
        sn = varargin{1};
        
        refImageDir = [glmDir, '\', subj_name{sn}]; %directories for functional reference image (beta 1)
        subcortImageDir = [subcorticalAreaDir, '\', subj_name{sn}]; %and for subcortical images
        
        for r=1:length(subcortStructs) %loop through subcortical regions and reslice them to beta image resolution
            
            matlabbatch{1}.spm.spatial.realign.write.data = {[refImageDir, '\beta_0001.nii']; [subcortImageDir, '\', subj_name{sn}, '_', subcortStructs{r}, '.nii']};
            matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
            
            spm_jobman('run',matlabbatch);
            disp([subcortStructs{r} ' resliced'])
        end
        
    case 'subcortical_voxel_counts' %provides number of voxels for each subcortical structure
        
        sn = varargin{1};
        voxCount = nan(1,length(subcortStructs));
        %         for i =anaSubj(1:end ~= 2) %for each subject...
        for j=1:length(subcortStructs) %and each subcortical structure...
            
            %cd into relevant subject directory
            cd([subcorticalDir, '\',  subj_name{sn}])
            
            %read in .nii file for each subcortical structure
            vol = spm_vol([subj_name{sn}, '_', subcortStructs{j}, '.nii']);
            readVol = spm_read_vols(vol);
            
            %count the number of non-zero elements
            voxelCount = nnz(readVol);
            
            %add the count to matrix
            voxCount(:,j) = voxelCount;
            
        end
        %         end
        
        %add names of structures to top row of matrix and save to excel
        voxCount = num2cell(voxCount);
        voxCountSave = [subcortStructs; voxCount];
        
        cd(['E:\projects\rhys\prepProd2\data\imaging\subcortical', '\' subj_name{sn}])
        xlswrite([subj_name{sn} '_voxelCounts'], voxCountSave)
        
    case 'subcortical_voxel_counts_group' %like case above, but produces one big group excel sheet. Quite slow...
        
        %set voxCount to size subj x structure for later loop
        groupVoxCount = nan(length(anaSubj),length(subcortStructs));
        %loop counter for subject rows in excel sheet
        loopCount = 1;
        
        for i =anaSubj(1:end ~= 2) %for each subject...
            for j=1:length(subcortStructs) %and each subcortical structure...
                
                %cd into relevant subject directory
                cd([subcorticalDir, '\',  subj_name{i}])
                
                %read in .nii file for each subcortical structure
                vol = spm_vol(['r' subj_name{i}, '_', subcortStructs{j}, '.nii']);
                readVol = spm_read_vols(vol);
                
                %count the number of non-zero elements
                voxelCount = nnz(readVol);
                
                %add the count to matrix
                groupVoxCount(loopCount,j) = voxelCount;
                
                
                
            end
            %loop counter increase
            loopCount = loopCount + 1;
        end
        
        %add names of structures to top row of matrix and save to excel
        groupVoxCount = num2cell(groupVoxCount);
        voxCountSave = [subcortStructs; groupVoxCount];
        
        cd(subcorticalDir')
        xlswrite('group_voxelCounts', voxCountSave)
        
        
        
        
    case 'make_mask_subcorticalSearch'  %SUBCORTICAL SEARCHLIGHT ANALYSIS - Makes restricted analysis mask for MVA
        
        s=varargin{1};
        
        for i=1:length(subcortStructs);
            funMask=fullfile(glmDir, subj_name{s},'mask.nii');
            omask=fullfile(subcorticalSearchDir, subj_name{s},['mr' subj_name{s} '_' subcortStructs{i} '.nii']); %output mask to be used in the future
            subcort = fullfile(subcorticalSearchDir, subj_name{s},['r' subj_name{s} '_' subcortStructs{i} '.nii']);
            
            spm_imcalc({funMask,subcort},omask,'i1 >0.01 & i2 > 0.1',{}); %recorded activity in brain (grey + white matter)
        end
        
    case 'MVA_search_subcorticalSearch' % Define the search lights for the MVA analysis
        
        s=varargin{1};
        
        for r = 1:length(subcortStructs);
            
            radius=16;
            numVox=160;
            cd(fullfile(subcorticalSearchDir, subj_name{s}));
            V=spm_vol(['mr' subj_name{s} '_' subcortStructs{r} '.nii']); %if preceded by case MVA_mask
            X=spm_read_vols(V);
            [i,j,k]=ind2sub(size(X),find(X~=0));
            vox=[i j k];
            [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[radius numVox],V.mat,V.dim,[],'mva160_numvox.nii');
            save ([subj_name{s} '_' subcortStructs{r} '_' 'volsearch160.mat'], 'vox', 'LI', 'voxmin', 'voxmax', 'n')
        end
        
    case 'MVA_do_overallMov_subcorticalSearch'                 % Conduct the classification analysis 4 sequences
        
        s=varargin{1};
        
        for r=1:length(subcortStructs);
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalSearchDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Comb_160_Mov.nii'])};
            
            
            % Generate column indices for Cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
                train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_volsearch160.mat']),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart)
        end
        
    case 'MVA_do_overallPrep_subcorticalSearch'                 % Conduct the classification analysis 4 sequences
        
        s=varargin{1};
        
        for r=1:length(subcortStructs);
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalSearchDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Comb_160_Prep.nii'])};
            
            
            % Generate column indices for Cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
                train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_volsearch160.mat']),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart)
        end
        
    case 'MVA_do_spatOneout_Mov_subcorticalSearch'
        
        s=varargin{1};
        
        for r=1:length(subcortStructs);
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalSearchDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Spat_160_Mov.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 0; 0 1];
            oneout=repmat(oneout,1,12);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_volsearch160.mat']),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_spatOneout_Prep_subcorticalSearch'
        
        s=varargin{1};
        
        for r=1:length(subcortStructs);
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalSearchDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Spat_160_Prep.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 0; 0 1];
            oneout=repmat(oneout,1,12);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_volsearch160.mat']),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_tempOneout_Mov_subcorticalSearch'
        
        s=varargin{1};
        
        for r=1:length(subcortStructs);
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalSearchDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Temp_160_Mov.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 1 0 0; 0 0 1 1];
            oneout=repmat(oneout,1,6);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_volsearch160.mat']),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_tempOneout_Prep_subcorticalSearch'
        
        s=varargin{1};
        
        for r=1:length(subcortStructs);
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalSearchDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Temp_160_Prep.nii'])};
            
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            %%%
            oneout=[1 1 0 0; 0 0 1 1];
            oneout=repmat(oneout,1,6);
            j=0;
            for i=1:2:nrruns*2
                j=j+1;
                test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
                test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
                
                
                train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
                train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
                
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_volsearch160.mat']),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart);
            
        end;
        
    case 'MVA_do_Int_Mov_subcorticalSearch'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        
        for r=1:length(subcortStructs);
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalSearchDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Int_160_Mov.nii'])};
            
            
            % Generate column indices for Cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
                train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_volsearch160.mat'])...
                ,Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
            telapsed = toc(tstart)
        end
        
    case 'MVA_do_Int_Prep_subcorticalSearch'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        
        for r=1:length(subcortStructs)
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalSearchDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Int_160_Prep.nii'])};
            
            
            % Generate column indices for Cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
                train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_volsearch160.mat'])...
                ,Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
            telapsed = toc(tstart)
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'MVA_zValue_subcorticalSearch'
        
        s=varargin{1};
        cd(fullfile(subcorticalSearchDir,subj_name{s}));
        
        for r=1:length(subcortStructs)
            
            numTests=6;
            numCat=4;
            mu=1/numCat; %mu=0.25;
            N=numTests*numCat;
            sigma=sqrt(mu*(1-mu)*1/N);
            
            images= {'_accuracy_Comb_160_Mov','_accuracy_Comb_160_Prep','_accuracy_Int_160_Mov','_accuracy_Int_160_Prep'};
            %         images= {'_accuracy_Comb_160_Mov','_accuracy_Comb_160_Prep'};
            
            outimages={'_zacc_Comb_160_Mov','_zacc_Comb_160_Prep','_zacc_Int_160_Mov','_zacc_Int_160_Prep'};
            %         outimages={'_zacc_Comb_160_Mov','_zacc_Comb_160_Prep'};
            
            
            for j=1:numel(images)
                input_image= fullfile([subj_name{s} '_' subcortStructs{r} images{j} '.nii']);
                output_image= fullfile([subj_name{s} '_' subcortStructs{r} outimages{j} '.nii']);
                spmj_imcalc_mtx(input_image, output_image,...
                    sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
            end
        end
        
    case 'MVA_zValue_oneOut_subcorticalSearch'
        
        s=varargin{1};
        cd(fullfile(subcorticalSearchDir,subj_name{s}));
        
        for r=1:length(subcortStructs)
            
            takeOneOutIter=2;
            numTests=6;
            numCat=2;
            mu=1/numCat; %mu=0.5;
            N=numTests*numCat*takeOneOutIter;
            sigma=sqrt(mu*(1-mu)*1/N);
            
            images= {'_accuracy_Spat_160_Mov','_accuracy_Spat_160_Prep','_accuracy_Temp_160_Mov','_accuracy_Temp_160_Prep'};
            
            outimages={'_zacc_Spat_160_Mov','_zacc_Spat_160_Prep','_zacc_Temp_160_Mov','_zacc_Temp_160_Prep'};
            
            
            for j=1:numel(images)
                input_image= fullfile([subj_name{s} '_' subcortStructs{r} images{j} '.nii']);
                output_image= fullfile([subj_name{s} '_' subcortStructs{r} outimages{j} '.nii']);
                spmj_imcalc_mtx(input_image, output_image,...
                    sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
            end
        end
        
    case 'MVA_smooth_subcorticalSearch' %%%Smoothing in subject space
        
        for r=1:length(subcortStructs)
            
            s=varargin{1};
            comb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Comb_160_Mov.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Comb_160_Mov.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Comb_160_Prep.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Comb_160_Prep.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Spat_160_Mov.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Spat_160_Mov.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Spat_160_Prep.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Spat_160_Prep.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Temp_160_Mov.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Temp_160_Mov.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Temp_160_Prep.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Temp_160_Prep.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Int_160_Mov.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Int_160_Mov.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Int_160_Prep.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalSearchDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Int_160_Prep.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            %smooth other images here as required
            
        end
        
    case 'MNI_normalisation_subcorticalSearch'
        
        s=varargin{1};
        mkdir(fullfile(subcorticalSearchGroupDir,'data')); % folder for each contrast
        
        %MVPA accuracy maps
        images= {'szacc_Comb_160_Mov.nii','szacc_Comb_160_Prep.nii',...
            'szacc_Int_160_Mov.nii','szacc_Int_160_Prep.nii',...
            'szacc_Spat_160_Mov.nii','szacc_Spat_160_Prep.nii',...
            'szacc_Temp_160_Mov.nii','szacc_Temp_160_Prep.nii'}; % please add other images as required, e.g. spmT_...
        
        for r=1:length(subcortStructs)
            for i=1:length(images)
                anaImages{1} = images{i};
                defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
                for j=1:numel(anaImages)
                    [~,name,ext]=spm_fileparts(anaImages{j});
                    sn_images{j}= fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_' anaImages{j}]);
                    out_images{j}= fullfile(subcorticalSearchGroupDir,[name '_' subj_name{s} '_' subcortStructs{r} '.nii']);
                    spmj_normalization_write(defor, sn_images,'outimages',out_images); %Trilinear interpolation
                end
            end
        end
        
    case 'MVA_group_subcorticalSearch'
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
        images = {'szacc_Comb_160_Mov';'szacc_Comb_160_Prep';'szacc_Int_160_Mov';'szacc_Int_160_Prep';'szacc_Spat_160_Mov';'szacc_Spat_160_Prep';'szacc_Temp_160_Mov';'szacc_Temp_160_Prep'};
        
        %         dataDir = {'MVA_comb_mov', 'MVA_comb_prep'}; %%Save folders for each contrast
        %         images = {'szacc_Comb_160_Mov';'szacc_Comb_160_Prep'};
        
        %ONLY PARTICIPANTS WHO MODULATED TIMING
        subNii = {'_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s09.nii','_s10.nii','_s13.nii','_s16.nii','_s17.nii','_s18.nii','_s20.nii','_s21.nii','_s22.nii','_s25.nii','_s26.nii'...
            '_s31.nii','_s32.nii','_s34.nii','_s36.nii','_s38.nii','_s39.nii','_s40.nii','_s41.nii','_s42.nii'}; %3 5 6 7 9 10 13 16 17 18 20 21 22 25 26 31 32 34 36 38 39 40 41 42
        
        %ALL PARTICIPANTS regardless of timing modulation ***TO DO***
        %         subNii = {'_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s09.nii','_s10.nii','_s13.nii','_s16.nii','_s17.nii','_s18.nii','_s20.nii','_s21.nii','_s22.nii','_s25.nii','_s26.nii'...
        %             '_s31.nii','_s32.nii','_s34.nii','_s36.nii','_s38.nii','_s39.nii','_s40.nii','_s41.nii'}; %3 4 5 6 7 9 10 11 12 13 15 16 17 18 20 21 22 25 26 31 32 33 34 36 37 38 39 40 41 42
        
        contrastN = length(dataDir);
        images = repmat(images,1,length(subNii));
        subNii = repmat (subNii,length(dataDir),1);
        fileName = strcat (images,subNii);  %%Concatenate contrast files and subject names
        
        for i=1:contrastN  %%Loop across contrasts
            glmscndDir = fullfile(subcorticalSearchGroupDir, dataDir(i));
            matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;  %Adjust directory
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (subcorticalSearchGroupDir, fileName(i,:))';  %%Select files from matrix
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
        %AVERAGE GROUP DATA
        %open spm fmri and select 'imcalc'
        %choose all images from participants of interest
        %(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+16+i17+i18+i19+i20+i21+i22+i23+i24)/24
        
        
        
        
        
        
        
        
        
    case 'subcortical_define_mva_area' %SUBCORTICAL AREA ANALYSIS - defines area in each subcortical region to be passed to later MVA cases
        
        s=varargin{1};
        
        % load respective subcortical .nii and extract x y z coordinates for voxels
        for r=1:length(subcortStructs)
            cd(fullfile(subcorticalAreaDir, subj_name{s}));
            V=spm_vol(['r' subj_name{s} '_' subcortStructs{r} '.nii']);
            X=spm_read_vols(V);
            [i,j,k]=ind2sub(size(X),find(X~=0));
            vox=[i j k];
            
            LI{1} = sub2ind(V.dim,vox(:,1),vox(:,2),vox(:,3)); %convert x y z coordinates into linear index for voxels
            voxmin = min(vox,[],1); %min x y z coordinates for each voxel
            voxmax = max(vox,[],1); %max as above
            n = size(vox,1); %number of voxels
            
            save(sprintf('%s_decodeArea_%s.mat',subj_name{s},subcortStructs{r}), 'vox', 'LI', 'voxmin', 'voxmax', 'n')
            disp([subcortStructs{r} ' done'])
        end
        
    case 'MVA_do_overallMov_subcorticalArea' %runs combMov LDA on all voxels within each subcortical region
        
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s})); %load SPM data
        load SPM;
        nrruns=length(SPM.nscan);
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [~,col] = find(prod>0); %identify relevant beta files and assign to pselect
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end
        
        predcat = zeros(length(subcortStructs),1);
        
        %%% loop across decoding areas for respective subcortical regions
        for r=1:length(subcortStructs)
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{r} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            
            pred = combinedclass(X, c, run, train, test); %classifier function which gives us accuracy value as a decimal
            % chance value = 0.25
            %             disp(['Prediction value ' num2str(pred) ' in ' subcortStructs{r}]) %display all in window
            
            predcat(r) = pred; %save prediction value into variable
        end
        
        combMov = [subcortStructs', num2cell(predcat)]; %save participant's decoding values to .mat file
        cd([subcorticalAreaDir, '\', subj_name{s}])
        save([subj_name{s}, '_combMov_subcortical'],'combMov')
        
    case 'MVA_do_overallPrep_subcorticalArea' %runs combPrep LDA on all voxels within each subcortical region
        
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        predcat = zeros(length(subcortStructs),1);
        
        %%% loop across decoding areas for respective subcortical regions
        for r=1:length(subcortStructs)
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{r} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            
            pred = combinedclass(X, c, run, train, test); %classifier function which gives us accuracy value as a decimal
            % chance value = 0.25
            %             disp(['Prediction value ' num2str(pred) ' in ' subcortStructs{r}]) %display all in window
            
            %%%TO DO%%%
            predcat(r) = pred; %save prediction value into variable
        end
        
        combPrep = [subcortStructs', num2cell(predcat)]; %save participant's decoding values to .mat file
        cd([subcorticalAreaDir, '\', subj_name{s}])
        save([subj_name{s}, '_combPrep_subcortical'],'combPrep')
        
    case 'MVA_do_spatOneout_Mov_subcorticalArea' %runs spatMov LDA on all voxels within each subcortical region
        
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        
        % Generate column indices for cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        %%%
        oneout=[1 0; 0 1];
        oneout=repmat(oneout,1,12);
        j=0;
        for i=1:2:nrruns*2
            j=j+1;
            test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
            test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
            
            
            train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
            train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
            
        end;
        
        [row,col] = find(prod>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        predcat = zeros(length(subcortStructs),1);
        
        %%% loop across decoding areas for respective subcortical regions
        for r=1:length(subcortStructs)
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{r} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            
            pred = combinedclass(X, c, run, train, test); %classifier function which gives us accuracy value as a decimal
            % chance value = 0.25
            %             disp(['Prediction value ' num2str(pred) ' in ' subcortStructs{r}]) %display all in window
            
            %%%TO DO%%%
            predcat(r) = pred; %save prediction value into variable
        end
        
        spatMov = [subcortStructs', num2cell(predcat)]; %save participant's decoding values to .mat file
        cd([subcorticalAreaDir, '\', subj_name{s}])
        save([subj_name{s}, '_spatMov_subcortical'],'spatMov')
        
    case 'MVA_do_spatOneout_Prep_subcorticalArea' %runs spatPrep LDA on all voxels within each subcortical region
        
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Spat_160_Prep.nii'])};
        
        
        % Generate column indices for cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        %%%
        oneout=[1 0; 0 1];
        oneout=repmat(oneout,1,12);
        j=0;
        for i=1:2:nrruns*2
            j=j+1;
            test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
            test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
            
            
            train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
            train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
            
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        predcat = zeros(length(subcortStructs),1);
        
        %%% loop across decoding areas for respective subcortical regions
        for r=1:length(subcortStructs)
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{r} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            
            pred = combinedclass(X, c, run, train, test); %classifier function which gives us accuracy value as a decimal
            % chance value = 0.25
            %             disp(['Prediction value ' num2str(pred) ' in ' subcortStructs{r}]) %display all in window
            
            %%%TO DO%%%
            predcat(r) = pred; %save prediction value into variable
        end
        
        spatPrep = [subcortStructs', num2cell(predcat)]; %save participant's decoding values to .mat file
        cd([subcorticalAreaDir, '\', subj_name{s}])
        save([subj_name{s}, '_spatPrep_subcortical'],'spatPrep')
        
    case 'MVA_do_tempOneout_Mov_subcorticalArea' %runs tempMov LDA on all voxels within each subcortical region
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        
        % Generate column indices for cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        %%%
        oneout=[1 1 0 0; 0 0 1 1];
        oneout=repmat(oneout,1,6);
        j=0;
        for i=1:2:nrruns*2
            j=j+1;
            test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
            test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
            
            
            train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
            train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
            
        end;
        
        [row,col] = find(prod>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        predcat = zeros(length(subcortStructs),1);
        
        %%% loop across decoding areas for respective subcortical regions
        for r=1:length(subcortStructs)
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{r} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            
            pred = combinedclass(X, c, run, train, test); %classifier function which gives us accuracy value as a decimal
            % chance value = 0.25
            %             disp(['Prediction value ' num2str(pred) ' in ' subcortStructs{r}]) %display all in window
            
            %%%TO DO%%%
            predcat(r) = pred; %save prediction value into variable
        end
        
        tempMov = [subcortStructs', num2cell(predcat)]; %save participant's decoding values to .mat file
        cd([subcorticalAreaDir, '\', subj_name{s}])
        save([subj_name{s}, '_tempMov_subcortical'],'tempMov')
        
    case 'MVA_do_tempOneout_Prep_subcorticalArea' %runs tempPrep LDA on all voxels within each subcortical region
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        
        % Generate column indices for cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        %%%
        oneout=[1 1 0 0; 0 0 1 1];
        oneout=repmat(oneout,1,6);
        j=0;
        for i=1:2:nrruns*2
            j=j+1;
            test{i}   =find(run==j & oneout(1,:)==1); % Classify S1 vs S2 with T1
            test{i+1} =find(run==j & oneout(2,:)==1); % Classify S1 vs S2 with T2
            
            
            train{i}  =find(run~=j & oneout(1,:)~=1); % Train on S1 vs S2 with T2 ....
            train{i+1}=find(run~=j & oneout(2,:)~=1); % Train on S1 vs S2 with T1 in different runs from testing
            
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        predcat = zeros(length(subcortStructs),1);
        
        %%% loop across decoding areas for respective subcortical regions
        for r=1:length(subcortStructs)
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{r} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            
            pred = combinedclass(X, c, run, train, test); %classifier function which gives us accuracy value as a decimal
            % chance value = 0.25
            %             disp(['Prediction value ' num2str(pred) ' in ' subcortStructs{r}]) %display all in window
            
            %%%TO DO%%%
            predcat(r) = pred; %save prediction value into variable
        end
        
        tempPrep = [subcortStructs', num2cell(predcat)]; %save participant's decoding values to .mat file
        cd([subcorticalAreaDir, '\', subj_name{s}])
        save([subj_name{s}, '_tempPrep_subcortical'],'tempPrep')
        
    case 'MVA_do_Int_Mov_subcorticalArea' %runs intMov LDA on all voxels within each subcortical region
        
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s})); %load SPM data
        load SPM;
        nrruns=length(SPM.nscan);
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [~,col] = find(prod>0); %identify relevant beta files and assign to pselect
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end
        
        predcat = zeros(length(subcortStructs),1);
        
        %%% loop across decoding areas for respective subcortical regions
        for r=1:length(subcortStructs)
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{r} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            
            pred = prepProd2_combinedclass_corrected4Main(X, c, run, train, test); %classifier function which gives us accuracy value as a decimal
            % chance value = 0.25
            %             disp(['Prediction value ' num2str(pred) ' in ' subcortStructs{r}]) %display all in window
            
            predcat(r) = pred; %save prediction value into variable
        end
        
        intMov = [subcortStructs', num2cell(predcat)]; %save participant's decoding values to .mat file
        cd([subcorticalAreaDir, '\', subj_name{s}])
        save([subj_name{s}, '_intMov_subcortical'],'intMov')
        
    case 'MVA_do_Int_Prep_subcorticalArea' %runs intPrep LDA on all voxels within each subcortical region
        
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end;
        
        [row,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end;
        
        predcat = zeros(length(subcortStructs),1);
        
        %%% loop across decoding areas for respective subcortical regions
        for r=1:length(subcortStructs)
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{r} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            
            pred = prepProd2_combinedclass_corrected4Main(X, c, run, train, test); %classifier function which gives us accuracy value as a decimal
            % chance value = 0.25
            %             disp(['Prediction value ' num2str(pred) ' in ' subcortStructs{r}]) %display all in window
            
            %%%TO DO%%%
            predcat(r) = pred; %save prediction value into variable
        end
        
        intPrep = [subcortStructs', num2cell(predcat)]; %save participant's decoding values to .mat file
        cd([subcorticalAreaDir, '\', subj_name{s}])
        save([subj_name{s}, '_intPrep_subcortical'],'intPrep')
        
    case 'MVA_group_subcorticalArea'
        
        classifiers = {'combMov', 'combPrep', 'spatMov', 'spatPrep', 'tempMov', 'tempPrep', 'intMov', 'intPrep'};
        combMovGroup = subcortStructs'; combPrepGroup = subcortStructs'; spatMovGroup = subcortStructs'; spatPrepGroup = subcortStructs';
        tempMovGroup = subcortStructs'; tempPrepGroup = subcortStructs'; intMovGroup = subcortStructs'; intPrepGroup = subcortStructs';
        
        for s = anaSubj(1:end ~= 2)
            cd(fullfile(subcorticalAreaDir, subj_name{s}))
            load([subj_name{s} '_combMov_subcortical.mat']); load([subj_name{s} '_combPrep_subcortical.mat'])
            load([subj_name{s} '_spatMov_subcortical.mat']); load([subj_name{s} '_spatPrep_subcortical.mat'])
            load([subj_name{s} '_tempMov_subcortical.mat']); load([subj_name{s} '_tempPrep_subcortical.mat'])
            load([subj_name{s} '_intMov_subcortical.mat']); load([subj_name{s} '_intPrep_subcortical.mat'])
            
            combMovGroup = [combMovGroup combMov(:,2)]; combPrepGroup = [combPrepGroup combPrep(:,2)];
            spatMovGroup = [spatMovGroup spatMov(:,2)]; spatPrepGroup = [spatPrepGroup spatPrep(:,2)];
            tempMovGroup = [tempMovGroup tempMov(:,2)]; tempPrepGroup = [tempPrepGroup tempPrep(:,2)];
            intMovGroup = [intMovGroup intMov(:,2)]; intPrepGroup = [intPrepGroup intPrep(:,2)];
            
        end
        
        cd(subcorticalAreaDir)
        save('combMov_subcortical', 'combMovGroup'); save('combPrep_subcortical', 'combPrepGroup')
        save('spatMov_subcortical', 'spatMovGroup'); save('spatPrep_subcortical', 'spatPrepGroup')
        save('tempMov_subcortical', 'tempMovGroup'); save('tempPrep_subcortical', 'tempPrepGroup')
        save('intMov_subcortical', 'intMovGroup'); save('intPrep_subcortical', 'intPrepGroup')
        
    case 'MVA_zacc_subcorticalArea'
        
        %%% 25% chance
        cd(subcorticalAreaDir);
        
        load('combMov_subcortical.mat'); load('combPrep_subcortical.mat')
        load('intMov_subcortical.mat'); load('intPrep_subcortical.mat')
        
        numTests=6;
        numCat=4;
        mu=1/numCat; %mu=0.25;
        N=numTests*numCat;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        %         sprintf('(X./X).*((X-%d)/%d)',mu, sigma); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
        combMovGroupZacc = (cell2mat(combMovGroup(:,2:24))-mu)/sigma;
        combPrepGroupZacc = (cell2mat(combPrepGroup(:,2:24))-mu)/sigma;
        intMovGroupZacc = (cell2mat(intMovGroup(:,2:24))-mu)/sigma;
        intPrepGroupZacc = (cell2mat(intPrepGroup(:,2:24))-mu)/sigma;
        
        combMovGroupZacc = [combMovGroup(:,1), num2cell(combMovGroupZacc)]';
        combPrepGroupZacc = [combPrepGroup(:,1), num2cell(combPrepGroupZacc)]';
        intMovGroupZacc = [intMovGroup(:,1), num2cell(intMovGroupZacc)]';
        intPrepGroupZacc = [intPrepGroup(:,1), num2cell(intPrepGroupZacc)]';
        
        save('zacc_combMov_subcortical', 'combMovGroupZacc'); save('zacc_combPrep_subcortical', 'combPrepGroupZacc')
        save('zacc_intMov_subcortical', 'intMovGroupZacc'); save('zacc_intPrep_subcortical', 'intPrepGroupZacc')
        
        xlswrite('zacc_combMov_subcortical', combMovGroupZacc); xlswrite('zacc_combPrep_subcortical', combPrepGroupZacc)
        xlswrite('zacc_intMov_subcortical', intMovGroupZacc); xlswrite('zacc_intPrep_subcortical', intPrepGroupZacc)
        
        %%% 50% chance
        cd(subcorticalAreaDir);
        
        load('spatMov_subcortical.mat'); load('spatPrep_subcortical.mat')
        load('tempMov_subcortical.mat'); load('tempPrep_subcortical.mat')
        
        takeOneOutIter=2;
        numTests=6;
        numCat=2;
        mu=1/numCat; %mu=0.5;
        N=numTests*numCat*takeOneOutIter;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        %         sprintf('(X./X).*((X-%d)/%d)',mu, sigma); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
        spatMovGroupZacc = (cell2mat(spatMovGroup(:,2:24))-mu)/sigma;
        spatPrepGroupZacc = (cell2mat(spatPrepGroup(:,2:24))-mu)/sigma;
        tempMovGroupZacc = (cell2mat(tempMovGroup(:,2:24))-mu)/sigma;
        tempPrepGroupZacc = (cell2mat(tempPrepGroup(:,2:24))-mu)/sigma;
        
        spatMovGroupZacc = [spatMovGroup(:,1), num2cell(spatMovGroupZacc)]';
        spatPrepGroupZacc = [spatPrepGroup(:,1), num2cell(spatPrepGroupZacc)]';
        tempMovGroupZacc = [tempMovGroup(:,1), num2cell(tempMovGroupZacc)]';
        tempPrepGroupZacc = [tempPrepGroup(:,1), num2cell(tempPrepGroupZacc)]';
        
        save('zacc_spatMov_subcortical', 'spatMovGroupZacc'); save('zacc_spatPrep_subcortical', 'spatPrepGroupZacc')
        save('zacc_tempMov_subcortical', 'tempMovGroupZacc'); save('zacc_tempPrep_subcortical', 'tempPrepGroupZacc')
        
        xlswrite('zacc_spatMov_subcortical', spatMovGroupZacc); xlswrite('zacc_spatPrep_subcortical', spatPrepGroupZacc)
        xlswrite('zacc_tempMov_subcortical', tempMovGroupZacc); xlswrite('zacc_tempPrep_subcortical', tempPrepGroupZacc)
        
        
        
        
        
        
    case 'RSA_dataPrep' %%%%%%%%%%%%%%%%%% Beginning of RSA analysis %%%%%%%%%%%%%%%%%%%%%%%
        
        %standard function to convert fMRI data to a format compatible with RSA toolbox
        userOptions = prepProdRSA_defineUserOptions(); %this function must be modified for each project
        betaCorrespondence = prepProdRSA_betaCorrespondence(); %likewise
        
        fullBrainVols = rsa.fmri.fMRIDataPreparation(betaCorrespondence, userOptions);
        binaryMasks_nS = rsa.fmri.fMRIMaskPreparation(userOptions);
        rsa.fmri.fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence, userOptions); %very long
        
    case 'RSA_calculateRDMs'
        cd([baseDir filesep 'imaging'])
        load([baseDir '\imaging\ImageData\prepProd_ImageData.mat']); load([baseDir '\imaging\ImageData\prepProd_Masks.mat'])
        load([baseDir '\imaging\ImageData\prepProd_responsePatterns.mat']);
        
        userOptions = prepProdRSA_defineUserOptions();
        betaCorrespondence = prepProdRSA_betaCorrespondence();
        
        RDMs  = rsa.constructRDMs(responsePatterns, betaCorrespondence, userOptions);
        sRDMs = rsa.rdm.averageRDMs_subjectSession(RDMs, 'session');
        RDMs  = rsa.rdm.averageRDMs_subjectSession(RDMs, 'session', 'subject');
        
        Models = rsa.constructModelRDMs(modelRDMs(), userOptions);
        rsa.figureRDMs(RDMs, userOptions, struct('fileName', 'RoIRDMs', 'figureNumber', 1));
        %         rsa.figureRDMs(Models, userOptions, struct('fileName', 'ModelRDMs', 'figureNumber', 2));
        
        
        
    case 'PCM_genModels' %%%%%%%%%%%%%%%%%% Beginning of PCM analysis %%%%%%%%%%%%%%%%%%
        
        prepProdSimu_genmodels('ordPrep')
        prepProdSimu_genmodels('ordProd')
        prepProdSimu_genmodels('ordSwitch')
        
        prepProdSimu_genmodels('tempPrep')
        prepProdSimu_genmodels('tempProd')
        prepProdSimu_genmodels('tempSwitch')
        
        prepProdSimu_genmodels('intPrep')
        prepProdSimu_genmodels('intProd')
        prepProdSimu_genmodels('intSwitch')
        
    case 'PCM_fitModels_PrepxProd' %main PCM script - fits models to data (4 sequences x prepProd)
        
        models = prepProdSimu_loadModels(modelDir);
        
        %         anaSubj = [3, 5, 6, 7, 9];
        
        for r=1:length(subcortStructs)
            %         for r=6;
            
            disp(['analysing ' subcortStructs{r}])
            
            %load subcortical data for analysis
            loopCounter = 1;
            for s=anaSubj
                
                cd(fullfile(glmDir, subj_name{s})); %load SPM data
                load SPM;
                nrruns=length(SPM.nscan);
                
                runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
                prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
                prep=[repmat(prep,1,nrruns) runBSL];
                
                prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
                prod=[repmat(prod,1,nrruns) runBSL];
                
                c=repmat(1:8,1,nrruns)'; % extract conditions (leave out errors!)
                run=[1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6]'; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
                
                [M, R] = prepProdSimu_prepModels(models,run,c);
                
                [~,prepCol] = find(prep>0); prepCol = reshape(prepCol,4,[]);
                [~,prodCol] = find(prod>0); prodCol = reshape(prodCol,4,[]);
                col = [prepCol; prodCol]; col = reshape(col,1,[]);
                
                
                P={SPM.Vbeta(1:126).fname}';
                Pselect=[];
                for i=1:length(col)
                    Pselect{i,1}= P{col(i)};
                end
                
                S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{r} '.mat']));
                Vin = spm_vol(char(Pselect));
                
                %extract voxels in decoding area from beta files and pass to MVPA function
                linVox=unique(cat(2,S.LI{:})');
                [I,J,K]=ind2sub(Vin(1).dim,linVox);
                X = zeros(length(linVox),length(Vin));
                
                for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                    X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
                end
                
                [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
                nanidx = unique(nanidx);
                X(nanidx,:) = [];
                Y{:,loopCounter}=X';
                
                loopCounter = loopCounter+1;
            end%subj loop
            
            %spatial pre-whitening
            %             for i=1:length(Y)
            %                 Y{i}(1:8,:) = prewhiten(Y{i}(1:8,:));
            %                 Y{i}(9:16,:) = prewhiten(Y{i}(9:16,:));
            %                 Y{i}(17:24,:) = prewhiten(Y{i}(17:24,:));
            %                 Y{i}(25:32,:) = prewhiten(Y{i}(25:32,:));
            %                 Y{i}(33:40,:) = prewhiten(Y{i}(33:40,:));
            %                 Y{i}(41:48,:) = prewhiten(Y{i}(41:48,:));
            %             end
            
            
            for i=1:length(Y)
                G_hat(:,:,i)=pcm_estGCrossval(Y{:,i},run,c);
            end;
            
            Gm = mean(G_hat,3); % Mean estimate
            
            figure
            subplot(3,7,[1, 8]);
            H = eye(8)-ones(8)/8;
            imagesc(H*Gm*H');
            title('Empirical Data')
            
            C= pcm_indicatorMatrix('allpairs',[1:8]');
            [COORD,l]=pcm_classicalMDS(Gm,'contrast',C);
            subplot(3,7,15);
            plot(COORD(:,1),COORD(:,2),'o');
            axis equal;
            
            % visualise models
            subplot(3,7,2);
            imagesc(R.orderPrepModel.G);
            title('Order Prep')
            
            subplot(3,7,9);
            imagesc(R.orderProdModel.G);
            title('Order Prod')
            
            subplot(3,7,16);
            imagesc(R.orderMaintModel.G);
            title('Order Maint')
            
            subplot(3,7,3);
            imagesc(R.timingPrepModel.G);
            title('Timing Prep')
            
            subplot(3,7,10);
            imagesc(R.timingProdModel.G);
            title('Timing Prod')
            
            subplot(3,7,17);
            imagesc(R.timingMaintModel.G);
            title('Timing Maint')
            
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
            
            % Fit the models on the group level
            [Tgroup,theta] = pcm_fitModelGroup(Y,M,run,c,'runEffect',runEffect,'fitScale',1);
            
            % Fit the models through cross-subject crossvalidation
            [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Y,M,run,c,'runEffect',runEffect,'groupFit',theta,'fitScale',1);
            
            % Provide a plot of the crossvalidated likelihoods
            subplot(3,7,[5 6 7 12 13 14 19 20 21]);
            T = pcm_plotModelLikelihood_RY(Tcross,M,'upperceil',Tgroup.likelihood(:,length(M)));
            title(subcortStructs{r})
            
        end
        
    case 'PCM_fitModels_Prep' %PCM during prep phase only (4 conditions)
        
        %load 9 representational models
        cd(modelDir)
        files = dir('*.mat');
        for i=1:length(files)
            load(files(i).name)
            models.(D.modelName) = D;
        end
        clear('D')
        
        %load subcortical data for analysis
        %         s=varargin{1};
        loopCounter = 1;
        for s=anaSubj
            
            cd(fullfile(glmDir, subj_name{s})); %load SPM data
            load SPM;
            nrruns=length(SPM.nscan);
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns)'; % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]'; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            [~,col] = find(prep>0); %identify relevant beta files and assign to pselect
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{6} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            Y{:,loopCounter}=X';
            
            loopCounter = loopCounter+1;
        end%subj loop
        
        for i=1:length(Y)
            G_hat(:,:,i)=pcm_estGCrossval(Y{:,i},run,c);
        end;
        
        Gm = mean(G_hat,3); % Mean estimate
        
        timingModel.G = pcm_estGCrossval(models.tempProd.data,run,c);
        orderModel.G = pcm_estGCrossval(models.ordProd.data,run,c);
        integratedModel.G = pcm_estGCrossval(models.intProd.data,run,c);
        
        figure
        subplot(2,4,1);
        H = eye(4)-ones(4)/4;
        imagesc(H*Gm*H');
        title('Empirical Data')
        
        C= pcm_indicatorMatrix('allpairs',[1:4]');
        [COORD,l]=pcm_classicalMDS(Gm,'contrast',C);
        subplot(2,4,2);
        plot(COORD(:,1),COORD(:,2),'o');
        axis equal;
        
        % visualise models
        subplot(2,4,5);
        imagesc(orderModel.G);
        title('Order control')
        
        subplot(2,4,6);
        imagesc(timingModel.G);
        title('Timing control')
        
        subplot(2,4,7);
        imagesc(integratedModel.G);
        title('Integrated control')
        
        % ----------------------------------------------------------------
        % Now build the models
        % Model 1: Null model for baseline: here we use a model which has all finger
        % Patterns be independent - i.e. all finger pairs are equally far away from
        % each other
        M{1}.type       = 'component';
        M{1}.numGparams = 1;
        M{1}.Gc         = ones(4);
        M{1}.name       = 'null';
        
        % Model 2: Order model, derived from simulations with
        % high order decoding
        M{2}.type       = 'component';
        M{2}.numGparams = 1;
        M{2}.Gc         = orderModel.G;
        M{2}.name       = 'order';
        
        % Model 3: Timing model, derived from simulations with
        % high timing decoding
        M{3}.type       = 'component';
        M{3}.numGparams = 1;
        M{3}.Gc         = timingModel.G;
        M{3}.name       = 'timing';
        
        % Model 4: Integrated model, derived from simulations with
        % high integrated decoding
        M{4}.type       = 'component';
        M{4}.numGparams = 1;
        M{4}.Gc         = integratedModel.G;
        M{4}.name       = 'integrated';
        
        % Model 5: Additive mixture between order and timing models
        M{5}.type       = 'component';
        M{5}.numGparams = 2;
        M{5}.Gc(:,:,1)  = orderModel.G;
        M{5}.Gc(:,:,2)  = timingModel.G;
        M{5}.name       = 'order + timing';
        
        % Model 6: Free model as Noise ceiling
        M{6}.type       = 'freechol';
        M{6}.numCond    = 4;
        M{6}.name       = 'noiseceiling';
        M{6}            = pcm_prepFreeModel(M{6});
        
        % Treat the run effect as random or fixed?
        % We are using a fixed run effect here, as we are not interested in the
        % activity relative the the baseline (rest) - so as in RSA, we simply
        % subtract out the mean patttern across all conditions.
        runEffect  = 'fixed';
        
        % Fit the models on the group level
        [Tgroup,theta] = pcm_fitModelGroup(Y,M,run,c,'runEffect',runEffect,'fitScale',1);
        
        % Fit the models through cross-subject crossvalidation
        [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Y,M,run,c,'runEffect',runEffect,'groupFit',theta,'fitScale',1);
        
        % Provide a plot of the crossvalidated likelihoods
        subplot(2,4,[4 8]);
        T = pcm_plotModelLikelihood(Tcross,M,'upperceil',Tgroup.likelihood(:,6));
        title('Preparation')
        
    case 'PCM_fitModels_Prod' %PCM during prod phase only (4 conditions)
        
        %load 9 representational models
        cd(modelDir)
        files = dir('*.mat');
        for i=1:length(files)
            load(files(i).name)
            models.(D.modelName) = D;
        end
        clear('D')
        
        %load subcortical data for analysis
        %         s=varargin{1};
        loopCounter = 1;
        for s=anaSubj
            
            cd(fullfile(glmDir, subj_name{s})); %load SPM data
            load SPM;
            nrruns=length(SPM.nscan);
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns)'; % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]'; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            [~,col] = find(prod>0); %identify relevant beta files and assign to pselect
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            S = load(fullfile(subcorticalAreaDir,subj_name{s},[subj_name{s} '_decodeArea_' subcortStructs{6} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            Y{:,loopCounter}=X';
            
            loopCounter = loopCounter+1;
        end%subj loop
        
        for i=1:length(Y)
            G_hat(:,:,i)=pcm_estGCrossval(Y{:,i},run,c);
        end;
        
        Gm = mean(G_hat,3); % Mean estimate
        
        timingModel.G = pcm_estGCrossval(models.tempProd.data,run,c);
        orderModel.G = pcm_estGCrossval(models.ordProd.data,run,c);
        integratedModel.G = pcm_estGCrossval(models.intProd.data,run,c);
        
        figure
        subplot(2,4,1);
        H = eye(4)-ones(4)/4;
        imagesc(H*Gm*H');
        title('Empirical Data')
        
        C= pcm_indicatorMatrix('allpairs',[1:4]');
        [COORD,l]=pcm_classicalMDS(Gm,'contrast',C);
        subplot(2,4,2);
        plot(COORD(:,1),COORD(:,2),'o');
        axis equal;
        
        % visualise models
        subplot(2,4,5);
        imagesc(orderModel.G);
        title('Order control')
        
        subplot(2,4,6);
        imagesc(timingModel.G);
        title('Timing control')
        
        subplot(2,4,7);
        imagesc(integratedModel.G);
        title('Integrated control')
        
        % ----------------------------------------------------------------
        % Now build the models
        % Model 1: Null model for baseline: here we use a model which has all finger
        % Patterns be independent - i.e. all finger pairs are equally far away from
        % each other
        M{1}.type       = 'component';
        M{1}.numGparams = 1;
        M{1}.Gc         = ones(4);
        M{1}.name       = 'null';
        
        % Model 2: Order model, derived from simulations with
        % high order decoding
        M{2}.type       = 'component';
        M{2}.numGparams = 1;
        M{2}.Gc         = orderModel.G;
        M{2}.name       = 'order';
        
        % Model 3: Timing model, derived from simulations with
        % high timing decoding
        M{3}.type       = 'component';
        M{3}.numGparams = 1;
        M{3}.Gc         = timingModel.G;
        M{3}.name       = 'timing';
        
        % Model 4: Integrated model, derived from simulations with
        % high integrated decoding
        M{4}.type       = 'component';
        M{4}.numGparams = 1;
        M{4}.Gc         = integratedModel.G;
        M{4}.name       = 'integrated';
        
        % Model 5: Additive mixture between order and timing models
        M{5}.type       = 'component';
        M{5}.numGparams = 2;
        M{5}.Gc(:,:,1)  = orderModel.G;
        M{5}.Gc(:,:,2)  = timingModel.G;
        M{5}.name       = 'order + timing';
        
        % Model 6: Free model as Noise ceiling
        M{6}.type       = 'freechol';
        M{6}.numCond    = 4;
        M{6}.name       = 'noiseceiling';
        M{6}            = pcm_prepFreeModel(M{6});
        
        % Treat the run effect as random or fixed?
        % We are using a fixed run effect here, as we are not interested in the
        % activity relative the the baseline (rest) - so as in RSA, we simply
        % subtract out the mean patttern across all conditions.
        runEffect  = 'fixed';
        
        % Fit the models on the group level
        [Tgroup,theta] = pcm_fitModelGroup(Y,M,run,c,'runEffect',runEffect,'fitScale',1);
        
        % Fit the models through cross-subject crossvalidation
        [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Y,M,run,c,'runEffect',runEffect,'groupFit',theta,'fitScale',1);
        
        % Provide a plot of the crossvalidated likelihoods
        subplot(2,4,[4 8]);
        T = pcm_plotModelLikelihood(Tcross,M,'upperceil',Tgroup.likelihood(:,6));
        title('Production')
        
    case 'PCM_normaliseBetaMaps' %normalises individual beta maps for cortical ROI analysis
        
        s = varargin{1};
        
        if ~isdir(pcmGroupDir)
            mkdir(pcmGroupDir) % folder for pcm group dir
        end
        
        load([glmDir, '\', subj_name{s}, '\', 'SPM.mat'])
        
        nrruns=length(SPM.nscan);
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        [~,prepCol] = find(prep>0); prepCol = reshape(prepCol,4,[]);
        [~,prodCol] = find(prod>0); prodCol = reshape(prodCol,4,[]);
        col = [prepCol; prodCol]; col = reshape(col,1,[]);
        
        P={SPM.Vbeta(1:126).fname}';
        P = P(col);
        %MVPA accuracy maps
        images= P; % the above loads all regressors for prep and prod
        images = strrep(images,'.nii','');%remove .nii from filename
        defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
        
        sn_images = strcat([glmDir, '\', subj_name{s}, '\'], images, '.nii');
        out_images = strcat([pcmGroupDir, '\'], images, ['_', subj_name{s}, '.nii']);
        spmj_normalization_write(defor, sn_images,'outimages',out_images); %Trilinear interpolation
        
    case 'PCM_corticalROItoNii' %Generates elife ROIs with marsbar and converts to .nii
        
        %%% from Elife 2014 bilateral:
        %         roiOUT=[baseDir, '\imaging\ROI\roi_Elife2014.mat'];
        
        refImageDir = [groupDir '\PCM']; %directories for functional reference image (beta 1)
        cortROIDir = [baseDir, '\imaging\ROI\']; %and for cortical ROIs
        cd(cortROIDir)
        
        mni=[-36 -22 53 ; 36 -22 53 ; -21 -12 60  ; 21 -12 60; -8 12 57;  8 18 49; -32 -54 56; 30 -59 46]; %peaks of overall classifier in elife
        roiName={'LM1','RM1','LPMd','RPMd','LSMA','RSMA','LSPC','RSPC'}; %corresponding names
        
        for i=1:length(roiName)
            
            sphereROI = maroi_sphere(struct('centre',mni(i,:),'radius',6,'label',roiName{i})); %constructs ROIs in marsbar format
            
            %%% Then save:
            saveroi(sphereROI, [roiName{i}, '_elife.mat']); %as .mat
            save_as_image(sphereROI, [roiName{i} '_elife.nii']) %as .nii too
            
            %then reslice into functional resolution:
            matlabbatch{1}.spm.spatial.realign.write.data = {[refImageDir, '\beta_0001_s03.nii']; [roiName{i}, '_elife.nii']};
            matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
            
            spm_jobman('run',matlabbatch);
            disp([roiName{i} ' resliced'])
            
        end%for rois
        
    case 'PCM_cortical_define_mva_area' %converts RoI from hand-drawn .nii to voxel coordinate .mat for PCM
        
        rois={'LM1','RM1','LPMd','RPMd','LSMA','RSMA','LSPC','RSPC'}; %corresponding names
        cd(roiDir);
        
        % load respective subcortical .nii and extract x y z coordinates for voxels
        for r=1:length(rois)
            V=spm_vol(['r' rois{r} '_elife.nii']);
            X=spm_read_vols(V);
            [i,j,k]=ind2sub(size(X),find(X~=0));
            vox=[i j k];
            
            LI{1} = sub2ind(V.dim,vox(:,1),vox(:,2),vox(:,3)); %convert x y z coordinates into linear index for voxels
            voxmin = min(vox,[],1); %min x y z coordinates for each voxel
            voxmax = max(vox,[],1); %max as above
            n = size(vox,1); %number of voxels
            
            save([pcmDir '\decodeArea_' rois{r} '.mat'], 'vox', 'LI', 'voxmin', 'voxmax', 'n')
            disp([rois{r} ' done'])
        end
        
    case 'PCM_fitModels_corticalROI' %fits PCM models to data from hand-drawn ROIs
        
        rois={'LM1','RM1','LPMd','RPMd','LSMA','RSMA','LSPC','RSPC'}; %corresponding names
        
        %load 9 representational models
        cd(modelDir)
        files = dir('*.mat');
        for i=1:length(files)
            load(files(i).name)
            models.(D.modelName) = D;
        end
        clear('D')
        
        for r=1:length(rois)%for each roi
            
            %load subcortical data for analysis
            loopCounter = 1;
            for s=anaSubj
                %         for s=3
                
                cd(fullfile(glmDir, subj_name{s})); %load SPM data
                load SPM;
                nrruns=length(SPM.nscan);
                
                runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
                prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
                prep=[repmat(prep,1,nrruns) runBSL];
                
                prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
                prod=[repmat(prod,1,nrruns) runBSL];
                
                c=repmat(1:8,1,nrruns)'; % extract conditions (leave out errors!)
                run=[1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6]'; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
                
                [~,prepCol] = find(prep>0); prepCol = reshape(prepCol,4,[]);
                [~,prodCol] = find(prod>0); prodCol = reshape(prodCol,4,[]);
                col = [prepCol; prodCol]; col = reshape(col,1,[]);
                
                
                P={SPM.Vbeta(1:126).fname}';
                Pselect=[];
                for i=1:length(col)
                    Pselect{i,1}= P{col(i)};
                end
                
                S = load(fullfile(pcmDir,['\decodeArea_' rois{r} '.mat']));
                
                Pselect = strcat([groupDir, '\PCM\'], strrep(Pselect,'.nii',''), '_', subj_name{s}, '.nii');
                Vin = spm_vol(char(Pselect));
                
                %extract voxels in decoding area from beta files and pass to MVPA function
                linVox=unique(cat(2,S.LI{:})');
                [I,J,K]=ind2sub(Vin(1).dim,linVox);
                X = zeros(length(linVox),length(Vin));
                
                for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                    X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
                end
                
                [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
                nanidx = unique(nanidx);
                X(nanidx,:) = [];
                Y{:,loopCounter}=X';
                
                loopCounter = loopCounter+1;
            end%subj loop
            
            for i=1:length(Y)
                G_hat(:,:,i)=pcm_estGCrossval(Y{:,i},run,c);
            end;
            
            Gm = mean(G_hat,3); % Mean estimate
            
            orderPrepModel.G = pcm_estGCrossval(models.ordPrep.data,run,c); orderPrepModel.G(isnan(orderPrepModel.G)) = 0;
            timingPrepModel.G = pcm_estGCrossval(models.tempPrep.data,run,c); timingPrepModel.G(isnan(timingPrepModel.G)) = 0;
            integratedPrepModel.G = pcm_estGCrossval(models.intPrep.data,run,c); integratedPrepModel.G(isnan(integratedPrepModel.G)) = 0;
            
            orderProdModel.G = pcm_estGCrossval(models.ordProd.data,run,c); orderProdModel.G(isnan(orderProdModel.G)) = 0;
            timingProdModel.G = pcm_estGCrossval(models.tempProd.data,run,c); timingProdModel.G(isnan(timingProdModel.G)) = 0;
            integratedProdModel.G = pcm_estGCrossval(models.intProd.data,run,c); integratedProdModel.G(isnan(integratedProdModel.G)) = 0;
            
            orderMaintModel.G = pcm_estGCrossval(models.ordMaint.data,run,c); orderMaintModel.G(isnan(orderMaintModel.G)) = 0;
            timingMaintModel.G = pcm_estGCrossval(models.tempMaint.data,run,c); timingMaintModel.G(isnan(timingMaintModel.G)) = 0;
            integratedMaintModel.G = pcm_estGCrossval(models.intMaint.data,run,c); integratedMaintModel.G(isnan(integratedMaintModel.G)) = 0;
            
            figure
            subplot(3,7,[1, 8]);
            H = eye(8)-ones(8)/8;
            imagesc(H*Gm*H');
            title('Empirical Data')
            
            C= pcm_indicatorMatrix('allpairs',[1:8]');
            [COORD,l]=pcm_classicalMDS(Gm,'contrast',C);
            subplot(3,7,15);
            plot(COORD(:,1),COORD(:,2),'o');
            axis equal;
            
            % visualise models
            subplot(3,7,2);
            imagesc(orderPrepModel.G);
            title('Order Prep')
            
            subplot(3,7,9);
            imagesc(orderProdModel.G);
            title('Order Prod')
            
            subplot(3,7,16);
            imagesc(orderMaintModel.G);
            title('Order Maint')
            
            subplot(3,7,3);
            imagesc(timingPrepModel.G);
            title('Timing Prep')
            
            subplot(3,7,10);
            imagesc(timingProdModel.G);
            title('Timing Prod')
            
            subplot(3,7,17);
            imagesc(timingMaintModel.G);
            title('Timing Maint')
            
            subplot(3,7,4);
            imagesc(integratedPrepModel.G);
            title('Integrated Prep')
            
            subplot(3,7,11);
            imagesc(integratedProdModel.G);
            title('Integrated Prod')
            
            subplot(3,7,18);
            imagesc(integratedMaintModel.G);
            title('Integrated Maint')
            
            % ----------------------------------------------------------------
            % ----------------------------------------------------------------
            % Now build the models
            % Model 1: Null model for baseline: here we use a model which has all finger
            % Patterns be independent - i.e. all finger pairs are equally far away from
            % each other
            M{1}.type       = 'component';
            M{1}.numGparams = 1;
            M{1}.Gc         = ones(8);
            M{1}.name       = 'null';
            
            % Models 2-4: Order models, derived from simulations with
            % high order decoding during prep, prod, & switch respectively
            M{2}.type       = 'component';
            M{2}.numGparams = 1;
            M{2}.Gc         = orderPrepModel.G;
            M{2}.name       = 'ordPrep';
            
            M{3}.type       = 'component';
            M{3}.numGparams = 1;
            M{3}.Gc         = orderProdModel.G;
            M{3}.name       = 'ordProd';
            
            M{4}.type       = 'component';
            M{4}.numGparams = 1;
            M{4}.Gc         = orderMaintModel.G;
            M{4}.name       = 'ordSwitch';
            
            % Models 5-7: Timing models, derived from simulations with
            % high timing decoding during prep, prod, & switch respectively
            M{5}.type       = 'component';
            M{5}.numGparams = 1;
            M{5}.Gc         = timingPrepModel.G;
            M{5}.name       = 'timPrep';
            
            M{6}.type       = 'component';
            M{6}.numGparams = 1;
            M{6}.Gc         = timingProdModel.G;
            M{6}.name       = 'timProd';
            
            M{7}.type       = 'component';
            M{7}.numGparams = 1;
            M{7}.Gc         = timingMaintModel.G;
            M{7}.name       = 'timSwitch';
            
            % Models 8-10: Integrated models, derived from simulations with
            % high integrated decoding
            M{8}.type       = 'component';
            M{8}.numGparams = 1;
            M{8}.Gc         = integratedPrepModel.G;
            M{8}.name       = 'intPrep';
            
            M{9}.type       = 'component';
            M{9}.numGparams = 1;
            M{9}.Gc         = integratedProdModel.G;
            M{9}.name       = 'intProd';
            
            M{10}.type       = 'component';
            M{10}.numGparams = 1;
            M{10}.Gc         = integratedMaintModel.G;
            M{10}.name       = 'intSwitch';
            
            % Model 11: Additive mixture between timing prep and prod models
            M{11}.type       = 'component';
            M{11}.numGparams = 2;
            M{11}.Gc(:,:,1)  = orderPrepModel.G;
            M{11}.Gc(:,:,2)  = orderProdModel.G;
            M{11}.name       = 'ordPrep+Prod';
            
            % Model 12: Additive mixture between order prep and prod models
            M{12}.type       = 'component';
            M{12}.numGparams = 2;
            M{12}.Gc(:,:,1)  = timingPrepModel.G;
            M{12}.Gc(:,:,2)  = timingProdModel.G;
            M{12}.name       = 'timPrep+Prod';
            
            % Model 12: Additive mixture between integrated prep and prod models
            M{13}.type       = 'component';
            M{13}.numGparams = 2;
            M{13}.Gc(:,:,1)  = integratedPrepModel.G;
            M{13}.Gc(:,:,2)  = integratedProdModel.G;
            M{13}.name       = 'intPrep+Prod';
            
            % Model 14: Free model as Noise ceiling
            M{14}.type       = 'freedirect';
            %         M{14}.type       = 'freechol';
            M{14}.numCond    = 8;
            M{14}.name       = 'noiseceiling';
            M{14}            = pcm_prepFreeModel(M{14});
            
            % Treat the run effect as random or fixed?
            % We are using a fixed run effect here, as we are not interested in the
            % activity relative the the baseline (rest) - so as in RSA, we simply
            % subtract out the mean patttern across all conditions.
            runEffect  = 'fixed';
            
            % Fit the models on the group level
            [Tgroup,theta] = pcm_fitModelGroup(Y,M,run,c,'runEffect',runEffect,'fitScale',1);
            
            % Fit the models through cross-subject crossvalidation
            [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Y,M,run,c,'runEffect',runEffect,'groupFit',theta,'fitScale',1);
            
            % Provide a plot of the crossvalidated likelihoods
            subplot(3,7,[5 6 7 12 13 14 19 20 21]);
            T = pcm_plotModelLikelihood(Tcross,M,'upperceil',Tgroup.likelihood(:,length(M)));
            title(rois{r})
        end
        
    case 'PCM_subcorticalSearchLight' % Define the search lights for the PCM analysis
        
        s=varargin{1};
        
        if ~isdir([pcmDir, '\searchlight'])
            mkdir([pcmDir, '\searchlight'])
        end
        
        for r = 1:length(subcortStructs);
            
            radius=16;
            numVox=160;
            cd(fullfile(subcorticalSearchDir, subj_name{s}));
            V=spm_vol(['mr' subj_name{s} '_' subcortStructs{r} '.nii']); %if preceded by case MVA_mask
            X=spm_read_vols(V);
            [i,j,k]=ind2sub(size(X),find(X~=0));
            vox=[i j k];
            [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[radius numVox],V.mat,V.dim,[],'mva160_numvox.nii');
            save ([pcmDir, '\searchlight\' subj_name{s} '_' subcortStructs{r} '_' 'volsearch160.mat'], 'vox', 'LI', 'voxmin', 'voxmax', 'n')
        end
        
    case 'PCM_do_subcorticalSearchLight' % Conduct the PCM in subcortical searchlights
        
        s=varargin{1};
        
        for r=1 %:length(subcortStructs);
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(pcmDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} 'PCMLikelihood.nii'])};
            
            
            % Generate column indices for Cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
                train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(subcorticalSearchDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_volsearch160.mat']),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart)
        end
        
    case 'PCM_resliceCorticalROI_byhand' %reslice hand-drawn ROIs to functional resolution
        
        sn = varargin{1};
        
        refImageDir = [glmDir, '\', subj_name{sn}]; %directories for functional reference image (beta 1)
        cortROIDir = [pcmDir, '\', subj_name{sn}]; %and for cortical ROIs
        
        rois = {'timing_premotorROI', 'integrated_parietalROI'};
        
        for r=1:length(rois) %loop through subcortical regions and reslice them to beta image resolution
            
            matlabbatch{1}.spm.spatial.realign.write.data = {[refImageDir, '\beta_0001.nii']; [cortROIDir, '\', subj_name{sn}, '_', rois{r}, '.nii']};
            matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
            
            spm_jobman('run',matlabbatch);
            disp([rois{r} ' resliced'])
        end
        
    case 'PCM_cortical_define_mva_area_byhand' %converts RoI from hand-drawn .nii to voxel coordinate .mat for PCM
        
        s=varargin{1};
        rois = {'timing_premotorROI', 'integrated_parietalROI'};
        
        % load respective subcortical .nii and extract x y z coordinates for voxels
        for r=1:length(rois)
            cd(fullfile(pcmDir, subj_name{s}));
            V=spm_vol(['r' subj_name{s} '_' rois{r} '.nii']);
            X=spm_read_vols(V);
            [i,j,k]=ind2sub(size(X),find(X~=0));
            vox=[i j k];
            
            LI{1} = sub2ind(V.dim,vox(:,1),vox(:,2),vox(:,3)); %convert x y z coordinates into linear index for voxels
            voxmin = min(vox,[],1); %min x y z coordinates for each voxel
            voxmax = max(vox,[],1); %max as above
            n = size(vox,1); %number of voxels
            
            save(sprintf('%s_decodeArea_%s.mat',subj_name{s},rois{r}), 'vox', 'LI', 'voxmin', 'voxmax', 'n')
            disp([rois{r} ' done'])
        end
        
    case 'PCM_fitModels_corticalROI_byhand' %fits PCM models to data from hand-drawn ROIs
        
        %Hand-generate cortical RoIs for each participant (?)
        %PMv, SPCa
        %Using MRIcron draw tool
        
        rois = {'timing_premotorROI', 'integrated_parietalROI'};
        
        %load 9 representational models
        cd(modelDir)
        files = dir('*.mat');
        for i=1:length(files)
            load(files(i).name)
            models.(D.modelName) = D;
        end
        clear('D')
        anaSubj = [3, 5, 6, 7, 9];
        
        %load subcortical data for analysis
        loopCounter = 1;
        for s=anaSubj
            %         for s=3
            
            cd(fullfile(glmDir, subj_name{s})); %load SPM data
            load SPM;
            nrruns=length(SPM.nscan);
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat(1:8,1,nrruns)'; % extract conditions (leave out errors!)
            run=[1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6]'; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            [~,prepCol] = find(prep>0); prepCol = reshape(prepCol,4,[]);
            [~,prodCol] = find(prod>0); prodCol = reshape(prodCol,4,[]);
            col = [prepCol; prodCol]; col = reshape(col,1,[]);
            
            
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            S = load(fullfile(pcmDir,subj_name{s},[subj_name{s} '_decodeArea_' rois{2} '.mat']));
            Vin = spm_vol(char(Pselect));
            
            %extract voxels in decoding area from beta files and pass to MVPA function
            linVox=unique(cat(2,S.LI{:})');
            [I,J,K]=ind2sub(Vin(1).dim,linVox);
            X = zeros(length(linVox),length(Vin));
            
            for i=1:length(Vin) %this gives us a voxelN x betaN matrix
                X(:,i)=spm_sample_vol(Vin(i),double(I),double(J),double(K),0);
            end
            
            [nanidx, ~] = find(isnan(X)); %remove nan values contained in beta image due to signal loss(?)
            nanidx = unique(nanidx);
            X(nanidx,:) = [];
            Y{:,loopCounter}=X';
            
            loopCounter = loopCounter+1;
        end%subj loop
        
        for i=1:length(Y)
            G_hat(:,:,i)=pcm_estGCrossval(Y{:,i},run,c);
        end;
        
        Gm = mean(G_hat,3); % Mean estimate
        
        orderPrepModel.G = pcm_estGCrossval(models.ordPrep.data,run,c); orderPrepModel.G(isnan(orderPrepModel.G)) = 0;
        timingPrepModel.G = pcm_estGCrossval(models.tempPrep.data,run,c); timingPrepModel.G(isnan(timingPrepModel.G)) = 0;
        integratedPrepModel.G = pcm_estGCrossval(models.intPrep.data,run,c); integratedPrepModel.G(isnan(integratedPrepModel.G)) = 0;
        
        orderProdModel.G = pcm_estGCrossval(models.ordProd.data,run,c); orderProdModel.G(isnan(orderProdModel.G)) = 0;
        timingProdModel.G = pcm_estGCrossval(models.tempProd.data,run,c); timingProdModel.G(isnan(timingProdModel.G)) = 0;
        integratedProdModel.G = pcm_estGCrossval(models.intProd.data,run,c); integratedProdModel.G(isnan(integratedProdModel.G)) = 0;
        
        orderSwitchModel.G = pcm_estGCrossval(models.ordSwitch.data,run,c); orderSwitchModel.G(isnan(orderSwitchModel.G)) = 0;
        timingSwitchModel.G = pcm_estGCrossval(models.tempSwitch.data,run,c); timingSwitchModel.G(isnan(timingSwitchModel.G)) = 0;
        integratedSwitchModel.G = pcm_estGCrossval(models.intSwitch.data,run,c); integratedSwitchModel.G(isnan(integratedSwitchModel.G)) = 0;
        
        figure
        subplot(3,7,[1, 8]);
        H = eye(8)-ones(8)/8;
        imagesc(H*Gm*H');
        title('Empirical Data')
        
        C= pcm_indicatorMatrix('allpairs',[1:8]');
        [COORD,l]=pcm_classicalMDS(Gm,'contrast',C);
        subplot(3,7,15);
        plot(COORD(:,1),COORD(:,2),'o');
        axis equal;
        
        % visualise models
        subplot(3,7,2);
        imagesc(orderPrepModel.G);
        title('Order Prep')
        
        subplot(3,7,9);
        imagesc(orderProdModel.G);
        title('Order Prod')
        
        subplot(3,7,16);
        imagesc(orderSwitchModel.G);
        title('Order Switch')
        
        subplot(3,7,3);
        imagesc(timingPrepModel.G);
        title('Timing Prep')
        
        subplot(3,7,10);
        imagesc(timingProdModel.G);
        title('Timing Prod')
        
        subplot(3,7,17);
        imagesc(timingSwitchModel.G);
        title('Timing Switch')
        
        subplot(3,7,4);
        imagesc(integratedPrepModel.G);
        title('Integrated Prep')
        
        subplot(3,7,11);
        imagesc(integratedProdModel.G);
        title('Integrated Prod')
        
        subplot(3,7,18);
        imagesc(integratedSwitchModel.G);
        title('Integrated Switch')
        
        % ----------------------------------------------------------------
        % Now build the models
        % Model 1: Null model for baseline: here we use a model which has all finger
        % Patterns be independent - i.e. all finger pairs are equally far away from
        % each other
        M{1}.type       = 'component';
        M{1}.numGparams = 1;
        M{1}.Gc         = ones(8);
        M{1}.name       = 'null';
        
        % Models 2-4: Order models, derived from simulations with
        % high order decoding during prep, prod, & switch respectively
        M{2}.type       = 'component';
        M{2}.numGparams = 1;
        M{2}.Gc         = orderPrepModel.G;
        M{2}.name       = 'orderPrep';
        
        M{3}.type       = 'component';
        M{3}.numGparams = 1;
        M{3}.Gc         = orderProdModel.G;
        M{3}.name       = 'orderProd';
        
        M{4}.type       = 'component';
        M{4}.numGparams = 1;
        M{4}.Gc         = orderSwitchModel.G;
        M{4}.name       = 'orderSwitch';
        
        % Models 5-7: Timing models, derived from simulations with
        % high timing decoding during prep, prod, & switch respectively
        M{5}.type       = 'component';
        M{5}.numGparams = 1;
        M{5}.Gc         = timingPrepModel.G;
        M{5}.name       = 'timingPrep';
        
        M{6}.type       = 'component';
        M{6}.numGparams = 1;
        M{6}.Gc         = timingProdModel.G;
        M{6}.name       = 'timingProd';
        
        M{7}.type       = 'component';
        M{7}.numGparams = 1;
        M{7}.Gc         = timingSwitchModel.G;
        M{7}.name       = 'timingSwitch';
        
        % Models 8-10: Integrated models, derived from simulations with
        % high integrated decoding
        M{8}.type       = 'component';
        M{8}.numGparams = 1;
        M{8}.Gc         = integratedPrepModel.G;
        M{8}.name       = 'integratedPrep';
        
        M{9}.type       = 'component';
        M{9}.numGparams = 1;
        M{9}.Gc         = integratedProdModel.G;
        M{9}.name       = 'integratedProd';
        
        M{10}.type       = 'component';
        M{10}.numGparams = 1;
        M{10}.Gc         = integratedSwitchModel.G;
        M{10}.name       = 'integratedSwitch';
        
        % Model 11: Additive mixture between order and timing models
        M{11}.type       = 'component';
        M{11}.numGparams = 2;
        M{11}.Gc(:,:,1)  = orderPrepModel.G;
        M{11}.Gc(:,:,2)  = orderProdModel.G;
        M{11}.name       = 'orderPrep + orderProd';
        
        % Model 12: Additive mixture between order and timing models
        M{12}.type       = 'component';
        M{12}.numGparams = 2;
        M{12}.Gc(:,:,1)  = timingPrepModel.G;
        M{12}.Gc(:,:,2)  = timingProdModel.G;
        M{12}.name       = 'timingPrep + timingProd';
        
        % Model 6: Free model as Noise ceiling
        M{13}.type       = 'freedirect';
        M{13}.numCond    = 8;
        M{13}.name       = 'noiseceiling';
        M{13}            = pcm_prepFreeModel(M{13});
        
        % Treat the run effect as random or fixed?
        % We are using a fixed run effect here, as we are not interested in the
        % activity relative the the baseline (rest) - so as in RSA, we simply
        % subtract out the mean patttern across all conditions.
        runEffect  = 'fixed';
        
        % Fit the models on the group level
        [Tgroup,theta] = pcm_fitModelGroup(Y,M,run,c,'runEffect',runEffect,'fitScale',1);
        
        % Fit the models through cross-subject crossvalidation
        [Tcross,thetaCr] = pcm_fitModelGroupCrossval(Y,M,run,c,'runEffect',runEffect,'groupFit',theta,'fitScale',1);
        
        % Provide a plot of the crossvalidated likelihoods
        subplot(3,7,[5 6 7 12 13 14 19 20 21]);
        T = pcm_plotModelLikelihood(Tcross,M,'upperceil',Tgroup.likelihood(:,length(M)));
        %         title('Preparation')
        
        
end%switch

end%function