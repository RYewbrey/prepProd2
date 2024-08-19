function prepProd2_imana_RY(what,varargin)


%PC paths Load relevant Matlab directories RY
% addpath(genpath('Z:/rhys/prepProd2/matlab')); %Adjust! loaded with subdirectories (genpath command)
% addpath('Z:/toolboxes/spm12');
% addpath(genpath('Z:/toolboxes/tools')); %joern's extensions for spm
% addpath(genpath('Z:/toolboxes/userfun')); %joern's util tools (open source)
% addpath(genpath('Z:/toolboxes/region-master')); %joern's region toolbox for spm
% addpath(genpath('Z:/toolboxes/spm12/toolbox/suit')); %SUIT Cerebellum
% addpath(genpath('Z:/toolboxes/spm12/toolbox/DARTEL')); %DARTEL deformation for suit reslice
% addpath(genpath('Z:/toolboxes/spm12/toolbox/OldSeg')); %for suit reslice
% addpath(genpath('Z:/toolboxes/permutest')); %permutest for crossSection analysis
% addpath('Z:/toolboxes/rsatoolbox_matlab'); %RSA toolbox
% addpath(genpath('Z:/toolboxes/pcm_toolbox')); %PCM toolbox
%
%
%PC paths Load relevant Matlab directories KK
% addpath(genpath('Z:\kornyshk-kornyshevalab\rhys\prepProd2\matlab')); %Adjust! loaded with subdirectories (genpath command)
% addpath('Z:\kornyshk-kornyshevalab\toolboxes\spm12');
% addpath(genpath('Z:\kornyshk-kornyshevalab\toolboxes\tools')); %joern's extensions for spm
% addpath(genpath('Z:\kornyshk-kornyshevalab\toolboxes\userfun')); %joern's util tools (open source)
% addpath(genpath('Z:\kornyshk-kornyshevalab\toolboxes\region-master')); %joern's region toolbox for spm
% addpath(genpath('Z:\kornyshk-kornyshevalab\toolboxes\spm12/toolbox\suit')); %SUIT Cerebellum
% addpath(genpath('Z:\kornyshk-kornyshevalab\toolboxes\spm12\toolbox\DARTEL')); %DARTEL deformation for suit reslice
% addpath(genpath('Z:\kornyshk-kornyshevalab\toolboxes\pm12/toolbox\OldSeg')); %for suit reslice
% addpath(genpath('Z:\kornyshk-kornyshevalab\toolboxes\permutest')); %permutest for crossSection analysis
% addpath('Z:\kornyshk-kornyshevalab\toolboxes\rsatoolbox_matlab'); %RSA toolbox
% addpath(genpath('Z:\kornyshk-kornyshevalab\toolboxes\pcm_toolbox')); %PCM toolbox
%
%
%BlueBear paths
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/rhys/prepProd2/matlab')); %Adjust! loaded with subdirectories (genpath command)
% addpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/spm12');
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/tools')); %joern's extensions for spm
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/userfun')); %joern's util tools (open source)
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/region-master')); %joern's region toolbox for spm
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/spm12/toolbox/suit')); %SUIT cerebellum toolbox
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/spm12/toolbox/DARTEL')); %DARTEL transformation
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/spm12/toolbox/OldSeg')); %For SUIT reslicing
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/permutest')); %permutest for crossSection analysis
% addpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/rsatoolbox_matlab'); %RSA toolbox
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/pcm_toolbox')); %PCM toolbox
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/surfing')); %surfing, used in RSA toolbox

%%%definition and variables

%%% Data paths:
% baseDir= 'Z:/rhys/prepProd2/data'; %PC
baseDir= '/rds/projects/k/kornyshk-kornyshevalab/rhys/prepProd2/data'; %BlueBear

rawDir=[baseDir filesep 'imaging' filesep 'raw']; %original files before conversion
anatDir=[baseDir filesep 'imaging' filesep 'anatomicals'];
epiDir=[baseDir filesep 'imaging' filesep 'epi'];
glmDir=[baseDir filesep 'imaging' filesep 'GLM_firstlevel'];
behDir=[baseDir filesep 'behavioural'];
groupDir=[baseDir filesep 'imaging' filesep 'GLM_secondlevel'];
scndDir=[baseDir filesep 'imaging' filesep 'GLM_secondlevel' filesep 'data'];
roiDir=[baseDir filesep 'imaging' filesep 'ROI' filesep 'cortical'];
roiSubDir=[baseDir filesep 'imaging' filesep 'ROI' filesep 'subcortical'];
roiCbDir=[baseDir filesep 'imaging' filesep 'ROI' filesep 'cerebellum'];
roiCorticalDir=[baseDir filesep 'imaging' filesep 'ROI' filesep 'cortical'];
suitDir=[baseDir filesep 'imaging' filesep 'suit'];
suitGroupDir=[baseDir filesep 'imaging' filesep 'suit_secondlevel'];
suitGroupRSADir=[suitGroupDir filesep 'rsa'];
suitScndDir=[baseDir filesep 'imaging' filesep 'suit_secondlevel' filesep 'data'];
subcorticalDir=[baseDir filesep 'imaging' filesep 'subcortical'];
subcorticalGroupDir=[baseDir filesep 'imaging' filesep 'subcortical_secondlevel'];
subcorticalGroupRSADir=[subcorticalGroupDir filesep 'rsa'];
cerebellumDir=[baseDir filesep 'imaging' filesep 'cerebellum'];
crossSectionDir=[baseDir filesep 'imaging' filesep 'crossSection'];
modelDir=[baseDir filesep 'imaging' filesep 'simulations' filesep 'models'];
rsaDir=[baseDir filesep 'imaging' filesep 'rsa'];
rsaCorticalDir=[rsaDir filesep 'cortical'];
rsaCorticalGroupDir=[baseDir filesep 'imaging' filesep 'rsa_secondlevel' filesep 'cortical'];
subcorticalPeaksDir=[subcorticalGroupDir filesep 'cluster_peaks'];
simuDir = [baseDir filesep 'imaging' filesep 'simulations'];

%%% Names and IDs
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

scanID={...
    'placeholder', 'ph', 'ph', 'ph', 'ph','ph'; ... %s01
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

run={'1','2','3','4','5','6'};

%%% Parameters
TR=2; %sec
nrVolumes=230; %number of volumes
sn8NrVolumes=218;
nrSlices=60; %different across subjects

%Slice Acquisition
MBsliceAcquisition = (1:nrSlices)';        %multiband acquisition
multiBand=2;       %N slices captured at one time
interSlice = TR/(nrSlices/multiBand);  %for multiband

%%% Generate vector for multiband acquisition. Comment if using standard acquisition
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
originAC=[...
    1, 1, 1, 1; ...      %s01 placeholder
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
    1, 1, 1, 1; ...      %s14 placeholder
    124,-131,-83, 1; ... %s15
    121,-139,-76, 1; ... %s16
    120,-128,-74, 1; ... %s17
    127,-136,-74, 1; ... %s18
    1, 1, 1, 1; ...      %s19 placeholder
    124,-135,-75, 1; ... %s20
    128,-131,-78, 1; ... %s21
    120,-130,-77, 1; ... %s22
    1, 1, 1, 1; ...      %s23 placeholder
    1, 1, 1, 1; ...      %s24 placeholder
    120,-132,-73, 1; ... %s25
    127,-130,-83, 1; ... %s26
    1, 1, 1, 1; ...      %s27 placeholder
    1, 1, 1, 1; ...      %s28 placeholder
    1, 1, 1, 1; ...      %s29 placeholder
    1, 1, 1, 1; ...      %s30 placeholder
    126,-132,-78, 1; ... %s31
    124,-131,-78, 1; ... %s32
    124,-135,-82, 1; ... %s33
    126,-135,-80, 1; ... %s34
    1, 1, 1, 1; ...      %s35 placeholder
    121,-128,-81, 1; ... %s36
    122,-129,-78, 1; ... %s37
    127,-128,-81, 1; ... %s38
    125,-132,-73, 1; ... %s39
    124,-134,-70, 1; ... %s40
    124,-131,-77, 1; ... %s41
    125,-141,-84, 1; ... %s42
    ];

subcortStructs = {... %subcortical structures of interest for later analysis, defined from freesurfer aseg
    10,              11,             12,             17,...                    %13,
    29,              50,             51,             53;...                    %52,
    'left_thalamus', 'left_caudate', 'left_putamen', 'left_hippocampus',...    %'left_pallidum',
    'right_thalamus','right_caudate','right_putamen','right_hippocampus',...   %'right_pallidum',
    };

suitCBRegions = {...
    1,               2,                3,               4                 5,               7,                8,             10              11,            13;
    'left_lobule_4', 'right_lobule_4', 'left_lobule_5', 'right_lobule_5', 'left_lobule_6', 'right_lobule_6', 'left_crus_1', 'right_crus_1', 'left_crus_2', 'right_crus_2'};

corticalRegions = ...
    {'LM1', 'LPMd', 'LSMA', 'LSPC', 'RM1', 'RPMd', 'RSMA', 'RSPC'};

%%%
switch(what)
    
    case 'WHOLE_BRAIN_DECODING_JNeurosci' %----------------- BEGINNING OF WHOLE-BRAIN -----------------%
    case 'make_nii' %pre-processing
        sn=varargin{1};
        cd(fullfile(rawDir,subjID{sn}));
        
        %         %%Fieldmaps:
        %         try
        %             for i=1:2
        %                 source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn}  '_' fieldMapID{sn,i} '_Fieldmap-'  num2str(i) '_te' fieldMapID2{sn,i} '.nii.gz']);
        %                 dest = fullfile(rawDir, subjID{sn});
        %                 gunzip(source,dest); %unzip
        %
        %                 source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn}  '_' fieldMapID{sn,i} '_Fieldmap-'  num2str(i) '_te' fieldMapID2{sn,i} '.nii']);
        %                 dest   = fullfile(epiDir, subj_name{sn}, [subj_name{sn} '_' fieldMapName{sn,i} '.nii']); %EPI directory
        %
        %                 destFolder   = fullfile(epiDir, subj_name{sn}); %EPI directory
        %                 k=exist(destFolder);
        %                 if k==0
        %                     mkdir(destFolder)  ;
        %                 end;
        %
        %                 movefile(source,dest); %move
        %
        %             end;
        %             disp('fieldmap done')
        %         catch
        %             disp('No field map. Most likely a typo in file name or missing destination directory!');
        %         end;
        
        %%Anatomical:
        try
            if sn <= 16
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} '_202_scombi.nii.gz']); %naming convention for s1-s16
            else
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_202.nii.gz']); %naming convention changed on Odin from S17 onwards
            end
            
            if sn == 20 %SS0094_25_2 had a different number assigned to scombi scan
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_302.nii.gz']); %naming convention changed on Odin from S17 onwards
            end
            
            if sn == 40 %SS0094_25_2 had a different number assigned to scombi scan
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_302.nii.gz']); %naming convention changed on Odin from S17 onwards
            end
            
            dest = fullfile(rawDir, subjID{sn});
            gunzip(source,dest); %unzip
            
            if sn <= 16
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} '_202_scombi.nii']); %naming convention for s1-s16
            else
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_202.nii']); %naming convention changed on Odin from S17 onwards
            end
            if sn == 20 %SS0094_25_2 and SS0094_38 had a different number assigned to scombi scan
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_302.nii']); %naming convention changed on Odin from S17 onwards
            end
            
            if sn == 40 %SS0094_25_2 and SS0094_38 had a different number assigned to scombi scan
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} 'scombi1_302.nii']); %naming convention changed on Odin from S17 onwards
            end
            
            dest   = fullfile(anatDir, subj_name{sn}, [subj_name{sn} '_anatomical_orig.nii']);
            
            destFolder   = fullfile(anatDir, subj_name{sn}); %EPI directory
            k=exist(destFolder,'dir');
            if k==0
                mkdir(destFolder)  ;
            end
            
            movefile(source,dest); %move
            disp('anatomical done')
        catch
            disp('No anatomicals. Most likely a typo in file name or missing destination directory!');
        end
        
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
                k=exist(destFolder,'dir');
                if k==0
                    mkdir(destFolder)  ;
                end
                
                movefile(source,dest); %move
            end
            disp('EPI done')
        catch
            disp('No EPI. Most likely a typo in file name or missing destination directory!');
        end
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
        for r= 1:length(run) %for each volume and each run, generate a path to the file and store it in N
            
            if sn == 8 && r == 4
                nrVolumes = sn8NrVolumes;
            else
                nrVolumes = nrVolumes;
            end
            
            for i=1:nrVolumes
                N{i,1} = [fullfile(epiDir, subj_name{sn}, [ subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
            end
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
            ['ua',
            subj_name{sn}, '_run5.nii,1']
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
        
        %%% *1. Coregister meanEPI to anatomical (no reslicing into anat space, dimension of EPI preserved!): %%%
        %         NOTE: If original image completely off in terms of alignment to anat - First coregister per hand via 'coregtool' command
        %         then run the spm algorithm:
        ref=fullfile(anatDir, subj_name{sn},[subj_name{sn}, '_anatomical.nii']);
        source=fullfile(epiDir, subj_name{sn},['meanepi_' subj_name{sn} '.nii']);
        matlabbatch{1}.spm.spatial.coreg.estimate.ref =  {ref};
        matlabbatch{1}.spm.spatial.coreg.estimate.source = {source};
        matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm_jobman('run',matlabbatch); %meanepi will be overwritten by a version that is corgegistered
        
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
    case 'glm_set' %GLM
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
                end
            end
            
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
            end
            
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
            end
            
            
            
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
            end
            
            
            
            
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
                J.sess(r).cond(c).duration = [( (B.tZero(idx)- B.tZeroCue(idx))+B.goCueDur(idx)+B.crossDur(idx)+B.feedbackCalcDur(idx) + 1000 )/1000 ];
                
            else
                J.sess(r).cond(c).onset = -99;
                J.sess(r).cond(c).duration = 1;
                
            end
            
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
        end
        
        
        T=addstruct(T,S);
        %J.fact = struct('name', {}, 'levels', {});
        
        %%%Derivatives:
        J.bases.hrf.derivs = [1 0];
        J.volt = 1;
        J.global = 'None';
        
        %%J.mask = {[fullfile(anatDir,  subj_name{sn},[subj_name{sn},'_mask4glm.img'])],1};
        J.mthresh = 0.8;%By default, SPM masks out voxels that are below 20% of the global mean image value for a given subject
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
        
        %%% Define contrasts:
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
        
        %%% Do the contrasts
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
        
        con=fullfile(glmDir, subj_name{s},'perc_0001.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_sperc_0001.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        con=fullfile(glmDir, subj_name{s},'perc_0002.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_sperc_0002.nii']);
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
            'Z:/toolboxes/spm12/tpm/OldSeg/grey.nii'
            'Z:/toolboxes/spm12/tpm/OldSeg/white.nii'
            'Z:/toolboxes/spm12/tpm/OldSeg/csf.nii'
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
            
            [~,col] = find(prep>0);
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
    case 'MVA_do_Int_Mov'     %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
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
    case 'MVA_do_Int_Mov_2x2'  %New addition 2023, due to 2x2 subtraction effect
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        %%%%%%% We change C to distinguish between 2 sequences, because 2x2
        %%%%%%% design causes subtraction to make two conditions the same.
        c=repmat([1 2 2 1],1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Int_160_Mov_2x2.nii'])};
        
        
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
    case 'MVA_do_Int_Prep_2x2' %New addition 2023, due to 2x2 subtraction effect
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
        prep=[repmat(prep,1,nrruns) runBSL];
        
        %%%%%%% We change C to distinguish between 2 sequences, because 2x2
        %%%%%%% design causes subtraction to make two conditions the same.
        c=repmat([1 2 2 1],1,nrruns); % extract conditions (leave out errors!)
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(glmDir, subj_name{s}, [subj_name{s}, '_accuracy_Int_160_Prep_2x2.nii'])};
        
        
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
        end
    case 'MVA_zValue_oneOut'
        
        s=varargin{1};
        cd(fullfile(glmDir,subj_name{s}));
        
        takeOneOutIter=2;
        numTests=6;
        numCat=2;
        mu=1/numCat; %mu=0.5;
        N=numTests*numCat*takeOneOutIter;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        images= {...
            '_accuracy_Spat_160_Mov',   '_accuracy_Spat_160_Prep',...
            '_accuracy_Temp_160_Mov',   '_accuracy_Temp_160_Prep',...
            '_accuracy_Int_160_Mov_2x2','_accuracy_Int_160_Prep_2x2',...
            };
        
        outimages={...
            '_zacc_Spat_160_Mov',   '_zacc_Spat_160_Prep',...
            '_zacc_Temp_160_Mov',   '_zacc_Temp_160_Prep',...
            '_zacc_Int_160_Mov_2x2','_zacc_Int_160_Prep_2x2',...
            };
        
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},[subj_name{s} images{j} '.nii']);
            output_image= fullfile(glmDir,subj_name{s},[subj_name{s} outimages{j} '.nii']);
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
        end
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
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Int_160_Mov_2x2.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Int_160_Mov_2x2.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        comb=fullfile(glmDir, subj_name{s},[subj_name{s} '_zacc_Int_160_Prep_2x2.nii']); %%MVPA smoother
        scomb=fullfile(glmDir, subj_name{s},[subj_name{s} '_szacc_Int_160_Prep_2x2.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        %smooth other images here as required
    case 'MNI_normalization'
        mkdir(fullfile(groupDir,'data')); % folder for each contrast
        
        %MVPA accuracy maps
        images= {...
            'szacc_Comb_160_Mov.nii',   'szacc_Comb_160_Prep.nii',...
            'szacc_Spat_160_Mov.nii',   'szacc_Spat_160_Prep.nii',...
            'szacc_Temp_160_Mov.nii',   'szacc_Temp_160_Prep.nii',...
            'szacc_Int_160_Mov.nii',    'szacc_Int_160_Prep.nii',...
            'szacc_Int_160_Mov_2x2.nii','szacc_Int_160_Prep_2x2.nii',...
            }; % please add other images as required, e.g. spmT_...
        
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
    case 'MNI_normalization_perc&con'
        if ~isfolder(fullfile(groupDir,'data'))
            mkdir(fullfile(groupDir,'data')); % folder for each contrast
        end
        
        %Activity maps
        images= {'scon_0001.nii','scon_0002.nii','scon_0003.nii','scon_0004.nii','scon_0005.nii','scon_0006.nii',...
            'sperc_0001.nii','sperc_0002.nii'}; % please add other images as required, e.g. spmT_...
        
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
    case 'group_avg' %group average percent signal change maps
        
        if ~isfolder(fullfile(groupDir, 'average'))
            mkdir(fullfile(groupDir, 'average'))
        end
        cd(fullfile(groupDir, 'average'))
        
        images = {...
            'sperc_0001';'sperc_0002';...
            'szacc_Int_160_Mov_2x2'; 'szacc_Int_160_Prep_2x2';...
            };
        
        conds = {...
            'perc_mov', 'perc_prep',...
            'int_mov_2x2', 'int_prep_2x2',...
            };
        for i=1:length(images)
            loopCount = 1;
            for s = anaSubj
                Vi(loopCount) = spm_vol(fullfile(groupDir, [images{i} '_' subj_name{s} '.nii']));
                loopCount = loopCount + 1;
            end
            Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
            Vo.fname = ['avg_' conds{i} '.nii'];
            Vo.n = [1 1];
            express = 'mean(X)';
            flags.dmtx = 1;
            
            spm_imcalc(Vi, Vo, express, flags)
        end
    case 'ROI_define' %%% 1. Define ROI via GUI:
        R=region_getregions; %opens gui
        
        %%% from Elife 2014:
        %Overall-cM1; %Overall-iPMd ( x coordinate flipped to left);
        %%Overall_cSMA; %Overall_cSupParietal
        %         roiOUT='D:/projects/rhys/prepProd/data/imaging/ROI/roi_Elife2014.mat';
        %             %%% Select/Copy paste parameters as follows:
        % % %         mni= [-36 -22 53; -21 -12 60; -8 -12 57; -32 -54 56];
        % % %         radius=6; %mm
        % % %         name={'L_M1','L_PMd','L_SMA','L_SPC'};
        
        %%% from Elife 2014 bilateral:
        roiOUT=[roiDir '/roi_Elife2014_bilateral.mat'];
        
        
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
    case 'ROI_run'        %%% 2. Extract decoding values via GUI
        %         roiIN='D:/projects/rhys/prepProd/data/imaging/ROI/roi_Elife2014.mat';
        roiIN = fullfile(roiDir, 'LDA', 'roi_Elife2014_bilateral.mat');
        
        load(roiIN);
        P={};
        %%%% Comb
        cd(roiDir);
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
                end
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
                end
                ROI=addstruct(ROI,S); %add MVA type by type
                S=[];
                %            fileOut=[MVA_ID{m} '_Mov_ROI_spss.mat'];
                fileOut=[MVA_ID{m} '_' Phase_ID{phase} '_ROI_spss.mat'];
                save(fileOut,'meanROI');
            end
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
    case 'ROI_runComb'        %%% 3. Extract ovr decoding values via GUI
        %         roiIN='D:/projects/rhys/prepProd/data/imaging/ROI/roi_Elife2014.mat';
        roiIN=[roiDir '/roi_Elife2014_bilateral.mat'];
        
        load(roiIN);
        P={};
        %%%% Comb
        cd(roiDir);
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
        roiIN=[roiDir '/roi_Elife2014_bilateral.mat'];
        
        load(roiIN);
        P={};
        %%%% Comb
        cd(roiDir);
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
    case 'ROI_crossSectionPermutest' %cross-sectional analysis
        
        cd(fullfile(crossSectionDir, 'data'))
        
        crossNames = {'Premotor', 'Profile'};
        
        for i=1:length(crossNames)
            load(['PrepProd_crossSection_lh' crossNames{i} '.mat'])
            
            null = zeros(24,120); %null results to compare to, acts as one-sample T test
            
            [clusters, p_values, t_sums, permutation_distribution] = permutest(spatPrep',null','true',.05,10^4,'false');
            %             [clusters, p_values, t_sums, permutation_distribution] = permutest(spatPrep',null','true',.001,10^4,'false');
            for j=1:length(clusters)
                clusterMean = mean(spatPrep(:,clusters{j})); clusterMean = mean(clusterMean);
                clusterStd = std(mean(spatPrep(:,clusters{j}),2));
                clusterd(j) = clusterMean / clusterStd;
            end
            spatPrepClust = {clusters; p_values; clusterd; t_sums; permutation_distribution; spatPrep};
            %             save(['spatPrepClusters_lh' crossNames{i}],'spatPrepClust')
            
            [clusters, p_values, t_sums, permutation_distribution] = permutest(spatMov',null','true',.05,10^4,'false');
            %             [clusters, p_values, t_sums, permutation_distribution] = permutest(spatMov',null','true',.001,10^4,'false');
            for j=1:length(clusters)
                clusterMean = mean(spatMov(:,clusters{j})); clusterMean = mean(clusterMean);
                clusterStd = std(mean(spatMov(:,clusters{j}),2));
                clusterd(j) = clusterMean / clusterStd;
            end
            spatMovClust = {clusters; p_values; clusterd; t_sums; permutation_distribution; spatMov};
            %             save(['spatMovClusters_lh' crossNames{i}],'spatMovClust')
            
            [clusters, p_values, t_sums, permutation_distribution] = permutest(tempPrep',null','true',.05,10^4,'false');
            %             [clusters, p_values, t_sums, permutation_distribution] = permutest(tempPrep',null','true',.001,10^4,'false');
            for j=1:length(clusters)
                clusterMean = mean(tempPrep(:,clusters{j})); clusterMean = mean(clusterMean);
                clusterStd = std(mean(tempPrep(:,clusters{j}),2));
                clusterd(j) = clusterMean / clusterStd;
            end
            tempPrepClust = {clusters; p_values; clusterd; t_sums; permutation_distribution; tempPrep};
            %             save(['tempPrepClusters_lh' crossNames{i}],'tempPrepClust')
            
            [clusters, p_values, t_sums, permutation_distribution] = permutest(tempMov',null','true',.05,10^4,'false');
            %             [clusters, p_values, t_sums, permutation_distribution] = permutest(tempMov',null','true',.001,10^4,'false');
            for j=1:length(clusters)
                clusterMean = mean(tempMov(:,clusters{j})); clusterMean = mean(clusterMean);
                clusterStd = std(mean(tempMov(:,clusters{j}),2));
                clusterd(j) = clusterMean / clusterStd;
            end
            tempMovClust = {clusters; p_values; clusterd; t_sums; permutation_distribution; tempMov};
            %             save(['tempMovClusters_lh' crossNames{i}],'tempMovClust')
            
            [clusters, p_values, t_sums, permutation_distribution] = permutest(intPrep',null','true',.05,10^4,'false');
            %             [clusters, p_values, t_sums, permutation_distribution] = permutest(intPrep',null','true',.001,10^4,'false');
            for j=1:length(clusters)
                clusterMean = mean(intPrep(:,clusters{j})); clusterMean = mean(clusterMean);
                clusterStd = std(mean(intPrep(:,clusters{j}),2));
                clusterd(j) = clusterMean / clusterStd;
            end
            intPrepClust = {clusters; p_values; clusterd; t_sums; permutation_distribution; intPrep};
            %             save(['intPrepClusters_lh' crossNames{i}],'intPrepClust')
            
            [clusters, p_values, t_sums, permutation_distribution] = permutest(intMov',null','true',.05,10^4,'false');
            %             [clusters, p_values, t_sums, permutation_distribution] = permutest(intMov',null','true',.05,10^4,'false');
            for j=1:length(clusters)
                clusterMean = mean(intMov(:,clusters{j})); clusterMean = mean(clusterMean);
                clusterStd = std(mean(intMov(:,clusters{j}),2));
                clusterd(j) = clusterMean / clusterStd;
            end
            intMovClust = {clusters; p_values; clusterd; t_sums; permutation_distribution; intMov};
            %             save(['intMovClusters_lh' crossNames{i}],'intMovClust')
            
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGINNING OF SUIT ANALYSIS (CEREBELLUM)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    case 'SUBCORTICAL_LDA&RSA' %----------------- BEGINNING OF SUBCORTICAL -----------------%
    case 'subcortical_make_nii' %%% subcortical analysis - unzip aseg file
        %unzips aseg.nii.gz into aseg.nii and moves to subcortical folder
        
        sn = varargin{1}; %requires recon-all to be completed through freesurfer for each participant, and for 'aseg' file to be moved to relevant directory
        
        %%Anatomical:
        try
            source = fullfile(subcorticalDir, subj_name{sn}, [subj_name{sn} '_' 'aseg.nii.gz']);
            dest = fullfile(subcorticalDir, subj_name{sn});
            
            gunzip(source,dest); %unzip
            
            disp('segmentation unzip done')
        catch
            disp('No segmentation files. Most likely a typo in file name or missing destination directory!');
        end
    case 'subcortical_make_structs' %extracts each subcortical structure from aseg file, and creates respective nii files
        
        sn = varargin{1};
        
        subcortValues = cell2mat(subcortStructs(1,:));%take values from first row
        subcortStructs = subcortStructs(2,:);%and names from second
        
        for i=1:length(subcortStructs)
            
            cd([subcorticalDir, '/', subj_name{sn}])
            
            matlabbatch{1}.spm.util.imcalc.input = {[subcorticalDir, '/' subj_name{sn} '/' subj_name{sn}, '_aseg.nii,1']};
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
        
        refImageDir = [glmDir, '/', subj_name{sn}]; %directories for functional reference image (beta 1)
        subcortImageDir = [subcorticalDir, '/', subj_name{sn}]; %and for subcortical images
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for r=1:length(subcortStructs) %loop through subcortical regions and reslice them to beta image resolution
            
            matlabbatch{1}.spm.spatial.realign.write.data = {[refImageDir, '/beta_0001.nii']; [subcortImageDir, '/', subj_name{sn}, '_', subcortStructs{r}, '.nii']};
            matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
            
            spm_jobman('run',matlabbatch);
            disp([subcortStructs{r} ' resliced'])
        end
    case 'subcortical_funcmask_structs'  %masks subcortical niftis using functional grey matter mask
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for i=1:length(subcortStructs)
            funMask=fullfile(glmDir, subj_name{s},'maskbrain.nii');
            omask=fullfile(subcorticalDir, subj_name{s},['mr' subj_name{s} '_' subcortStructs{i} '.nii']); %output mask to be used in the future
            subcort = fullfile(subcorticalDir, subj_name{s},['r' subj_name{s} '_' subcortStructs{i} '.nii']);
            
            %spm_imcalc({funMask,subcort},omask,'i1 >0.01 & i2 > 0.1',{}); %recorded activity in brain (grey + white matter)
            spm_imcalc({funMask,subcort},omask,'i1 >0.01 & i2 > 0.3',{}); %recorded activity in brain (grey + white matter)
            %0.3 is the best threshold to make the resliced maps similar to
            %the anatomical segmentations
        end
    case 'subcortical_make_ROIs' %Area decoding - uses region toolbox to define R struct for prewhitening
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        R=cell(length(subcortStructs), 1);
        cd(fullfile(subcorticalDir, subj_name{s}))
        
        for i=1:length(subcortStructs)%for each subcort region
            R{i} = region('roi_image',['mr' subj_name{s} '_' subcortStructs{i} '.nii'],1,subcortStructs{i});
            R{i}.name = [subj_name{s} '_' subcortStructs{i}]; %use toolbox to define, then add name
        end%for each subcort region
        
        R = region_calcregions(R);%because we reslice, we don't need voxelspace option
        %R = region_calcregions(R, 'voxelspace', fullfile(glmDir, subj_name{s}, 'beta_0001.nii'));
        
        out = [roiSubDir '/' subj_name{s} '_subcortical_roi'];
        save(strjoin(out,''), 'R'); %save as participant file which holds all regions
    case 'subcortical_make_search' %make searchlights within each subcort struct
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        radius=16; %searchlight max radius (mm)
        numVox=160;%searchlight max voxel number
        
        saveDir = fullfile(subcorticalDir, subj_name{s});
        
        if ~isfolder(saveDir)
            mkdir(saveDir)
        end
        cd(saveDir);
        
        for r = 1:length(subcortStructs)%for subcort region
            V=spm_vol(['mr' subj_name{s} '_' subcortStructs{r} '.nii']); %if preceded by case MVA_mask
            X=spm_read_vols(V);
            [i,j,k]=ind2sub(size(X),find(X~=0));
            vox=[i j k];
            
            L = [];
            [L.LI,L.voxmin,L.voxmax,L.n]=lmva_voxelselection(vox(:,:)',vox',[radius numVox],V.mat,V.dim,[],[subcortStructs{r} '_mva160_numvox.nii']);
            L.voxel = vox;
            save ([subcortStructs{r} '_' 'volsearch160.mat'], 'L') %saves to subcortical subject directory
        end%for subcort region
    case 'subcortical_preWhiten' %pre-whiten data from subcortical ROIs ready for RSA
        
        T = []; blueBear = varargin{1}; %s=varargin{1};
        
        for s=anaSubj
            fprintf('%d.',s); fprintf('\n')
            load([glmDir, '/', subj_name{s}, '/', 'SPM.mat'],            'SPM')
            load([roiSubDir, '/', subj_name{s}, '_subcortical_roi.mat'], 'R')
            cd(fullfile(subcorticalDir, subj_name{s}))
            
            V = SPM.xY.VY;
            
            %replace SPM fnames for bluebear compatibility
            if blueBear == 1
                for i=1:length(V)
                    V(i).fname = strrep(V(i).fname,'\','/');
                    V(i).fname = strrep(V(i).fname,'Z:','/rds/projects/k/kornyshk-kornyshevalab');
                end
            else
            end
            
            for r=1:length(R)
                Y = region_getdata(V,R{r});
                
                percMove = region_getdata(spm_vol(fullfile(glmDir, subj_name{s}, 'perc_0001.nii')), R{r}); %percent signal change
                percPrep = region_getdata(spm_vol(fullfile(glmDir, subj_name{s}, 'perc_0002.nii')), R{r});
                
                [betaW,resMS,~,beta,~,~,snr] = noiseNormalizeBeta_RY(Y,SPM);
                
                S.betaW    = {betaW}; %multivariate pre-whiten
                S.betaUW   = {bsxfun(@rdivide,beta,sqrt(resMS))}; %univariate pre-whiten
                S.betaRAW  = {beta}; %Raw
                S.resMS    = {resMS}; %residual
                S.snr      = snr;
                S.percMov  = {percMove};
                S.percPrep = {percPrep};
                S.SN       = s;
                S.region   = r;
                
                T = addstruct(T,S);
                fprintf('%d.',r)
            end
            fprintf('\n');
        end
        
        save(fullfile(roiSubDir, 'preWhitened_betas.mat'),'-struct','T');
    case 'subcortical_calculate' %calculates overall distance & factorial classification acc for each subcort struct
        
        cd(roiSubDir)
        R = load('preWhitened_betas.mat');
        saveDir = fullfile(rsaDir, 'subcortical');
        
        nrruns = length(run); nCond = 8; nClassifiers = 10;
        
        %%% Distances: prepare condition and partition (run) vectors
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 0 0]; %prep
        prep=[repmat(prep,1,nrruns) runBSL];%1 x nBeta, 1 2 3 4 = prep sequences 1:4
        
        prod      =[5 0 0 0 6 0 0 0 7 0 0 0 8 0 0 0 0 0 0 0]; %prod
        prod=[repmat(prod,1,nrruns) runBSL];%1 x nBeta, 5 6 7 8 = prod sequences 1:4
        
        condVec = prep + prod; condVec = condVec'; % conditions, including no interest regressors as 0
        partVec = double(condVec > 0); %assign run number to beta
        for i=1:nCond %assign run numbers to conds, ignore no-interest betas
            partIdx = find(condVec == i);
            for j=1:length(run)
                partVec(partIdx(j)) = j;
            end
        end
        
        %%% Classifiers: prepare beta selections (prep and prod separately)
        prepBetas = condVec < 5 & condVec > 0;
        prodBetas = condVec > 4;
        
        %%% loop through prewhitened data, calculate crossnobis dissimilarities
        %pre-allocate output variables
        R.d         = NaN(length(R.SN), nCond * (nCond - 1) / 2); %pairwise distance measures
        R.Sig       = cell(length(R.SN), 1); %covariance matrix of the beta estimates across different imaging runs
        R.G         = cell(length(R.SN), 1); %second moment matrix
        R.matrix    = cell(length(R.SN), 1); % Pairwise contrast matrix
        R.acc       = NaN(length(R.SN), nClassifiers); %LDA accuracy - 1ordPrep, 2timPrep, 3intPrep, 4intPrepOneOut, 5ordProd, 6timProd, 7intProd, 8intProdOneOut
        
        for s = anaSubj%for subj
            load(fullfile(glmDir, subj_name{s}, 'SPM'), 'SPM') %load SPM design matrix for distance function
            for r = unique(R.region)'%for region
                B      = R.betaW{R.region == r & R.SN == s};
                BRAW   = R.betaRAW{R.region == r & R.SN == s};
                resMS = R.resMS{R.region == r & R.SN == s};
                
                %%%Distances
                [d, Sig] = rsa.distanceLDC(B, partVec, condVec, SPM.xX.X);
                [G,~]    = pcm_estGCrossval(B,partVec,condVec, 'X', SPM.xX.X);
                matrix   = indicatorMatrix('allpairs',1:nCond); % Pairwise contrast matrix
                
                %%%Classification
                [acc(1), acc(2), acc(3), acc(4), acc(5)]  = prepProdSimu_classify(BRAW(prepBetas, :)); %factorial classify prep sequences
                [acc(6), acc(7), acc(8), acc(9), acc(10)] = prepProdSimu_classify(BRAW(prodBetas, :)); %and prod sequences
                %acc: 1=ovrPrep, 2=ordPrep, 3=tempPrep, 4=intPrep, 5=intPrepSubtract
                %     6=ovrProd, 7=ordProd, 8=tempProd, 9=intProd, 10=intProdSubtract
                
                %%%Percent signal change
                meanPercMov  = mean(R.percMov{R.region == r & R.SN == s});
                meanPercPrep = mean(R.percPrep{R.region == r & R.SN == s});
                
                %%%Variable assignment
                R.d     (R.region == r & R.SN == s, :) = d;
                R.Sig   {R.region == r & R.SN == s}    = Sig;
                R.G     {R.region == r & R.SN == s}    = G;
                R.matrix{R.region == r & R.SN == s}    = matrix;
                R.acc   (R.region == r & R.SN == s, :) = acc;
                R.mov   (R.region == r & R.SN == s, :) = meanPercMov;
                R.prep  (R.region == r & R.SN == s, :) = meanPercPrep;
            end%for subj
        end%for region
        
        %%% Identify and extract distance values for preparation, production, and cross-phase
        phaseVector = [...
            1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...  %pairwise contrast index (1s are contrasts within preparation)
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1; ...  %^ index (1s are contrasts within production)
            0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0  ...  %^ index (1s are contrasts across phases)
            ];
        
        %%%Loop over subjs, conditions, and phases, storing distances and
        %%%accuracies and signal change as a struct to plot later
        loopCounter = 1;
        
        %%% Distance
        for i=1:length(R.betaW)%for data points
            for j=1:size(phaseVector, 1)%for phase
                Rdist.dist(loopCounter,1)   = mean(R.d(i, phaseVector(j,:) == 1));
                Rdist.phase(loopCounter,1)  = j;
                Rdist.cond(loopCounter,1)   = 1;
                Rdist.SN(loopCounter,1)     = R.SN(i);
                Rdist.region(loopCounter,1) = R.region(i);
                
                loopCounter = loopCounter + 1;
            end%for phase
        end%for data points
        
        loopCounter = 1;
        
        %%% Decoding
        for i=1:length(R.betaW)%for data points
            for j=1:size(R.acc, 2)%for classifiers
                Racc.acc(loopCounter,1)   = R.acc(i, j);
                if j < 6%if prep classifier
                    Racc.phase(loopCounter,1) = 1;
                    Racc.cond(loopCounter,1)  = j;
                elseif j > 5%if prod classifier
                    Racc.phase(loopCounter,1) = 2;
                    Racc.cond(loopCounter,1)  = j - 5; %keeps conds 1:5
                end
                Racc.SN(loopCounter,1)    = R.SN(i);
                Racc.region(loopCounter,1) = R.region(i);
                
                loopCounter = loopCounter + 1;
            end%for classifiers
        end%for subj
        
        R;
        loopCounter = 1;
        
        %%% Distance
        for i=1:length(R.betaW)%for data points
            for j=1:2%for phase
                
                if j == 1
                    Rperc.perc(loopCounter,1) = mean(R.percPrep{i});
                elseif j == 2
                    Rperc.perc(loopCounter,1) = mean(R.percMov{i});
                end
                Rperc.phase(loopCounter,1)    = j;
                Rperc.cond(loopCounter,1)     = 1;
                Rperc.SN(loopCounter,1)       = R.SN(i);
                Rperc.region(loopCounter,1)   = R.region(i);
                
                loopCounter = loopCounter + 1;
            end%for phase
        end%for data points
        
        if ~isfolder(saveDir)
            mkdir(saveDir)
        end
        
        save(fullfile(saveDir, 'subRoiDistances.mat'), 'R', 'Rdist', 'Racc', 'Rperc')
    case 'subcortical_mds_distances' %cross-phase MDS distance when PC1 is ignored - saves to subRoiDistances_withmds.mat
        
        %%%Load
        dataFile = fullfile(rsaDir, 'subcortical', 'subRoiDistances.mat');
        load(dataFile, 'R', 'Rdist', 'Racc', 'Rperc')%load it
        
        %Set eucScaling and pcs variables to same size as Rdist
        Rdist.eucScaling = nan(length(Rdist.SN),1);
        Rdist.pcs        = cell(size(Rdist.SN));
        Rdist.eigVals    = cell(size(Rdist.SN));
        
        %%%Set save directory - we re-save the distances to include
        %%%cross-phase MDS
        saveDir = fullfile(rsaDir, 'subcortical');
        
        %%% RDMs and multi-dimensional scaling plots (for visualisation)
        labels = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        
        %Extract data
        for s=1:length(R.G)%for subj & region
            G(:,:,s)= R.G{s}; %extract representational dissimilarity matrix from cells into 3D matrix
        end%for subj * region
        
        for i=unique(R.region)'%for region
            COORD = [];
            GRegion = G(:,:,R.region == i); %extract region variance/covariance matrices
            
            for j=1:length(GRegion) %get MDS for each subj in the region
                [COORD(:,:,j),eigVals(:,j)] = pcm_classicalMDS(GRegion(:,:,j));
            end
            
            PCs = COORD(:,2:3,:); %PCs 2 and 3
            
            %Euclidean distance between phases within sequences for PC2
            %and PC3
            for j=1:length(PCs)%for subj
                for k=1:4%for within-sequence prep vs prod
                    eucPcDist(k) = pdist([PCs(k,:,j)' PCs(k+4,:,j)'], 'euclidean');
                end%for within-sequence prep vs prod
                
                %pairwise distance measures k(k-1) / 2
                %within prep
                prepDists(1) = pdist([PCs(1,:,j)', PCs(2,:,j)'],'euclidean');
                prepDists(2) = pdist([PCs(1,:,j)', PCs(3,:,j)'],'euclidean');
                prepDists(3) = pdist([PCs(1,:,j)', PCs(4,:,j)'],'euclidean');
                prepDists(4) = pdist([PCs(2,:,j)', PCs(3,:,j)'],'euclidean');
                prepDists(5) = pdist([PCs(2,:,j)', PCs(4,:,j)'],'euclidean');
                prepDists(6) = pdist([PCs(3,:,j)', PCs(4,:,j)'],'euclidean');
                %within prod
                prodDists(1) = pdist([PCs(5,:,j)', PCs(6,:,j)'],'euclidean');
                prodDists(2) = pdist([PCs(5,:,j)', PCs(7,:,j)'],'euclidean');
                prodDists(3) = pdist([PCs(5,:,j)', PCs(8,:,j)'],'euclidean');
                prodDists(4) = pdist([PCs(6,:,j)', PCs(7,:,j)'],'euclidean');
                prodDists(5) = pdist([PCs(6,:,j)', PCs(8,:,j)'],'euclidean');
                prodDists(6) = pdist([PCs(7,:,j)', PCs(8,:,j)'],'euclidean');
                %%%FIND A BETTER WAY TO CODE THIS ^
                
                prepSize = mean(prepDists);
                prodSize = mean(prodDists);
                avgSize  = mean([prepSize, prodSize]);
                
                %Collect variables to later add to RDist
                Rdist.dist(end+1,1)        = mean(eucPcDist);
                Rdist.eucScaling(end+1,1)  = avgSize;
                Rdist.phase(end+1,1)       = 0;
                Rdist.cond(end+1,1)        = 2;
                Rdist.SN(end+1,1)          = anaSubj(j);
                Rdist.region(end+1,1)      = i;
                Rdist.pcs{end+1,1}         = COORD(:,:,j);
                Rdist.eigVals{end+1,1}     = eigVals(:,j);
            end%for subj
        end%for region
        
        save(fullfile(saveDir, 'subRoiDistances_withmds.mat'), 'R', 'Rdist', 'Racc', 'Rperc')
        
    case 'subcortical_run_search_RSA'
        
        s=varargin{1}; blueBear=varargin{2};
        subcortStructs = subcortStructs(2,:);%names from second row
        nrruns = length(run);
        
        cd(fullfile(subcorticalDir, subj_name{s}))
        
        for r=1:length(subcortStructs)
            
            %%% Searchlight file
            load([subcortStructs{r} '_volsearch160.mat'], 'L');
            
            %%% SPM file
            spmDir = fullfile(glmDir, subj_name{s});
            load(fullfile(spmDir, 'SPM'), 'SPM');
            
            %replace fnames for bluebear compatibility
            if blueBear == 1
                for i=1:length(SPM.xY.VY)
                    SPM.xY.VY(i).fname = strrep(SPM.xY.VY(i).fname,'\','/');
                    SPM.xY.VY(i).fname = strrep(SPM.xY.VY(i).fname,'Z:','/rds/projects/k/kornyshk-kornyshevalab');
                end
            else
            end
            
            %%% prepare condition and partition (run) vectors
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 0 0]; %prep
            prep=[repmat(prep,1,nrruns) runBSL];%1 x nBeta, 1 2 3 4 = prep sequences 1:4
            
            prod      =[5 0 0 0 6 0 0 0 7 0 0 0 8 0 0 0 0 0 0 0]; %prod
            prod=[repmat(prod,1,nrruns) runBSL];%1 x nBeta, 5 6 7 8 = prod sequences 1:4
            
            condVec = prep + prod; condVec = condVec'; % conditions, including no interest regressors as 0
            
            %%% Run searchlight function on subcortical structure
            rsa.runSearchlightLDC_RY(L, SPM, 'spmDir', spmDir, 'conditionVec', condVec, ...
                'analysisName', [subj_name{s}, '_' subcortStructs{r} '_RSA_All'], 'outDir', fullfile(subcorticalDir, subj_name{s}))
        end
    case 'subcortical_run_overallMov_LDA'  % Conduct the classification analysis 4 sequences
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for r=1:length(subcortStructs)
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(subcorticalDir, subj_name{s}, [subcortStructs{r} '_volsearch160.mat']), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Comb_160_Mov.nii'])};
            
            
            % Generate column indices for Cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
                train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
            end
            
            [~,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
        end
    case 'subcortical_run_overallPrep_LDA' % Conduct the classification analysis 4 sequences
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for r=1:length(subcortStructs)
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(subcorticalDir, subj_name{s}, [subcortStructs{r} '_volsearch160.mat']), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
            nrruns=length(SPM.nscan);
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Comb_160_Prep.nii'])};
            
            
            % Generate column indices for Cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
                train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
            end
            
            [~,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
        end
    case 'subcortical_run_spatMov_LDA'
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for r=1:length(subcortStructs)
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(subcorticalDir, subj_name{s}, [subcortStructs{r} '_volsearch160.mat']), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
            nrruns=length(SPM.nscan);
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Spat_160_Mov.nii'])};
            
            
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
            end
            
            [~,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
        end
    case 'subcortical_run_spatPrep_LDA'
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for r=1:length(subcortStructs)
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(subcorticalDir, subj_name{s}, [subcortStructs{r} '_volsearch160.mat']), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 1 2 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Spat_160_Prep.nii'])};
            
            
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
                
            end
            
            [~,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
        end
    case 'subcortical_run_tempMov_LDA'
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        
        for r=1:length(subcortStructs)
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(subcorticalDir, subj_name{s}, [subcortStructs{r} '_volsearch160.mat']), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Temp_160_Mov.nii'])};
            
            
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
                
            end
            
            [~,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
            
        end
    case 'subcortical_run_tempPrep_LDA'
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for r=1:length(subcortStructs)
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(subcorticalDir, subj_name{s}, [subcortStructs{r} '_volsearch160.mat']), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat([1 2 1 2],1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Temp_160_Prep.nii'])};
            
            
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
                
            end
            
            [~,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
        end
    case 'subcortical_run_intMov_LDA'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for r=1:length(subcortStructs)
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(subcorticalDir, subj_name{s}, [subcortStructs{r} '_volsearch160.mat']), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
            nrruns=length(SPM.nscan);
            
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
            prod=[repmat(prod,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Int_160_Mov.nii'])};
            
            
            % Generate column indices for Cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
                train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
            end
            
            [~,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
        end
    case 'subcortical_run_intPrep_LDA'   %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for r=1:length(subcortStructs)
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(subcorticalDir, subj_name{s}, [subcortStructs{r} '_volsearch160.mat']), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
            nrruns=length(SPM.nscan);
            
            runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
            prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
            prep=[repmat(prep,1,nrruns) runBSL];
            
            c=repmat(1:4,1,nrruns); % extract conditions (leave out errors!)
            run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            out = {fullfile(subcorticalDir, subj_name{s}, [subj_name{s}, '_' subcortStructs{r} '_accuracy_Int_160_Prep.nii'])};
            
            
            % Generate column indices for Cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
                train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
            end
            
            [~,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
        end
    case 'subcortical_zValue_LDA'
        
        s=varargin{1};
        cd(fullfile(subcorticalDir,subj_name{s}));
        subcortStructs = subcortStructs(2,:);%names from second row
        
        %zValue the order and timing maps
        for r=1:length(subcortStructs)
            
            takeOneOutIter=2;
            numTests=6;
            numCat=2;
            mu=1/numCat; %mu=0.5;
            N=numTests*numCat*takeOneOutIter;
            sigma=sqrt(mu*(1-mu)*1/N);
            
            images= {...
                '_accuracy_Spat_160_Mov','_accuracy_Spat_160_Prep',...
                '_accuracy_Temp_160_Mov','_accuracy_Temp_160_Prep',...
                };
            outimages={...
                '_zacc_Spat_160_Mov','_zacc_Spat_160_Prep',...
                '_zacc_Temp_160_Mov','_zacc_Temp_160_Prep',...
                };
            
            for j=1:numel(images)
                input_image= fullfile([subj_name{s} '_' subcortStructs{r} images{j} '.nii']);
                output_image= fullfile([subj_name{s} '_' subcortStructs{r} outimages{j} '.nii']);
                spmj_imcalc_mtx(input_image, output_image,...
                    sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
            end
        end
        
        %zValue the integrated maps
        for r=1:length(subcortStructs)
            
            numTests=6;
            numCat=4;
            mu=1/numCat; %mu=0.25;
            N=numTests*numCat;
            sigma=sqrt(mu*(1-mu)*1/N);
            
            images= {...
                '_accuracy_Comb_160_Mov', '_accuracy_Comb_160_Prep',...
                '_accuracy_Int_160_Mov', '_accuracy_Int_160_Prep',...
                };
            outimages={...
                '_zacc_Comb_160_Mov', '_zacc_Comb_160_Prep',...
                '_zacc_Comb_160_Mov', '_zacc_Comb_160_Prep',...
                };
            
            for j=1:numel(images)
                input_image= fullfile([subj_name{s} '_' subcortStructs{r} images{j} '.nii']);
                output_image= fullfile([subj_name{s} '_' subcortStructs{r} outimages{j} '.nii']);
                spmj_imcalc_mtx(input_image, output_image,...
                    sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
            end
        end
    case 'subcortical_wholebrain_percent_signal' %generates whole brain percent signal change maps (perc_0001.nii and perc_0002.nii)
        %%%READ ME%%%
        %The below functions do not provide accurate PSC maps (on
        %blueBear). The correct maps, using the same code, were generated
        %on the office Mac used to perform the surface analysis. The
        %perc_000* files that are used for subsequent analysis in this
        %script are taken from there. Do not overwrite them, or you will
        %have to retrieve them again. I intend to troubleshoot why the same
        %code does not produce the same results, but for now do not
        %overwrite. RY
        
        
        %%%Calculate the same way as JNeurosci paper (https://www.jneurosci.org/content/43/10/1742/):
        %sn=varargin{1};
        %for s=sn
        %    cd(fullfile(glmDir, subj_name{s}));
        %    load('SPM.mat', 'SPM');
        %    Vrunmean=SPM.Vbeta(SPM.xX.iB);
        %    P={Vrunmean.fname};
        %    spmj_imcalc_mtx(P,'meanRestEPI.nii','mean(X)');
        %    spm_imcalc({'meanRestEPI.nii','con_0001.nii'},'perc_0001.nii','i2./i1*100');
        %    spm_imcalc({'meanRestEPI.nii','con_0002.nii'},'perc_0002.nii','i2./i1*100');
        %end
        
        %%%Calculate the same way as Berlot:
        %sn=varargin{1};
        %for s=sn
        %   cd(fullfile(glmDir, subj_name{s}));
        %   load('SPM.mat', 'SPM');
        %   Vrunmean=SPM.Vbeta(SPM.xX.iB); %rest betas
        %   P={Vrunmean.fname};
        %   spmj_imcalc_mtx(P,'meanRestEPI.nii','mean(X)'); %average rest betas
        %
        %   X=(SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
        %   h=median(max(X));               % Height of response;
        %
        %   spm_imcalc_ui({'meanRestEPI.nii','con_0001.nii'}, 'psc_0001.nii', ['100 .* ' num2str(h) ' .* i2 ./ i1']); %percent signal change - movement
        %   spm_imcalc_ui({'meanRestEPI.nii','con_0002.nii'}, 'psc_0002.nii', ['100 .* ' num2str(h) ' .* i2 ./ i1']); %percent signal change - rest
        %end
    case 'subcortical_segment_percent_signal' %separates whole brain perc map into subcortical structures using masks
        
        sn = varargin{1};
        cd([subcorticalDir, '/', subj_name{sn}])
        
        subcortValues = cell2mat(subcortStructs(1,:));%take values from first row
        subcortStructs = subcortStructs(2,:);%and names from second
        fileName = {'perc_0001', 'perc_0002'};
        
        for i=1:length(subcortStructs)
            for j=1:length(fileName)
                Vi(1) = spm_vol(fullfile(glmDir, subj_name{sn}, [fileName{j} '.nii'])); %percent signal change map
                Vi(2) = spm_vol(fullfile(subcorticalDir, subj_name{sn}, ['/mr' subj_name{sn} '_' subcortStructs{i} '.nii'])); %subcortical mask
                
                Vo = Vi(1);
                Vo.fname = [subj_name{sn}, '_', subcortStructs{i} '_' fileName{j} '.nii'];
                Vo = rmfield(Vo, 'pinfo');
                
                f = 'i1 .* (i2 > 0)';
                flags.mask   = 1; %0s treated as NaNs
                flags.interp = 2; %gives us good interpolation to whole brain maps
                flags.dtype  = 'float';
                
                
                Vo = spm_imcalc(Vi, Vo, f, flags);
            end
        end
    case 'subcortical_segment_contrasts' %separates whole brain con maps into subcortical structures using masks
        
        sn = varargin{1};
        cd([subcorticalDir, '/', subj_name{sn}])
        
        subcortValues = cell2mat(subcortStructs(1,:));%take values from first row
        subcortStructs = subcortStructs(2,:);%and names from second
        fileName = {'con_0001', 'con_0002'};
        
        for i=1:length(subcortStructs)
            for j=1:length(fileName)
                Vi(1) = spm_vol(fullfile(glmDir, subj_name{sn}, [fileName{j} '.nii'])); %percent signal change map
                Vi(2) = spm_vol(fullfile(subcorticalDir, subj_name{sn}, ['/mr' subj_name{sn} '_' subcortStructs{i} '.nii'])); %subcortical mask
                
                Vo = Vi(1);
                Vo.fname = [subj_name{sn}, '_', subcortStructs{i} '_' fileName{j} '.nii'];
                Vo = rmfield(Vo, 'pinfo');
                
                f = 'i1 .* (i2 > 0)';
                flags.mask   = 1; %0s treated as NaNs
                flags.interp = 2; %gives us good interpolation to whole brain maps
                flags.dtype  = 'float';
                
                
                Vo = spm_imcalc(Vi, Vo, f, flags);
            end
        end
    case 'subcortical_smooth' %Smoothing in subject space
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for r=1:length(subcortStructs)
            
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_RSA_All_LDC.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_RSA_All_sLDC.nii']);
            spm_smooth(comb,scomb,[2 2 2]); %smooth with 2mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Comb_160_Mov.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Comb_160_Mov.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Comb_160_Prep.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Comb_160_Prep.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Spat_160_Mov.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Spat_160_Mov.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Spat_160_Prep.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Spat_160_Prep.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Temp_160_Mov.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Temp_160_Mov.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Temp_160_Prep.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Temp_160_Prep.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Int_160_Mov.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Int_160_Mov.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_zacc_Int_160_Prep.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_szacc_Int_160_Prep.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_perc_0001.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_sperc_0001.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_perc_0002.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_sperc_0002.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_con_0001.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_scon_0001.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            s=varargin{1};
            comb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_con_0002.nii']); %%MVPA smoother
            scomb=fullfile(subcorticalDir, subj_name{s},[subj_name{s} '_' subcortStructs{r} '_scon_0002.nii']);
            spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
            
            %smooth other images here as required
        end
    case 'subcortical_calc_dissimilarity_maps' %extract overall distances from RSA_ALL_sLDC.nii
        %each volume of the searchlight corresponds to a pairwise
        %dissimilarity measure between sequences.
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        cd(fullfile(subcorticalDir, subj_name{s}))
        
        %%% Identify and extract values for overall, order, and timing
        %%% within preparation, production, and cross-phase
        prepCols       = [1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';   %^ index (1s are contrasts within preparation)
        prodCols       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1]';   %^ index (1s are contrasts within production)
        crossPhaseCols = [0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0]';   %^ index (1s are contrasts across phases within sequences)
        
        conds = {...
            'overall_prep', 'overall_prod', 'overall_cross'; ...
            prepCols == 1,  prodCols == 1,  crossPhaseCols == 1, ...
            };
        
        for r=1:length(subcortStructs)
            for j=1:length(conds)
                vol = spm_vol([subj_name{s}, '_' subcortStructs{r} '_RSA_All_sLDC.nii']);
                
                Vi = vol(conds{2,j}, :);
                Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
                Vo.fname = [subj_name{s} '_' subcortStructs{r} '_LDC_' conds{1,j} '.nii'];
                Vo.n = [1 1];
                express = 'mean(X)';
                flags.dmtx = 1;
                
                spm_imcalc(Vi, Vo, express, flags)
            end
        end
    case 'subcortical_normalise_RSA' %normalisation into MNI space - RSA
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        if ~isfolder(fullfile(subcorticalGroupRSADir,'data'))
            mkdir(fullfile(subcorticalGroupRSADir,'data')); % folder for each contrast
        end
        
        %RSA distance maps
        images = {...
            'overall_prep.nii',    'overall_prod.nii',    'overall_cross.nii', ...
            };
        
        for r=1:length(subcortStructs)%for subcort region
            for i=1:length(images)%for decoder
                defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
                [~,name,~]=spm_fileparts(images{i});
                sn_images{1}= fullfile(subcorticalDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_LDC_' images{i}]);
                out_images{1}= fullfile(subcorticalGroupRSADir,[name '_LDC_' subcortStructs{r} '_' subj_name{s} '.nii']);
                spmj_normalization_write(defor, sn_images,'outimages',out_images); %Trilinear interpolation
            end%for decoder
        end%for subcort region
    case 'subcortical_normalise_LDA'
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        if ~isfolder(fullfile(subcorticalGroupDir,'data'))
            mkdir(fullfile(subcorticalGroupDir,'data')); % folder for each contrast
        end
        
        %MVPA accuracy maps
        images= {...
            'szacc_Comb_160_Mov.nii','szacc_Comb_160_Prep.nii',...
            'szacc_Spat_160_Mov.nii','szacc_Spat_160_Prep.nii',...
            'szacc_Temp_160_Mov.nii','szacc_Temp_160_Prep.nii',...
            'szacc_Int_160_Mov.nii','szacc_Int_160_Prep.nii',...
            }; % add other images as required
        for r=1:length(subcortStructs)%for subcort region
            for i=1:length(images)%for decoder
                defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
                [~,name,~]=spm_fileparts(images{i});
                sn_images{1}= fullfile(subcorticalDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_' images{i}]);
                out_images{1}= fullfile(subcorticalGroupDir,[name '_' subcortStructs{r} '_' subj_name{s} '.nii']);
                spmj_normalization_write(defor, sn_images,'outimages',out_images); %Trilinear interpolation
            end%for decoder
        end%for subcort region
    case 'subcortical_normalise_conperc'
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        if ~isfolder(fullfile(subcorticalGroupDir,'data'))
            mkdir(fullfile(subcorticalGroupDir,'data')); % folder for each contrast
        end
        
        %MVPA accuracy maps
        images= {...
            'sperc_0001.nii','sperc_0002.nii', ...
            'scon_0001.nii', 'scon_0002.nii', ...
            }; % add other images as required
        for r=1:length(subcortStructs)%for subcort region
            for i=1:length(images)%for decoder
                defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
                [~,name,~]=spm_fileparts(images{i});
                sn_images{1}= fullfile(subcorticalDir,subj_name{s},[subj_name{s} '_' subcortStructs{r} '_' images{i}]);
                out_images{1}= fullfile(subcorticalGroupDir,[name '_' subcortStructs{r} '_' subj_name{s} '.nii']);
                spmj_normalization_write(defor, sn_images,'outimages',out_images); %Trilinear interpolation
            end%for decoder
        end%for subcort region
    case 'subcortical_group_avg_RSA' %group average RSA maps
        
        subcortStructs = subcortStructs(2,:);%names from second row
        
        if ~isfolder(fullfile(subcorticalGroupRSADir, 'average'))
            mkdir(fullfile(subcorticalGroupRSADir, 'average'))
        end
        cd(fullfile(subcorticalGroupRSADir, 'average'))
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        for r=1:length(subcortStructs)%for subcort region
            for i=1:length(conds)
                loopCount = 1;
                for s = anaSubj
                    Vi(loopCount) = spm_vol(fullfile(subcorticalGroupRSADir, [conds{i} '_LDC_' subcortStructs{r} '_' subj_name{s} '.nii']));
                    loopCount = loopCount + 1;
                end
                Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
                Vo.fname = ['avg_' conds{i} '_' subcortStructs{r} '_LDC.nii'];
                Vo.n = [1 1];
                express = 'mean(X)';
                flags.dmtx = 1;
                
                spm_imcalc(Vi, Vo, express, flags)
            end
        end
    case 'subcortical_group_avg_LDA' %group average LDA maps
        
        subcortStructs = subcortStructs(2,:);%names from second row
        
        if ~isfolder(fullfile(subcorticalGroupDir, 'average'))
            mkdir(fullfile(subcorticalGroupDir, 'average'))
        end
        cd(fullfile(subcorticalGroupDir, 'average'))
        
        images = {...
            'szacc_Comb_160_Mov';'szacc_Comb_160_Prep';...
            'szacc_Spat_160_Mov';'szacc_Spat_160_Prep';...
            'szacc_Temp_160_Mov';'szacc_Temp_160_Prep';...
            'szacc_Int_160_Mov'; 'szacc_Int_160_Prep'...
            };
        
        conds = {...
            'comb_mov', 'comb_prep',...
            'spat_mov', 'spat_prep',...
            'temp_mov', 'temp_prep', ...
            'int_mov', 'int_prep', ...
            };
        for r=1:length(subcortStructs)%for subcort region
            for i=1:length(images)
                loopCount = 1;
                for s = anaSubj
                    Vi(loopCount) = spm_vol(fullfile(subcorticalGroupDir, [images{i} '_' subcortStructs{r} '_' subj_name{s} '.nii']));
                    loopCount = loopCount + 1;
                end
                Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
                Vo.fname = ['avg_' conds{i} '_' subcortStructs{r} '_LDA.nii'];
                Vo.n = [1 1];
                express = 'mean(X)';
                flags.dmtx = 1;
                
                spm_imcalc(Vi, Vo, express, flags)
            end
        end
    case 'subcortical_group_avg_perc' %group average LDA maps
        
        subcortStructs = subcortStructs(2,:);%names from second row
        
        if ~isfolder(fullfile(subcorticalGroupDir, 'average'))
            mkdir(fullfile(subcorticalGroupDir, 'average'))
        end
        cd(fullfile(subcorticalGroupDir, 'average'))
        
        images = {...
            'sperc_0001';'sperc_0002';...
            };
        
        conds = {...
            'perc_mov', 'perc_prep',...
            };
        for r=1:length(subcortStructs)%for subcort region
            for i=1:length(images)
                loopCount = 1;
                for s = anaSubj
                    Vi(loopCount) = spm_vol(fullfile(subcorticalGroupDir, [images{i} '_' subcortStructs{r} '_' subj_name{s} '.nii']));
                    loopCount = loopCount + 1;
                end
                Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
                Vo.fname = ['avg_' conds{i} '_' subcortStructs{r} '_perc.nii'];
                Vo.n = [1 1];
                express = 'mean(X)';
                flags.dmtx = 1;
                
                spm_imcalc(Vi, Vo, express, flags)
            end
        end
    case 'subcortical_group_randomeffects_RSA' %random effects - RSA
        
        subcortStructs = subcortStructs(2,:);%names from second row
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        
        dataDir = strcat('RSA_', conds);
        
        contrastN = length(dataDir);
        images = repmat(conds', 1, length(anaSubj));
        subNii = repmat (subj_name(anaSubj), length(dataDir), 1);
        
        for r=1:length(subcortStructs)%for subcort regions
            fileName = strcat (images, '_LDC_', subcortStructs{r}, '_', subNii, '.nii');  %%Concatenate contrast files and subject names
            subcortDataDir = strcat(dataDir, '_', subcortStructs{r});
            
            for i=1:contrastN%for contrasts
                glmscndDir = fullfile(subcorticalGroupRSADir, subcortDataDir(i));
                matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;  %Adjust directory
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (subcorticalGroupRSADir, fileName(i,:))';  %%Select files from matrix
                matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
                matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
                matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
                matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
                matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
                matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
                matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
                matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
                
                spm_jobman('run',matlabbatch);
            end%for contrasts
        end%for subcort regions
    case 'subcortical_group_randomeffects_LDA' %random effects - LDA
        
        subcortStructs = subcortStructs(2,:);%names from second row
        images = {...
            'szacc_Comb_160_Mov';'szacc_Comb_160_Prep';...
            'szacc_Spat_160_Mov';'szacc_Spat_160_Prep';...
            'szacc_Temp_160_Mov';'szacc_Temp_160_Prep';...
            'szacc_Int_160_Mov'; 'szacc_Int_160_Prep'...
            };
        
        conds = {...
            'comb_mov', 'comb_prep',...
            'spat_mov', 'spat_prep',...
            'temp_mov', 'temp_prep', ...
            'int_mov', 'int_prep', ...
            };
        
        dataDir = strcat('MVA_', conds);
        
        contrastN = length(dataDir);
        images    = repmat(images, 1, length(anaSubj));
        subNii    = repmat (subj_name(anaSubj), length(dataDir), 1);
        
        for r=1:length(subcortStructs)%for subcort regions
            fileName       = strcat (images, '_', subcortStructs{r}, '_', subNii, '.nii');  %%Concatenate contrast files and subject names
            subcortDataDir = strcat(dataDir, '_', subcortStructs{r});
            
            for i=1:contrastN%for contrasts
                glmscndDir = fullfile(subcorticalGroupDir, subcortDataDir(i));
                matlabbatch{1}.spm.stats.factorial_design.dir                    = glmscndDir;  %Adjust directory
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans           = fullfile (subcorticalGroupDir, fileName(i,:))';  %%Select files from matrix
                matlabbatch{1}.spm.stats.factorial_design.cov                    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
                matlabbatch{1}.spm.stats.factorial_design.multi_cov              = struct('files', {}, 'iCFI', {}, 'iCC', {});
                matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
                matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;
                matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};
                matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
                matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
                matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;
                
                spm_jobman('run',matlabbatch);
            end%for contrasts
        end%for subcort regions
    case 'subcortical_group_randomeffects_perc'%random effects - %SC
        
        subcortStructs = subcortStructs(2,:);%names from second row
        images = {...
            'sperc_0001';'sperc_0002';...
            };
        
        conds = {...
            'perc_mov', 'perc_prep',...
            };
        
        dataDir = conds;
        
        contrastN = length(dataDir);
        images    = repmat(images, 1, length(anaSubj));
        subNii    = repmat (subj_name(anaSubj), length(dataDir), 1);
        
        for r=1:length(subcortStructs)%for subcort regions
            fileName       = strcat (images, '_', subcortStructs{r}, '_', subNii, '.nii');  %%Concatenate contrast files and subject names
            subcortDataDir = strcat(dataDir, '_', subcortStructs{r});
            
            for i=1:contrastN%for contrasts
                glmscndDir = fullfile(subcorticalGroupDir, subcortDataDir(i));
                matlabbatch{1}.spm.stats.factorial_design.dir                    = glmscndDir;  %Adjust directory
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans           = fullfile (subcorticalGroupDir, fileName(i,:))';  %%Select files from matrix
                matlabbatch{1}.spm.stats.factorial_design.cov                    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
                matlabbatch{1}.spm.stats.factorial_design.multi_cov              = struct('files', {}, 'iCFI', {}, 'iCC', {});
                matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
                matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;
                matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};
                matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
                matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
                matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;
                
                spm_jobman('run',matlabbatch);
            end%for contrasts
        end%for subcort regions
    case 'subcortical_group_randomeffects_con' %random effects - con
        
        subcortStructs = subcortStructs(2,:);%names from second row
        
        dataDir = {'Mov', 'Prep'}; %%Save folders for each contrast
        images = {'scon_0001';'scon_0002'};
        
        contrastN = length(dataDir);
        images    = repmat(images, 1, length(anaSubj));
        subNii    = repmat (subj_name(anaSubj), length(dataDir), 1);
        
        for r=1:length(subcortStructs)%for subcort regions
            fileName       = strcat (images, '_', subcortStructs{r}, '_', subNii, '.nii');  %%Concatenate contrast files and subject names
            subcortDataDir = strcat('data/', dataDir, '_', subcortStructs{r});
            
            for i=1:contrastN%for contrasts
                glmscndDir = fullfile(subcorticalGroupDir, subcortDataDir(i));
                matlabbatch{1}.spm.stats.factorial_design.dir                    = glmscndDir;  %Adjust directory
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans           = fullfile (subcorticalGroupDir, fileName(i,:))';  %%Select files from matrix
                matlabbatch{1}.spm.stats.factorial_design.cov                    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
                matlabbatch{1}.spm.stats.factorial_design.multi_cov              = struct('files', {}, 'iCFI', {}, 'iCC', {});
                matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
                matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;
                matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};
                matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
                matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
                matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;
                
                spm_jobman('run',matlabbatch);
            end%for contrasts
        end%for subcort regions
    case 'subcortical_group_estimate_RSA' %estimate RSA group
        
        subcortStructs = subcortStructs(2,:);%names from second row
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        
        contrastN = length(conds);
        for r=1:length(subcortStructs)%for subcort regions
            for i=1:contrastN%for contrasts
                matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile (subcorticalGroupRSADir, ['RSA_' conds{i} '_' subcortStructs{r}], 'SPM.mat');  %Adjust directory
                matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
                matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                
                spm_jobman('run',matlabbatch);
            end%for contrasts
        end%for subcort regions
    case 'subcortical_group_estimate_LDA' %estimate LDA group
        
        subcortStructs = subcortStructs(2,:);%names from second row
        conds = {...
            'comb_prep',   'comb_mov',...
            'int_prep',    'int_mov',...
            'spat_prep',   'spat_mov',...
            'temp_prep',   'temp_mov',...
            };
        
        contrastN = length(conds);
        for r=1:length(subcortStructs)%for subcort regions
            for i=1:contrastN%for contrasts
                matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile (subcorticalGroupDir, ['MVA_' conds{i} '_' subcortStructs{r}], 'SPM.mat');  %Adjust directory
                matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
                matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                
                spm_jobman('run',matlabbatch);
            end%for contrasts
        end%for subcort regions
    case 'subcortical_group_estimate_con' %estimate LDA group
        
        subcortStructs = subcortStructs(2,:);%names from second row
        conds = {...
            'Mov',    'Prep',...
            };
        
        contrastN = length(conds);
        for r=1:length(subcortStructs)%for subcort regions
            for i=1:contrastN%for contrasts
                matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile (subcorticalGroupDir, 'data', [conds{i} '_' subcortStructs{r}], 'SPM.mat');  %Adjust directory
                matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
                matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                
                spm_jobman('run',matlabbatch);
            end%for contrasts
        end%for subcort regions
        
        
    case 'subcortical_normalise_anat_masks' %Generate underlay to display results: normalise structure masks...
        
        subcortStructs = subcortStructs(2,:);%names from second row
        sn = varargin{1};
        inDir  = fullfile(subcorticalDir,  subj_name{sn});
        outDir = fullfile(subcorticalGroupDir, 'atlas');
        
        if ~isfolder(outDir)
            mkdir(outDir)
        end
        
        %%%Assemble file names
        defor = fullfile(anatDir, subj_name{sn}, [subj_name{sn}, '_anatomical_seg_sn.mat']); %defor
        
        %images    = cell(length(subcortStructs), 1); %preallocate
        %outimages = cell(length(subcortStructs), 1); %preallocate
        for i=1:length(subcortStructs)
            images{i}    = fullfile(inDir,  [subj_name{sn} '_' subcortStructs{i} '.nii']); %segmented subcort regions
            outimages{i} = fullfile(outDir, [subj_name{sn} '_' subcortStructs{i} '.nii']); %segmented subcort regions
        end
        
        %%%Run normalisation
        spmj_normalization_write(defor, images', 'outimages', outimages') %normalise images
    case 'subcortical_average_anat_masks'   %then average to use to display results
        
        subcortStructs = subcortStructs(2,:);%names from second row
        fileDir = fullfile(subcorticalGroupDir, 'atlas');
        cd(fileDir)
        
        %%%Average normalised structure maps and save
        for i=1:length(subcortStructs)
            in  = dir(['*_' subcortStructs{i} '.nii']); in = {in.name};
            
            out = spm_vol(in{1}); out = rmfield(out, 'pinfo');
            out.fname = [subcortStructs{i} '.nii'];
            
            spmj_imcalc_mtx(in,out,'mean(X)');
        end
        
        %%%Sum maps to produce 'all.nii' which includes all structures
        in  = cell(length(subcortStructs), 1);
        
        for i=1:length(subcortStructs)
            in{i} = [subcortStructs{i} '.nii'];
        end
        
        out       = spm_vol(in{1});  out = rmfield(out, 'pinfo');
        out.fname = 'all.nii';
        
        spmj_imcalc_mtx(in,out,'sum(X)');
    case 'subcortical_plot' %plot area distance and decoding results
        %%%Also requires case subcortical_mds_distances to run
        subcortStructs = subcortStructs(2, :);%names from second row
        subcortStructs = strrep(subcortStructs, '_', ' '); %replace _ with space
        
        %Text
        titleFontSize = 12;
        fontSize      = 12;
        
        %Colour
        decodeBRG = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880], [0.4660 0.6740 0.1880]}; %blue red green
        %decodeBRG = {[1, 0.65, 0], [0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880], [0.4660 0.6740 0.1880]}; %with yellow for overall
        rsaGB     = {[0.8 0.8 0.8], [0 0 0]};
        
        %Style
        spl = [zeros(4,1); ones(4,1)];
        label = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        ms = 12;
        ls = 20;
        lw = 3;
        
        %y axis scales
        pscY     = [-0.5,  1.0];
        rsaY     = [-0.003, 0.0045];
        ldaY     = [-4.1,  4.5];
        crossY   = [0,      0.11];
        mdsDistY = [0       0.07];
        eigY     = [-0.0003 0.29];
        
        %%%Load
        dataFile = fullfile(rsaDir, 'subcortical', 'subRoiDistances_withmds.mat');
        load(dataFile, 'R', 'Rdist', 'Racc', 'Rperc')%load it
        
        %%%Collate data
        Tperc = tapply(Rperc,{'SN', 'region', 'phase'},{'perc', 'mean', 'name', 'perc'});
        Tdist = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset', ...
            Rdist.phase < 3 & Rdist.cond == 1);
        Tacc = tapply(Racc,{'SN', 'region', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'});
        Tcross = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, ...
            'subset', Rdist.phase == 3 & Rdist.cond == 1);
        TmdsDist = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, ...
            'subset', Rdist.cond == 2);
        
        %%%Plot
        
        %         figure %%%General percent signal change overview (zoom to regions of interest)
        %         Tperc = tapply(Rperc,{'SN', 'region', 'phase'},{'perc', 'mean', 'name', 'perc'});
        %         colour={[0 0 0], [1 1 1]};
        %         regions = repmat({'lTha', 'lCau', 'lPut', 'lHip', 'rTha', 'rCau', 'rPut', 'rHip'}, 1, 12);
        %         barplot([Tperc.phase, Tperc.region], Tperc.perc, 'split', Tperc.phase, 'facecolor', colour);
        %         ylim(pscY)
        %         xticklabels(regions)
        %         ylabel('% signal change')
        %         title('Overview - activity increases')
        %         %-------------------------------------------------------------------------------%
        %
        %
        %         figure %%% General distance overview (zoom to regions of interest)
        %         Tdist = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset', ...
        %         Rdist.phase < 3 & Rdist.cond == 1);
        %         colour=rsaGB;
        %         regions = repmat({'lTha', 'lCau', 'lPut', 'lHip', 'rTha', 'rCau', 'rPut', 'rHip'}, 1, 12);
        %         barplot([Tdist.phase, Tdist.region], Tdist.dist, 'split', Tdist.phase, 'facecolor', colour);
        %         ylim(rsaY)
        %         xticklabels(regions)
        %         ylabel('Crossnobis dissimilarity')
        %         title('Overview - representational similarity analysis')
        %         %-------------------------------------------------------------------------------%
        %
        %
        %         figure %%% General decoding overview (zoom to regions of interest)
        %         Tacc = tapply(Racc,{'SN', 'region', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'});
        %         colour=decodeBRG;
        %         regions = repmat({'lTha', 'lCau', 'lPut', 'lHip', 'rTha', 'rCau', 'rPut', 'rHip'}, 1, 12);
        %         barplot([Tacc.phase, Tacc.cond, Tacc.region], Tacc.acc, 'split', Tacc.cond, 'facecolor', colour);
        %         drawline(0, 'dir', 'horz', 'linestyle', '- -')
        %         ylabel('Decoding accuracy')
        %         xticklabels(regions)
        %         title('Overview - Linear decoding accuracy')
        %         %-------------------------------------------------------------------------------%
        %
        %
        %         figure %%% General cross-phase distance overview
        %         Tcross = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, ...
        %             'subset', Rdist.phase == 3 & Rdist.cond == 1);
        %         colour={[0 0 0]};
        %         regions = repmat({'lTha', 'lCau', 'lPut', 'lHip', 'rTha', 'rCau', 'rPut', 'rHip'}, 1, 12);
        %         barplot(Tcross.region, Tcross.dist, 'facecolor', colour);
        %         ylim(crossY)
        %         xticklabels(regions)
        %         ylabel('Crossnobis dissimilarity')
        %         title('Overview - cross-phase RSA')
        %         %-------------------------------------------------------------------------------%
        %
        %
        %         figure %%% General MDS (PC2 & PC3) distance overview
        %         TmdsDist = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, ...
        %             'subset', Rdist.cond == 2);
        %         colour={[0 0 0]};
        %         regions = repmat({'lTha', 'lCau', 'lPut', 'lHip', 'rTha', 'rCau', 'rPut', 'rHip'}, 1, 12);
        %         barplot(TmdsDist.region, TmdsDist.dist, 'facecolor', colour);
        %         ylim(mdsDistY)
        %         xticklabels(regions)
        %         ylabel('Euclidean distance')
        %         title('Overview - MDS PC2 & PC3 RSA')
        %         %-------------------------------------------------------------------------------%
        
        
        %         figure %%% PSC Region subplots for prep/prod activity
        %         loopCount = 1;
        %         for i=[1 5 2 6 3 7 4 8]%plots left hem on the left, right hem on the right
        %             subplot(4,2,loopCount)
        %             loopCount = loopCount + 1;
        %             T = tapply(Rperc,{'SN', 'cond', 'phase'},{'perc', 'mean', 'name', 'perc'}, 'subset', ...
        %                 Rperc.region == i);
        %             colour={[0 0 0]}; %black %{[0 0.545 0.545], [1 0.647 0]}; %blue & orange
        %             lineplot([T.phase], T.perc, ...
        %                 'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
        %                 'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
        %             ylim(pscY)
        %             drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
        %             if i==1
        %                 ylabel('% signal change')
        %                 set(gca,'xticklabel',{'Prep', 'Prod'})
        %             else
        %                 ylabel('')
        %                 yticklabels({''})
        %                 set(gca,'xticklabel',{''})
        %             end
        %
        %             set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
        %             title(subcortStructs{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        %         end%for subcort region
        %         %-------------------------------------------------------------------------------%
        
        
        figure %%% PSC Region subplots for prep/prod activity - box plots
        loopCount = 1;
        for i=[1 5 2 6 3 7 4 8]%plots left hem on the left, right hem on the right
            subplot(4,2,loopCount)
            loopCount = loopCount + 1;
            T = tapply(Rperc,{'SN', 'cond', 'phase'},{'perc', 'mean', 'name', 'perc'}, 'subset', ...
                Rperc.region == i);
            colour={[0 0 0]}; %black %{[0 0.545 0.545], [1 0.647 0]}; %blue & orange
            myboxplot([T.phase], T.perc, ...
                'linewidth', 2, 'markersize', 5, 'markertype', 'o')%, ...
            %'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
            %'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
            ylim(pscY)
            drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
            if i==1
                ylabel('% signal change')
                set(gca,'xticklabel',{'Prep', 'Prod'})
            else
                ylabel('')
                yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(subcortStructs{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        end%for subcort region
        %-------------------------------------------------------------------------------%
        
        
        %         figure %%% RSA Region subplots for prep/prod distances
        %         loopCount = 1;
        %         for i=[1 5 2 6 3 7 4 8]%plots left hem on the left, right hem on the right
        %             subplot(4,2,loopCount)
        %             loopCount = loopCount + 1;
        %             T = tapply(Rdist,{'SN', 'cond', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset', ...
        %                 Rdist.region == i & Rdist.phase < 3 & Rdist.cond == 1);
        %             colour={[0 0 0]}; %black %{[0 0.545 0.545], [1 0.647 0]}; %blue & orange
        %             lineplot([T.phase], T.dist, ...
        %                 'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
        %                 'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
        %             ylim(rsaY)
        %             drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
        %             if i==1
        %                 ylabel('Crossnobis Distance')
        %                 set(gca,'xticklabel',{'Prep', 'Prod'})
        %             else
        %                 ylabel('')
        %                 yticklabels({''})
        %                 set(gca,'xticklabel',{''})
        %             end
        %
        %             set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
        %             title(subcortStructs{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        %         end%for subcort region
        %         %-------------------------------------------------------------------------------%
        
        figure %%% LDA Overall classifier Region subplots for prep/prod activity - box plots
        loopCount = 1;
        for i=[1 5 2 6 3 7 4 8]%plots left hem on the left, right hem on the right
            subplot(4,2,loopCount)
            loopCount = loopCount + 1;
            T = tapply(Racc,{'SN', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'}, 'subset', ...
                Racc.region == i & Racc.cond == 1);
            colour={[0 0 0]}; %black %{[0 0.545 0.545], [1 0.647 0]}; %blue & orange
            myboxplot([T.phase], T.acc, ...
                'linewidth', 2, 'markersize', 5, 'markertype', 'o', 'xtickoff')%, ...
            %'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
            %'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
            ylim([-4.1, 5.7])
            drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
            if i==1
                ylabel('Decoding accuracy (Z)')
                set(gca,'xticklabel',{'Prep', 'Prod'})
            else
                ylabel('')
                yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(subcortStructs{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        end%for subcort region
        %-------------------------------------------------------------------------------%
        
        
        figure %%% LDA Region subplots for prep/prod order/timing/integrated
        loopCount = 1;
        for i=[1 5 2 6 3 7 4 8]%plots left hem on the left, right hem on the right
            subplot(4,2,loopCount)
            loopCount = loopCount + 1;
            T = tapply(Racc,{'SN', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'},...
                'subset', Racc.region == i & Racc.cond ~= 5 & Racc.cond ~=1);
            colour=decodeBRG;
            lineplot([T.phase], T.acc, 'split', T.cond, ...
                'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
                'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
            ylim(ldaY)
            drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
            if i==1
                ylabel('Decoding accuracy (Z)')
                set(gca,'xticklabel',{'Prep', 'Prod'})
            else
                ylabel('')
                yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(subcortStructs{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        end%for subcort region
        %-------------------------------------------------------------------------------%
        
        
        figure %%% LDA Region subplots for prep/prod order/timing/integrated - box plots
        loopCount = 1;
        for i=[1 5 2 6 3 7 4 8]%plots left hem on the left, right hem on the right
            subplot(4,2,loopCount)
            loopCount = loopCount + 1;
            T = tapply(Racc,{'SN', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'}, ...
                'subset', Racc.region == i & Racc.cond ~= 5 & Racc.cond ~=1);
            colour=decodeBRG;
            myboxplot([T.phase], T.acc, 'split', T.cond, 'fillcolor', colour)%...
            %'markertype', 'o', 'markercolor', colour, 'markerfill', colou%, 'markersize', 5, ...
            %'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
            ylim(ldaY)
            drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
            if i==1
                ylabel('Decoding accuracy (Z)')
                set(gca,'xticklabel',{'Prep', '', '', 'Prod', '', ''})
            else
                ylabel('')
                %yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(subcortStructs{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        end%for subcort region
        %-------------------------------------------------------------------------------%
        
        
        %         figure %%% LDA Region subplots for prep/prod order/timing/integrated - box plots with old and new integrated decoders
        %         loopCount = 1;
        %         for i=[1 5 2 6 3 7 4 8]%plots left hem on the left, right hem on the right
        %             subplot(4,2,loopCount)
        %             loopCount = loopCount + 1;
        %             T = tapply(Racc,{'SN', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'}, 'subset', Racc.region == i);
        %             colour=decodeBRG;
        %             myboxplot([T.phase], T.acc, 'split', T.cond, 'fillcolor', colour)%...
        %             %'markertype', 'o', 'markercolor', colour, 'markerfill', colou%, 'markersize', 5, ...
        %             %'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
        %             ylim(ldaY)
        %             drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
        %             if i==1
        %                 ylabel('Decoding accuracy (Z)')
        %                 set(gca,'xticklabel',{'Prep', '', 'Old', 'New', 'Prod', '', 'Old', 'New'})
        %             else
        %                 ylabel('')
        %                 %yticklabels({''})
        %                 set(gca,'xticklabel',{''})
        %             end
        %             set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
        %             title(subcortStructs{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        %         end%for subcort region
        %         %-------------------------------------------------------------------------------%
        
        %%% RDMs and multi-dimensional scaling plots (for visualisation)
        labels = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        
        %Extract data
        for s=1:length(R.G)%for subj * region
            G(:,:,s)= R.G{s}; %extract representational dissimilarity matrix into 3D matrix
        end%for subj * region
        
        figure %Plot RDMs
        loopCount = 1;
        for i=[1 3 5 7 2 4 6 8]%plots left hem on the left, right hem on the right
            
            GRegion = G(:,:,R.region == loopCount); %extract region variance/covariance matrices
            GmRegion = mean(GRegion, 3); %mean across subjs
            
            subplot(4, 2, i) %%%RDM for each subcortical region
            ind=indicatorMatrix('allpairs',1:8); %matrix for all pairwise distances (k*(k-1))
            imagesc(rsa.rdm.squareRDM(diag(ind*GmRegion*ind'))); %display
            %multiplying variance/covariance by indicator matrix results in
            %dissimilarity values (crossnobis)
            title([subcortStructs{loopCount} ' RDM (crossnobis)'])
            colorbar
            axis image
            loopCount = loopCount + 1;
        end
        %----------------------------------------------------------------------------------------------%
        
        
        
        loopCount = 1; %Multi-dimensional scaling
        for i=[1 3 5 7 2 4 6 8]%plots left hem on the left, right hem on the right
            
            figure %plot MDS
            
            GRegion = G(:,:,R.region == loopCount); %extract region variance/covariance matrices
            
            lc=1;
            for j=1:length(GRegion) %get MDS for each subj in the region
                [SubjCOORD(:,:,lc),~] = pcm_classicalMDS(GRegion(:,:,lc));
                lc = lc+1;
            end
            COORD = mean(SubjCOORD, 3);
            
            %3D scatter plot
            scatterplot3(COORD(1:end,2),COORD(1:end,3),COORD(1:end,1),'split',spl, ... %here we plot PC 2&3 first because
                'markersize',ms, 'markercolor',rsaGB, 'markerfill',rsaGB, 'label',label);%PC 1 is just activity differences
            
            %%%Draw coloured lines between distinct order & timing conditions
            %Prep
            colors = decodeBRG{1}; %blue for order
            indx=[1 3]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[2 4]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            colors = decodeBRG{2}; %red for timing
            indx=[1 2]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[3 4]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            
            %Prod
            colors = decodeBRG{1}; %blue for order
            indx=[5 7]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[6 8]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            colors = decodeBRG{2}; %red for timing
            indx=[5 6]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[7 8]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            
            grid off
            hold on; plot3(0,0,0,'+','MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',ms+3, 'LineWidth',lw);
            hold off; xlabel('PC 2'); ylabel('PC 3'); zlabel('PC 1'); set(gca,'fontsize',12);
            %axis equal
            title([subcortStructs{loopCount} ' multi-dimensional scaling'])
            %xlim([-0.035, 0.035])
            %ylim([-0.027, 0.029])
            %zlim([-0.02, 0.3])
            view(30,30)
            loopCount = loopCount + 1;
        end
        
        %%%Plot Eigenvalues for PC1, PC2, and PC3
        eigen.vals   = [];%preallocate struct vars
        eigen.pcs    = [];
        eigen.region = [];
        eigen.sn     = [];
        
        figure
        loopCount = 1;
        for i = [1 5 2 6 3 7 4 8]%plots left hem on the left, right hem on the right %unique(Rdist.region)'
            eigValsAll = [Rdist.eigVals{Rdist.region == i & Rdist.cond == 2}]';
            eigVals    = reshape(eigValsAll(:,1:3), 72,1);%takes the first 3 PCs from eigValsAll and makes them one column
            eigPcs     = [ones(24,1); ones(24,1)*2; ones(24,1)*3];%variable to track which PC
            eigRegion  = ones(72,1)*i;%track region num
            eigSn      = repmat((1:24)', 3, 1);
            
            eigen.vals   = [eigen.vals; eigVals];%concat struct
            eigen.pcs    = [eigen.pcs; eigPcs];
            eigen.region = [eigen.region; eigRegion];
            eigen.sn     = [eigen.sn; eigSn];
            
            subplot(4,2,loopCount)
            myboxplot(eigen.pcs, eigen.vals, 'subset', eigen.region == i, 'xtickoff')
            ylim(eigY)
            if loopCount==1
                ylabel('Eigenvalues (a.u.)')
                set(gca,'xticklabel',{'PC1', 'PC2', 'PC3'})
            else
                ylabel('')
                %yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(subcortStructs{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
            loopCount = loopCount + 1;
        end
        %Save eigen values to excel
        varnames = {...
            'subj', ...
            'PC1_LThal', 'PC2_LThal', 'PC3_LThal', ...
            'PC1_LCaud', 'PC2_LCaud', 'PC3_LCaud', ...
            'PC1_LPut',  'PC2_LPut',  'PC3_LPut', ...
            'PC1_LHip',  'PC2_LHip',  'PC3_LHip', ...
            'PC1_RThal', 'PC2_RThal', 'PC3_RThal', ...
            'PC1_RCaud', 'PC2_RCaud', 'PC3_RCaud', ...
            'PC1_RPut',  'PC2_RPut',  'PC3_RPut', ...
            'PC1_RHip',  'PC2_RHip',  'PC3_RHip', ...
            };
        
        eigenData = [];
        for i = unique(eigen.region)'
            for j = unique(eigen.pcs)'
                eigenData(:,end+1) = eigen.vals(eigen.region == i & eigen.pcs == j);
            end
        end
        
        eigenTable = array2table([anaSubj' eigenData], 'VariableNames', varnames);
        writetable(eigenTable, fullfile(roiSubDir, 'eigenValsROI.xlsx'))
        
        %%%Save all results to respective, formatted excel files for SPSS stats
        %psc
        varnames = {...
            'subj', ...
            'prepLThal', 'prepLCaud', 'prepLPut', 'prepLHip', 'prepRThal', 'prepRCaud', 'prepRPut', 'prepRHip', ...
            'prodLThal', 'prodLCaud', 'prodLPut', 'prodLHip', 'prodRThal', 'prodRCaud', 'prodRPut', 'prodRHip', ...
            };
        percData = [];
        for i=unique(Tperc.phase)'
            for j=unique(Tperc.region)'
                percData(:,end+1) = Tperc.perc(Tperc.phase == i & Tperc.region == j);
            end
        end
        percTable = array2table([anaSubj' percData], 'VariableNames', varnames);
        writetable(percTable, fullfile(roiSubDir, 'percROIspss.xlsx'))
        
        %dist
        varnames = {...
            'subj', ...
            'prepLThal', 'prepLCaud', 'prepLPut', 'prepLHip', 'prepRThal', 'prepRCaud', 'prepRPut', 'prepRHip', ...
            'prodLThal', 'prodLCaud', 'prodLPut', 'prodLHip', 'prodRThal', 'prodRCaud', 'prodRPut', 'prodRHip', ...
            };
        distData = [];
        for i=unique(Tdist.phase)'
            for j=unique(Tdist.region)'
                distData(:,end+1) = Tdist.dist(Tdist.phase == i & Tdist.region == j);
            end
        end
        distTable = array2table([anaSubj' distData], 'VariableNames', varnames);
        writetable(distTable, fullfile(roiSubDir, 'distROIspss.xlsx'))
        
        %acc
        varnames = {...
            'subj', ...
            'ovrPrepLThal',  'ovrPrepLCaud',  'ovrPrepLPut',  'ovrPrepLHip',  'ovrPrepRThal',  'ovrPrepRCaud',  'ovrPrepRPut',  'ovrPrepRHip', ...
            'ovrProdLThal',  'ovrProdLCaud',  'ovrProdLPut',  'ovrProdLHip',  'ovrProdRThal',  'ovrProdRCaud',  'ovrProdRPut',  'ovrProdRHip', ...
            'ordPrepLThal',  'ordPrepLCaud',  'ordPrepLPut',  'ordPrepLHip',  'ordPrepRThal',  'ordPrepRCaud',  'ordPrepRPut',  'ordPrepRHip', ...
            'ordProdLThal',  'ordProdLCaud',  'ordProdLPut',  'ordProdLHip',  'ordProdRThal',  'ordProdRCaud',  'ordProdRPut',  'ordProdRHip', ...
            'tempPrepLThal', 'tempPrepLCaud', 'tempPrepLPut', 'tempPrepLHip', 'tempPrepRThal', 'tempPrepRCaud', 'tempPrepRPut', 'tempPrepRHip', ...
            'tempProdLThal', 'tempProdLCaud', 'tempProdLPut', 'tempProdLHip', 'tempProdRThal', 'tempProdRCaud', 'tempProdRPut', 'tempProdRHip', ...
            'intPrepLThal',  'intPrepLCaud',  'intPrepLPut',  'intPrepLHip',  'intPrepRThal',  'intPrepRCaud',  'intPrepRPut',  'intPrepRHip', ...
            'intProdLThal',  'intProdLCaud',  'intProdLPut',  'intProdLHip',  'intProdRThal',  'intProdRCaud',  'intProdRPut',  'intProdRHip', ...
            'intPrepOldLThal',  'intPrepOldLCaud',  'intPrepOldLPut',  'intPrepOldLHip',  'intPrepOldRThal',  'intPrepOldRCaud',  'intPrepOldRPut',  'intPrepOldRHip', ...
            'intProdOldLThal',  'intProdOldLCaud',  'intProdOldLPut',  'intProdOldLHip',  'intProdOldRThal',  'intProdOldRCaud',  'intProdOldRPut',  'intProdOldRHip', ...
            };
        accData = [];
        for i=unique(Tacc.cond)'
            for j=unique(Tacc.phase)'
                for k=unique(Tacc.region)'
                    accData(:,end+1) = Tacc.acc(Tacc.cond == i & Tacc.phase == j & Tacc.region == k);
                end
            end
        end
        accTable = array2table([anaSubj' accData], 'VariableNames', varnames);
        writetable(accTable, fullfile(roiSubDir, 'accROIspss.xlsx'))
        
        %cross
        varnames = {...
            'subj', ...
            'LTha', 'LCau', 'LPut', 'LHip', 'RTha', 'RCau', 'RPut', 'RHip', ...
            };
        crossData = [];
        for i=unique(Tcross.region)'
            crossData(:,end+1) = Tcross.dist(Tcross.region == i);
        end
        crossTable = array2table([anaSubj' crossData], 'VariableNames', varnames);
        writetable(crossTable, fullfile(roiSubDir, 'crossROIspss.xlsx'))
        
        %MDS Dist
        varnames = {...
            'subj', ...
            'LTha', 'LCau', 'LPut', 'LHip', 'RTha', 'RCau', 'RPut', 'RHip', ...
            };
        mdsDistData = [];
        for i=unique(TmdsDist.region)'
            mdsDistData(:,end+1) = TmdsDist.dist(Tcross.region == i);
        end
        mdsDistTable = array2table([anaSubj' mdsDistData], 'VariableNames', varnames);
        writetable(mdsDistTable, fullfile(roiSubDir, 'mdsDistROIspss.xlsx'))
        
        
    case 'subcortical_voxel_counts' %provides number of voxels for each subcortical structure
        
        sn = varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        voxCount = nan(1,length(subcortStructs));
        
        for j=1:length(subcortStructs) %and each subcortical structure...
            
            %cd into relevant subject directory
            cd([subcorticalDir, '/',  subj_name{sn}])
            
            %read in .nii file for each subcortical structure
            vol = spm_vol(['mr' subj_name{sn}, '_', subcortStructs{j}, '.nii']);
            readVol = spm_read_vols(vol);
            
            %count the number of non-zero elements
            voxelCount = nnz(readVol);
            
            %add the count to matrix
            voxCount(:,j) = voxelCount;
            
        end
        
        %add names of structures to top row of matrix and save to excel
        voxCount = num2cell(voxCount);
        voxCountSave = [subcortStructs; voxCount];
        
        writecell(voxCountSave, [subj_name{sn} '_voxelCounts'], 'FileType', 'spreadsheet')
    case 'subcortical_voxel_counts_group' %like case above, but produces one big group excel sheet. Quite slow...
        
        %set voxCount to size subj x structure for later loop
        subcortStructs = subcortStructs(2,:);%names from second row
        groupVoxCount = nan(length(anaSubj),length(subcortStructs));
        %loop counter for subject rows in excel sheet
        loopCount = 1;
        
        for i =anaSubj %for each subject...
            for j=1:length(subcortStructs) %and each subcortical structure...
                
                %cd into relevant subject directory
                cd([subcorticalDir, '/',  subj_name{i}])
                
                %read in .nii file for each subcortical structure
                vol = spm_vol(['mr' subj_name{i}, '_', subcortStructs{j}, '.nii']);
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
        
        cd(subcorticalDir)
        writecell(voxCountSave, 'group_voxelCounts', 'FileType', 'spreadsheet')
        
        
    case 'cerebellum_make_nii' %----------------- BEGINNING OF SUIT RSA&LDA (CEREBELLUM) -----------------%
        %Isolate the cerebellum of each participant - produces 'c_<source>_pcereb' (cerebellar mask) and '<source>_seg1/2' (grey and white matter respectively) images
        
        s=varargin{1};
        cd([baseDir '/imaging/anatomicals/' subj_name{s}]);
        disp(['suit isolating ' subj_name{s}])
        
        anatomical = {[subj_name{s} '_anatomical.nii']};
        
        suit_isolate_seg(anatomical, 'maskp', 0.2) %change maskp for probability value. Higher = tighter mask. Hand-correct using MRIcron if necessary
    case 'cerebellum_suit_normalise' %normalise the isolated cerebellum to the suit atlas - produces 'affine_<source>.mat' and 'u_a_<name>.nii'
        
        s=varargin{1};
        cd([baseDir '/imaging/anatomicals/' subj_name{s}]);
        disp(['suit_normalizing ' subj_name{s}])
        
        gray = {[subj_name{s} '_anatomical_seg1.nii']}; %grey and white matter images from previous suit stage
        white = {[subj_name{s} '_anatomical_seg2.nii']};
        isoMask = {['c_' subj_name{s} '_anatomical_pcereb.nii']}; %isolated cerebellum
        
        job.subjND.gray = gray; %put them all into a struct...
        job.subjND.white = white;
        job.subjND.isolation = isoMask;
        
        suit_normalize_dartel(job) %run the function with the struct as the input
    case 'cerebellum_make_mask'      %restrict area of analysis to grey matter - produces 'maskbrainSUIT.nii'
        
        s=varargin{1};
        
        if isfolder([baseDir '/imaging/suit/' subj_name{s}]) == 0
            mkdir([baseDir '/imaging/suit/' subj_name{s}])
        end
        
        mask=fullfile(glmDir, subj_name{s},'mask.nii');
        %suit=fullfile(anatDir, subj_name{s},['c_', subj_name{s},'_anatomical_pcereb.nii']); %pcereb holds all cerebellum-related regions to a value of 1...
        suit=fullfile(anatDir, subj_name{s},[subj_name{s}, '_anatomical_seg1.nii']); %whereas _seg1 is only grey matter and sets extra-cerebellar regions (e.g. pons) to values other than 1...
        omask=fullfile(suitDir, subj_name{s},'maskbrainSUIT.nii');
        
        spm_imcalc_ui({mask,suit},omask,'i1>0 & i2>0.999',{}); %so including a mask of 0.999 makes sure we only include cerebellar regions.
    case 'cerebellum_make_structs'      %extracts each cerebellar lobule from Lobules-SUIT.nii atlas file, and creates respective nii files
        
        regionValues = cell2mat(suitCBRegions(1,:));%take values from first row
        regionNames = suitCBRegions(2,:);%and names from second
        cd(fullfile(suitDir, 'atlasSUIT'))
        
        lobuleAtlasFile = 'Lobules-SUIT.nii';
        
        for i=1:length(regionNames)
            
            matlabbatch{1}.spm.util.imcalc.input = {lobuleAtlasFile};
            matlabbatch{1}.spm.util.imcalc.output = regionNames{i};
            matlabbatch{1}.spm.util.imcalc.expression = ['i1 == ', num2str(regionValues(i))];
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            
            spm_jobman('run',matlabbatch);
        end
    case 'cerebellum_structs_to_anat'   %uses SUIT to reslice segmented Lobules-SUIT atlas ROIs into native space
        
        sn = varargin{1};
        regionNames = suitCBRegions(2,:);%names from second row of variable
        cd(fullfile(suitDir, subj_name{sn}))
        
        inFiles  = strcat(fullfile(suitDir, 'atlasSUIT'), '/', regionNames, '.nii')';
        outFiles = strcat(fullfile(suitDir, subj_name{sn}), '/', subj_name{sn}, '_', regionNames, '.nii')';
        
        job.Affine     = {fullfile(anatDir, subj_name{sn}, ['Affine_' subj_name{sn} '_anatomical_seg1.mat'])};
        job.flowfield  = {fullfile(anatDir, subj_name{sn}, ['u_a_' subj_name{sn} '_anatomical_seg1.nii'])};
        job.resample   = inFiles;
        job.ref        = {fullfile(anatDir, subj_name{sn}, [subj_name{sn} '_anatomical.nii'])};
        job.out        = outFiles;
        
        suit_reslice_dartel_inv_RY(job)
    case 'cerebellum_reslice_structs'   %reslices cerebellar region niftis into functional scan resolution
        
        sn = varargin{1};
        
        refImageDir = [glmDir, '/', subj_name{sn}]; %directories for functional reference image (beta 1)
        suitImageDir = [suitDir, '/', subj_name{sn}]; %and for subcortical images
        regionNames = suitCBRegions(2,:);%names from second row of variable
        
        for r=1:length(regionNames) %loop through subcortical regions and reslice them to beta image resolution
            
            matlabbatch{1}.spm.spatial.realign.write.data = {[refImageDir, '/beta_0001.nii']; fullfile(suitImageDir, [subj_name{sn} '_' regionNames{r} '.nii'])};
            matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
            
            spm_jobman('run',matlabbatch);
            disp([subcortStructs{r} ' resliced'])
        end
    case 'cerebellum_funcmask_structs'  %masks cerebellar region niftis using functional grey matter mask
        
        s=varargin{1};
        regionNames = suitCBRegions(2,:);%names from second row of variable
        
        for i=1:length(regionNames)
            funMask=fullfile(suitDir, subj_name{s},'maskbrainSUIT.nii');
            omask=fullfile(suitDir, subj_name{s},['mr' subj_name{s} '_' regionNames{i} '.nii']); %output mask to be used in the future
            cbregion = fullfile(suitDir, subj_name{s},['r' subj_name{s} '_' regionNames{i} '.nii']);
            
            %spm_imcalc({funMask,cbregion},omask,'i1 >0.01 & i2 > 0.1',{}); %recorded activity in brain (grey + white matter)
            spm_imcalc({funMask,cbregion},omask,'i1 >0.01 & i2 > 0.3',{}); %recorded activity in brain (grey + white matter)
            %0.3 is the best threshold to make the resliced maps similar to
            %the anatomical segmentations
        end
    case 'cerebellum_make_ROIs'         %Area decoding - uses region toolbox to define R struct for prewhitening
        
        s=varargin{1};
        regionNames = suitCBRegions(2,:);%names from second row of variable
        R=cell(length(regionNames), 1);
        cd(fullfile(suitDir, subj_name{s}))
        
        for i=1:length(regionNames)%for each subcort region
            R{i} = region('roi_image',['mr' subj_name{s} '_' regionNames{i} '.nii'],1,regionNames{i});
            R{i}.name = [subj_name{s} '_' regionNames{i}]; %use toolbox to define, then add name
        end%for each subcort region
        
        R = region_calcregions(R);%because we reslice, we don't need voxelspace option
        %R = region_calcregions(R, 'voxelspace', fullfile(glmDir, subj_name{s}, 'beta_0001.nii'));
        
        out = [roiCbDir '/' subj_name{s} '_cerebellum_roi'];
        save(strjoin(out,''), 'R'); %save as participant file which holds all regions
    case 'cerebellum_preWhiten' %pre-whiten data from CB ROIs ready for RSA & LDA
        
        T = []; blueBear = varargin{1}; %s=varargin{1};
        
        for s=anaSubj
            fprintf('%d.',s); fprintf('\n')
            load([glmDir, '/', subj_name{s}, '/', 'SPM.mat'],            'SPM')
            load([roiCbDir, '/', subj_name{s}, '_cerebellum_roi.mat'], 'R')
            cd(fullfile(suitDir, subj_name{s}))
            
            V = SPM.xY.VY;
            
            %replace fnames for bluebear compatibility
            if blueBear == 1
                for i=1:length(V)
                    V(i).fname = strrep(V(i).fname,'\','/');
                    V(i).fname = strrep(V(i).fname,'Z:','/rds/projects/k/kornyshk-kornyshevalab');
                end
            else
            end
            
            for r=1:length(R)
                Y = region_getdata(V,R{r});
                
                percMove = region_getdata(spm_vol(fullfile(glmDir, subj_name{s}, 'perc_0001.nii')), R{r}); %percent signal change
                percPrep = region_getdata(spm_vol(fullfile(glmDir, subj_name{s}, 'perc_0002.nii')), R{r});
                
                [betaW,resMS,~,beta,~,~,snr] = noiseNormalizeBeta_RY(Y,SPM);
                
                S.betaW    = {betaW}; %multivariate pre-whiten
                S.betaUW   = {bsxfun(@rdivide,beta,sqrt(resMS))}; %univariate pre-whiten
                S.betaRAW  = {beta}; %Raw
                S.resMS    = {resMS}; %residual
                S.snr      = snr;
                S.percMov  = {percMove};
                S.percPrep = {percPrep};
                S.SN       = s;
                S.region   = r;
                
                T = addstruct(T,S);
                fprintf('%d.',r)
            end
            fprintf('\n');
        end
        
        save(fullfile(roiCbDir, 'preWhitened_betas.mat'),'-struct','T');
    case 'cerebellum_calculate' %calculates overall distance & factorial classification acc for each CB ROI
        
        cd(roiCbDir)
        R = load('preWhitened_betas.mat');
        saveDir = fullfile(rsaDir, 'cerebellum');
        
        nrruns = length(run); nCond = 8; nClassifiers = 10;
        
        %%% Distances: prepare condition and partition (run) vectors
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 0 0]; %prep
        prep=[repmat(prep,1,nrruns) runBSL];%1 x nBeta, 1 2 3 4 = prep sequences 1:4
        
        prod      =[5 0 0 0 6 0 0 0 7 0 0 0 8 0 0 0 0 0 0 0]; %prod
        prod=[repmat(prod,1,nrruns) runBSL];%1 x nBeta, 5 6 7 8 = prod sequences 1:4
        
        condVec = prep + prod; condVec = condVec'; % conditions, including no interest regressors as 0
        partVec = double(condVec > 0); %assign run number to beta
        for i=1:nCond %assign run numbers to conds, ignore no-interest betas
            partIdx = find(condVec == i);
            for j=1:length(run)
                partVec(partIdx(j)) = j;
            end
        end
        
        %%% Classifiers: prepare beta selections (prep and prod separately)
        prepBetas = condVec < 5 & condVec > 0;
        prodBetas = condVec > 4;
        
        %%% loop through prewhitened data, calculate crossnobis dissimilarities
        %pre-allocate output variables
        R.d         = NaN(length(R.SN), nCond * (nCond - 1) / 2); %pairwise distance measures
        R.Sig       = cell(length(R.SN), 1); %covariance matrix of the beta estimates across different imaging runs
        R.G         = cell(length(R.SN), 1); %second moment matrix
        R.matrix    = cell(length(R.SN), 1); % Pairwise contrast matrix
        R.acc       = NaN(length(R.SN), nClassifiers); %LDA accuracy - 1ordPrep, 2timPrep, 3intPrep, 4ordProd, 5timProd, 6intProd
        
        for s = anaSubj%for subj
            load(fullfile(glmDir, subj_name{s}, 'SPM'), 'SPM') %load SPM design matrix for distance function
            for r = unique(R.region)'%for region
                B    = R.betaW{R.region == r & R.SN == s};
                BRAW = R.betaRAW{R.region == r & R.SN == s};
                resMS = R.resMS{R.region == r & R.SN == s};
                
                %%%Distances
                [d, Sig] = rsa.distanceLDC(B, partVec, condVec, SPM.xX.X);
                [G,~]    = pcm_estGCrossval(B,partVec,condVec, 'X', SPM.xX.X);
                matrix   = indicatorMatrix('allpairs',1:nCond); % Pairwise contrast matrix
                
                %%%Classification
                [acc(1), acc(2), acc(3), acc(4), acc(5)]  = prepProdSimu_classify(BRAW(prepBetas, :)); %factorial classify prep sequences
                [acc(6), acc(7), acc(8), acc(9), acc(10)] = prepProdSimu_classify(BRAW(prodBetas, :)); %and prod sequences
                %acc: 1=ovrPrep, 2=ordPrep, 3=tempPrep, 4=intPrep, 5=intPrepSubtract
                %     6=ovrProd, 7=ordProd, 8=tempProd, 9=intProd, 10=intProdSubtract
                
                %%%Percent signal change
                meanPercMov  = mean(R.percMov{R.region == r & R.SN == s});
                meanPercPrep = mean(R.percPrep{R.region == r & R.SN == s});
                
                %%%Variable assignment
                R.d     (R.region == r & R.SN == s, :) = d;
                R.Sig   {R.region == r & R.SN == s}    = Sig;
                R.G     {R.region == r & R.SN == s}    = G;
                R.matrix{R.region == r & R.SN == s}    = matrix;
                R.acc   (R.region == r & R.SN == s, :) = acc;
                R.mov   (R.region == r & R.SN == s, :) = meanPercMov;
                R.prep  (R.region == r & R.SN == s, :) = meanPercPrep;
            end%for subj
        end%for region
        
        %%% Identify and extract distance values for preparation, production, and cross-phase
        phaseVector = [...
            1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...  %pairwise contrast index (1s are contrasts within preparation)
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1; ...  %^ index (1s are contrasts within production)
            0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0  ...  %^ index (1s are contrasts across phases)
            ];
        
        %%%Loop over subjs, conditions, and phases, storing distances and
        %%%accuracies and signal change as a struct to plot later
        loopCounter = 1;
        
        %%% Distance
        for i=1:length(R.betaW)%for data points
            for j=1:size(phaseVector, 1)%for phase
                Rdist.dist(loopCounter,1)   = mean(R.d(i, phaseVector(j,:) == 1));
                Rdist.phase(loopCounter,1)  = j;
                Rdist.cond(loopCounter,1)   = 1;
                Rdist.SN(loopCounter,1)     = R.SN(i);
                Rdist.region(loopCounter,1) = R.region(i);
                
                loopCounter = loopCounter + 1;
            end%for phase
        end%for data points
        
        loopCounter = 1;
        
        %%% Decoding
        for i=1:length(R.betaW)%for data points
            for j=1:size(R.acc, 2)%for classifiers
                Racc.acc(loopCounter,1)   = R.acc(i, j);
                if j < 6%if prep classifier
                    Racc.phase(loopCounter,1) = 1;
                    Racc.cond(loopCounter,1)  = j;
                elseif j > 5%if prod classifier
                    Racc.phase(loopCounter,1) = 2;
                    Racc.cond(loopCounter,1)  = j - 5;
                end
                Racc.SN(loopCounter,1)    = R.SN(i);
                Racc.region(loopCounter,1) = R.region(i);
                
                loopCounter = loopCounter + 1;
            end%for classifiers
        end%for subj
        
        R;
        loopCounter = 1;
        
        %%% Distance
        for i=1:length(R.betaW)%for data points
            for j=1:2%for phase
                
                if j == 1
                    Rperc.perc(loopCounter,1) = mean(R.percPrep{i});
                elseif j == 2
                    Rperc.perc(loopCounter,1) = mean(R.percMov{i});
                end
                Rperc.phase(loopCounter,1)    = j;
                Rperc.cond(loopCounter,1)     = 1;
                Rperc.SN(loopCounter,1)       = R.SN(i);
                Rperc.region(loopCounter,1)   = R.region(i);
                
                loopCounter = loopCounter + 1;
            end%for phase
        end%for data points
        
        if ~isfolder(saveDir)
            mkdir(saveDir)
        end
        
        save(fullfile(saveDir, 'cbRoiDistances.mat'), 'R', 'Rdist', 'Racc', 'Rperc')
    case 'cerebellum_mds_distances' %cross-phase MDS distance when PC1 is ignored - saves to subRoiDistances_withmds.mat
        
        %%%Load
        dataFile = fullfile(rsaDir, 'cerebellum', 'cbRoiDistances.mat');
        load(dataFile, 'R', 'Rdist', 'Racc', 'Rperc')%load it
        
        %Set eucScaling and pcs variables to same size as Rdist
        Rdist.eucScaling = nan(length(Rdist.SN),1);
        Rdist.pcs        = cell(size(Rdist.SN));
        Rdist.eigVals    = cell(size(Rdist.SN));
        
        %%%Set save directory - we re-save the distances to include
        %%%cross-phase MDS
        saveDir = fullfile(rsaDir, 'cerebellum');
        
        %%% RDMs and multi-dimensional scaling plots (for visualisation)
        labels = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        
        %Extract data
        for s=1:length(R.G)%for subj & region
            G(:,:,s)= R.G{s}; %extract representational dissimilarity matrix from cells into 3D matrix
        end%for subj * region
        
        for i=unique(R.region)'%for region
            COORD = [];
            GRegion = G(:,:,R.region == i); %extract region variance/covariance matrices
            
            for j=1:length(GRegion) %get MDS for each subj in the region
                [COORD(:,:,j),eigVals(:,j)] = pcm_classicalMDS(GRegion(:,:,j));
            end
            
            PCs = COORD(:,2:3,:); %PCs 2 and 3
            
            %Euclidean distance between phases within sequences for PC2
            %and PC3
            for j=1:length(PCs)%for subj
                for k=1:4%for within-sequence prep vs prod
                    eucPcDist(k) = pdist([PCs(k,:,j)' PCs(k+4,:,j)'], 'euclidean');
                end%for within-sequence prep vs prod
                
                %pairwise distance measures k(k-1) / 2
                %within prep
                prepDists(1) = pdist([PCs(1,:,j)', PCs(2,:,j)'],'euclidean');
                prepDists(2) = pdist([PCs(1,:,j)', PCs(3,:,j)'],'euclidean');
                prepDists(3) = pdist([PCs(1,:,j)', PCs(4,:,j)'],'euclidean');
                prepDists(4) = pdist([PCs(2,:,j)', PCs(3,:,j)'],'euclidean');
                prepDists(5) = pdist([PCs(2,:,j)', PCs(4,:,j)'],'euclidean');
                prepDists(6) = pdist([PCs(3,:,j)', PCs(4,:,j)'],'euclidean');
                %within prod
                prodDists(1) = pdist([PCs(5,:,j)', PCs(6,:,j)'],'euclidean');
                prodDists(2) = pdist([PCs(5,:,j)', PCs(7,:,j)'],'euclidean');
                prodDists(3) = pdist([PCs(5,:,j)', PCs(8,:,j)'],'euclidean');
                prodDists(4) = pdist([PCs(6,:,j)', PCs(7,:,j)'],'euclidean');
                prodDists(5) = pdist([PCs(6,:,j)', PCs(8,:,j)'],'euclidean');
                prodDists(6) = pdist([PCs(7,:,j)', PCs(8,:,j)'],'euclidean');
                %%%FIND A BETTER WAY TO CODE THIS ^
                
                prepSize = mean(prepDists);
                prodSize = mean(prodDists);
                avgSize  = mean([prepSize, prodSize]);
                
                %Collect variables to later add to RDist
                Rdist.dist(end+1,1)        = mean(eucPcDist);
                Rdist.eucScaling(end+1,1)  = avgSize;
                Rdist.phase(end+1,1)       = 0;
                Rdist.cond(end+1,1)        = 2;
                Rdist.SN(end+1,1)          = anaSubj(j);
                Rdist.region(end+1,1)      = i;
                Rdist.pcs{end+1,1}         = COORD(:,:,j);
                Rdist.eigVals{end+1,1}     = eigVals(:,j);
            end%for subj
        end%for region
        
        save(fullfile(saveDir, 'cbRoiDistances_withmds.mat'), 'R', 'Rdist', 'Racc', 'Rperc')
    case 'cerebellum_make_search'      %makes 160 voxel searchlight contained to CB grey matter
        
        s=varargin{1};
        
        radius=16;
        numVox=160;
        
        cd(fullfile(suitDir, subj_name{s}));
        V=spm_vol('maskbrainSUIT.nii'); %preceded by case suit_make_mask
        X=spm_read_vols(V);
        [i,j,k]=ind2sub(size(X),find(X~=0));
        vox=[i j k];
        
        L = [];
        [L.LI,L.voxmin,L.voxmax,L.n]=lmva_voxelselection(vox(:,:)',vox',[radius numVox],V.mat,V.dim,[],'mva160_numvoxSUIT.nii');
        L.voxel = vox;
        save('volsearch160SUIT.mat', 'L')
    case 'cerebellum_run_search_RSA'   %runs RSA searchlight, generates crossnobis dissimilarity values
        
        s=varargin{1}; blueBear=varargin{2};
        nrruns = length(run);
        
        cd(fullfile(suitDir, subj_name{s}))
        
        %%% Searchlight file
        load(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'), 'L');
        
        %%% SPM file
        spmDir = fullfile(glmDir, subj_name{s});
        load(fullfile(spmDir, 'SPM'), 'SPM');
        
        %replace fnames for bluebear compatibility
        if blueBear == 1
            for i=1:length(SPM.xY.VY)
                SPM.xY.VY(i).fname = strrep(SPM.xY.VY(i).fname,'\','/');
                SPM.xY.VY(i).fname = strrep(SPM.xY.VY(i).fname,'Z:','/rds/projects/k/kornyshk-kornyshevalab');
            end
        else
        end
        
        %%% prepare condition and partition (run) vectors
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 0 0]; %prep
        prep=[repmat(prep,1,nrruns) runBSL];%1 x nBeta, 1 2 3 4 = prep sequences 1:4
        
        prod      =[5 0 0 0 6 0 0 0 7 0 0 0 8 0 0 0 0 0 0 0]; %prod
        prod=[repmat(prod,1,nrruns) runBSL];%1 x nBeta, 5 6 7 8 = prod sequences 1:4
        
        condVec = prep + prod; condVec = condVec'; % conditions, including no interest regressors as 0
        
        %%% Run searchlight function on whole CB
        rsa.runSearchlightLDC_RY(L, SPM, 'spmDir', spmDir, 'conditionVec', condVec, ...
            'analysisName', [subj_name{s}, 'RSA_All'], 'outDir', fullfile(suitDir, subj_name{s}))
    case 'cerebellum_run_spatMov_LDA'  %LDA searchlight order during movement
        sn=varargin{1};
        
        for s=sn
            
            cd(fullfile(glmDir, subj_name{s}));
            load('SPM.mat', 'SPM');
            load(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
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
            end
            
            [~,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
        end
    case 'cerebellum_run_spatPrep_LDA' %LDA searchlight order during prep
        sn=varargin{1};
        
        for s=sn
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
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
                
            end
            
            [~,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
        end
    case 'cerebellum_run_tempMov_LDA'  %LDA searchlight timing during movement
        sn=varargin{1};
        
        for s=sn
            
            s=varargin{1};
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
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
                
            end
            
            [~,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
        end
    case 'cerebellum_run_tempPrep_LDA' %LDA searchlight timing during prep
        sn=varargin{1};
        
        for s=sn
            
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            load(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'), 'L') %load searchlight
            L.vox = L.voxel; L = rmfield(L, 'voxel');
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
                
            end
            
            [~,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            lmva_spm(L,Pselect,out,@combinedclass,'params',{c,run,train,test});
        end
    case 'cerebellum_run_intMov_LDA'   %LDA 'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, O1, O2 and classifies residual)
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        load(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'), 'L') %load searchlight
        L.vox = L.voxel; L = rmfield(L, 'voxel');
        nrruns=length(SPM.nscan);
        
        
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
        prod=[repmat(prod,1,nrruns) runBSL];
        
        c=repmat(1:4,1,nrruns); 
        run=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6]; % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
        out = {fullfile(suitDir, subj_name{s}, [subj_name{s}, '_accuracy_Int_160_Mov.nii'])};
        
        
        % Generate column indices for Cross-validation, where
        % cell i contains column indices of the respective test and
        % train set
        for i=1:nrruns
            test{i}=find(run==i);  %fprintf('test:'); display(test{i}');
            train{i}=find(run~=i);  %fprintf('train:'); display(train{i}');
        end
        
        [~,col] = find(prod>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end
        
        lmva_spm(L,Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
    case 'cerebellum_run_intPrep_LDA'  %LDA 'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, O1, O2 and classifies residual)
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
        load SPM;
        load(fullfile(suitDir,subj_name{s},'volsearch160SUIT.mat'), 'L') %load searchlight
        L.vox = L.voxel; L = rmfield(L, 'voxel');
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
        end
        
        [~,col] = find(prep>0);
        P={SPM.Vbeta(1:126).fname}';
        Pselect=[];
        for i=1:length(col)
            Pselect{i,1}= P{col(i)};
        end
        
        lmva_spm(L,Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
    case 'cerebellum_zValue_LDA'       %z value at 50% chance level (include int because of 2x2 design subtraction feature, see paper)
        
        s=varargin{1};
        cd(fullfile(suitDir,subj_name{s}));
        
        %zValue the order and timing maps
        takeOneOutIter=2;
        numTests=6;
        numCat=2;
        mu=1/numCat; %mu=0.5;
        N=numTests*numCat*takeOneOutIter;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        images= {...
            '_accuracy_Spat_160_Mov','_accuracy_Spat_160_Prep',...
            '_accuracy_Temp_160_Mov','_accuracy_Temp_160_Prep',...
            };
        
        outimages={...
            '_zacc_Spat_160_Mov','_zacc_Spat_160_Prep',...
            '_zacc_Temp_160_Mov','_zacc_Temp_160_Prep',...
            };
        
        
        for j=1:numel(images)
            input_image= fullfile(suitDir,subj_name{s},[subj_name{s} images{j} '.nii']);
            output_image= fullfile(suitDir,subj_name{s},[subj_name{s} outimages{j} '.nii']);
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
        end
        
        %zValue the integrated map
        numTests=6;
        numCat=4;
        mu=1/numCat; %mu=0.25;
        N=numTests*numCat;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        images= {...
            '_accuracy_Int_160_Mov', '_accuracy_Int_160_Prep'...
            };
        
        outimages={...
            '_zacc_Int_160_Mov', '_zacc_Int_160_Prep'...
            };
        
        
        for j=1:numel(images)
            input_image= fullfile(suitDir,subj_name{s},[subj_name{s} images{j} '.nii']);
            output_image= fullfile(suitDir,subj_name{s},[subj_name{s} outimages{j} '.nii']);
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask! z_accuracy=(accuracy-mu)/sigma;
        end
    case 'cerebellum_reslice_contrast' %reslice individual contrast maps into SUIT space
        
        sn=varargin{1};
        cd([baseDir '/imaging/suit/' subj_name{sn}]);
        disp(['suit_reslicing_contrast ' subj_name{sn}])
        
        inDir = fullfile(glmDir, subj_name{sn}); %path to where data is stored (to be normalised)
        
        contrasts = {...
            'con_0001', 'con_0002', 'con_0003', 'con_0004', 'con_0005', 'con_0006', ...
            'perc_0001', 'perc_0002'...
            };
        
        outDir = fullfile(suitDir, subj_name{sn});
        
        % prepare files for input
        affine = {[anatDir '/' subj_name{sn} '/' 'Affine_' subj_name{sn} '_anatomical_seg1.mat']};
        flowfield = {[anatDir '/' subj_name{sn} '/' 'u_a_' subj_name{sn} '_anatomical_seg1.nii']};
        
        dataFiles = cell(length(contrasts),1);
        for i=1:length(contrasts)
            dataFiles{i} = fullfile(inDir, [contrasts{i} '.nii']);
        end
        
        mask = {[anatDir '/' subj_name{sn} '/' 'c_' subj_name{sn} '_anatomical_pcereb.nii']};
        
        outFiles = cell(length(contrasts),1);
        for i=1:length(contrasts)
            outFiles{i} = fullfile(outDir, [contrasts{i} '.nii']);
        end
        
        %%% prepare struct for function
        job.subj.affineTr = affine; %fill job.subj. struct with respective items
        job.subj.flowfield = flowfield;
        job.subj.resample = dataFiles;
        job.subj.mask = mask;
        job.subj.outname = outFiles;
        
        %function
        suit_reslice_dartel(job)
    case 'cerebellum_smooth'           %smooth dissimilarity, decoding, and psc/contrast suit maps
        
        s=varargin{1};
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} 'RSA_All_LDC.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} 'RSA_ALL_sLDC.nii']);
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
        
        comb=fullfile(suitDir, subj_name{s}, 'con_0001.nii'); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0001.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},'con_0002.nii'); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0002.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s}, 'con_0003.nii'); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0003.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},'con_0004.nii'); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0004.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},'con_0005.nii'); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0005.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},'con_0006.nii'); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_scon_0006.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},'perc_0001.nii'); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_sperc_0001.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
        
        comb=fullfile(suitDir, subj_name{s},'perc_0002.nii'); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} '_sperc_0002.nii']);
        spm_smooth(comb,scomb,[4 4 4]); %smooth with 4mm kernel
    case 'cerebellum_normalise_contrast' %normalise contrast & psc into suit space
        
        sn=varargin{1};
        cd([baseDir '/imaging/suit/']);
        
        if isfolder(suitGroupDir) == 0
            mkdir(suitGroupDir)
        end
        disp(['suit_reslicing ' subj_name{sn}])
        
        inDir = [suitDir '/' subj_name{sn} '/']; %path to where data is stored (to be normalised)
        outDir = suitGroupDir;
        filenames = {'scon_0001', 'scon_0002', 'scon_0003', 'scon_0004', 'scon_0005', 'scon_0006', ...
            'sperc_0001', 'sperc_0002'};
        
        % prepare files for input
        affine = {[anatDir '/' subj_name{sn} '/' 'Affine_' subj_name{sn} '_anatomical_seg1.mat']};
        flowfield = {[anatDir '/' subj_name{sn} '/' 'u_a_' subj_name{sn} '_anatomical_seg1.nii']};
        
        dataFiles = cell(length(filenames),1); %loop to put all full input file directories into a cell
        for i=1:length(filenames)
            dataFiles{i} = [inDir, subj_name{sn}, '_', filenames{i}, '.nii'];
        end
        
        mask = {[anatDir '/' subj_name{sn} '/' 'c_' subj_name{sn} '_anatomical_pcereb.nii']};
        
        outFiles = cell(length(filenames),1);
        for i=1:length(filenames)
            outFiles{i} = [outDir, '/', filenames{i}, '_', subj_name{sn}, '.nii'];
        end
        
        %%% prepare struct for function
        job.subj.affineTr = affine; %fill job.subj. struct with respective items
        job.subj.flowfield = flowfield;
        job.subj.resample = dataFiles;
        job.subj.mask = mask;
        job.subj.outname = outFiles;
        
        %function
        suit_reslice_dartel(job)
    case 'cerebellum_calc_dissimilarity_maps' %extract overall distances from RSA_ALL_sLDC.nii
        %each volume of the searchlight corresponds to a pairwise
        %dissimilarity measure between sequences.
        
        s=varargin{1};
        cd(fullfile(suitDir, subj_name{s}))
        
        %%% Identify and extract values for overall, order, and timing
        %%% within preparation, production, and cross-phase
        prepCols       = [1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';   %^ index (1s are contrasts within preparation)
        prodCols       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1]';   %^ index (1s are contrasts within production)
        crossPhaseCols = [0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0]';   %^ index (1s are contrasts across phases within sequences)
        
        conds = {...
            'overall_prep', 'overall_prod', 'overall_cross', ...
            prepCols == 1,  prodCols == 1,  crossPhaseCols == 1, ...
            };
        
        for j=1:length(conds)
            vol = spm_vol([subj_name{s} 'RSA_ALL_sLDC.nii']);
            
            Vi = vol(conds{2,j}, :);
            Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
            Vo.fname = [subj_name{s} '_LDC_' conds{1,j} '.nii'];
            Vo.n = [1 1];
            express = 'mean(X)';
            flags.dmtx = 1;
            
            spm_imcalc(Vi, Vo, express, flags)
        end
    case 'cerebellum_normalise_RSA' %normalisation into suit space - RSA
        
        sn=varargin{1};
        cd([baseDir '/imaging/suit/']);
        
        if isfolder(suitGroupDir) == 0
            mkdir(suitGroupDir)
        end
        disp(['suit_reslicing ' subj_name{sn}])
        
        inDir = [suitDir '/' subj_name{sn} '/']; %path to where data is stored (to be normalised)
        outDir = suitGroupRSADir;
        filenames = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        
        % prepare files for input
        affine = {[anatDir '/' subj_name{sn} '/' 'Affine_' subj_name{sn} '_anatomical_seg1.mat']};
        flowfield = {[anatDir '/' subj_name{sn} '/' 'u_a_' subj_name{sn} '_anatomical_seg1.nii']};
        
        dataFiles = cell(length(filenames),1); %loop to put all full input file directories into a cell
        for i=1:length(filenames)
            dataFiles{i} = [inDir, subj_name{sn}, '_LDC_' filenames{i}, '.nii'];
        end
        
        mask = {[anatDir '/' subj_name{sn} '/' 'c_' subj_name{sn} '_anatomical_pcereb.nii']};
        
        outFiles = cell(length(filenames),1);
        for i=1:length(filenames)
            outFiles{i} = [outDir, '/', filenames{i}, '_', subj_name{sn}, '.nii'];
        end
        
        %%% prepare struct for function
        job.subj.affineTr = affine; %fill job.subj. struct with respective items
        job.subj.flowfield = flowfield;
        job.subj.resample = dataFiles;
        job.subj.mask = mask;
        job.subj.outname = outFiles;
        
        %function
        suit_reslice_dartel(job)
    case 'cerebellum_normalise_LDA' %normalisation into suit space - LDA
        
        sn=varargin{1};
        cd([baseDir '/imaging/suit/']);
        
        if isfolder(suitGroupDir) == 0
            mkdir(suitGroupDir)
        end
        disp(['suit_reslicing ' subj_name{sn}])
        
        inDir = [suitDir '/' subj_name{sn} '/']; %path to where data is stored (to be normalised)
        outDir = suitGroupDir;
        filenames = {...
            'szacc_Spat_160_Prep', 'szacc_Spat_160_Mov',...
            'szacc_Temp_160_Prep', 'szacc_Temp_160_Mov',...
            'szacc_Int_160_Prep',  'szacc_Int_160_Mov'...
            };
        
        % prepare files for input
        affine = {[anatDir '/' subj_name{sn} '/' 'Affine_' subj_name{sn} '_anatomical_seg1.mat']};
        flowfield = {[anatDir '/' subj_name{sn} '/' 'u_a_' subj_name{sn} '_anatomical_seg1.nii']};
        
        dataFiles = cell(length(filenames),1); %loop to put all full input file directories into a cell
        for i=1:length(filenames)
            dataFiles{i} = [inDir, subj_name{sn}, '_', filenames{i}, '.nii'];
        end
        
        mask = {[anatDir '/' subj_name{sn} '/' 'c_' subj_name{sn} '_anatomical_pcereb.nii']};
        
        outFiles = cell(length(filenames),1);
        for i=1:length(filenames)
            outFiles{i} = [outDir, '/', filenames{i}, '_', subj_name{sn}, '.nii'];
        end
        
        %%% prepare struct for function
        job.subj.affineTr = affine; %fill job.subj. struct with respective items
        job.subj.flowfield = flowfield;
        job.subj.resample = dataFiles;
        job.subj.mask = mask;
        job.subj.outname = outFiles;
        
        %function
        suit_reslice_dartel(job)
    case 'cerebellum_group_avg_RSA'       %group average RSA maps
        
        if ~isfolder(fullfile(suitGroupRSADir, 'average'))
            mkdir(fullfile(suitGroupRSADir, 'average'))
        end
        cd(fullfile(suitGroupRSADir, 'average'))
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        
        for i=1:length(conds)
            loopCount = 1;
            for s = anaSubj
                Vi(loopCount) = spm_vol(fullfile(suitGroupRSADir, [conds{i} '_' subj_name{s} '.nii']));
                loopCount = loopCount + 1;
            end
            Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
            Vo.fname = ['avg_' conds{i} '_LDC.nii'];
            Vo.n = [1 1];
            express = 'mean(X)';
            flags.dmtx = 1;
            
            spm_imcalc(Vi, Vo, express, flags)
        end
    case 'cerebellum_group_avg_LDA'       %group average LDA maps
        
        if ~isfolder(fullfile(suitGroupDir, 'average'))
            mkdir(fullfile(suitGroupDir, 'average'))
        end
        cd(fullfile(suitGroupDir, 'average'))
        
        images = {...
            'szacc_Spat_160_Mov';'szacc_Spat_160_Prep';...
            'szacc_Temp_160_Mov';'szacc_Temp_160_Prep';...
            'szacc_Int_160_Mov'; 'szacc_Int_160_Prep'...
            };
        
        conds = {...
            'spat_mov', 'spat_prep',...
            'temp_mov', 'temp_prep', ...
            'int_mov', 'int_prep', ...
            };
        
        for i=1:length(images)
            loopCount = 1;
            for s = anaSubj
                Vi(loopCount) = spm_vol(fullfile(suitGroupDir, [images{i} '_' subj_name{s} '.nii']));
                loopCount = loopCount + 1;
            end
            Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
            Vo.fname = ['avg_' conds{i} '_LDA.nii'];
            Vo.n = [1 1];
            express = 'mean(X)';
            flags.dmtx = 1;
            
            spm_imcalc(Vi, Vo, express, flags)
        end
    case 'cerebellum_group_avg_psc'       %group average percent signal change & contrast maps
        
        if ~isfolder(fullfile(suitGroupDir, 'average'))
            mkdir(fullfile(suitGroupDir, 'average'))
        end
        cd(fullfile(suitGroupDir, 'average'))
        
        images = {...
            'scon_0001';'scon_0002';'scon_0003';'scon_0004';'scon_0005';'scon_0006';...
            'sperc_0001';'sperc_0002'...
            };
        
        conds = {...
            'con_mov', 'con_prep', 'con_error', 'con_prepvprod', 'con_prodvprep', 'con_rest', ...
            'perc_mov', 'perc_prep'...
            };
        
        for i=1:length(images)
            loopCount = 1;
            for s = anaSubj
                Vi(loopCount) = spm_vol(fullfile(suitGroupDir, [images{i} '_' subj_name{s} '.nii']));
                loopCount = loopCount + 1;
            end
            Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
            Vo.fname = ['avg_' conds{i} '.nii'];
            Vo.n = [1 1];
            express = 'mean(X)';
            flags.dmtx = 1;
            
            spm_imcalc(Vi, Vo, express, flags)
        end
        
    case 'cerebellum_group_randomeffects_RSA' %random effects - RSA
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        
        dataDir = strcat('RSA_', conds);
        
        contrastN = length(dataDir);
        images = repmat(conds', 1, length(anaSubj));
        subNii = repmat (subj_name(anaSubj), length(dataDir), 1);
        fileName = strcat (images, '_', subNii, '.nii');  %%Concatenate contrast files and subject names
        
        for i=1:contrastN  %%Loop across contrasts
            glmscndDir = fullfile(suitGroupRSADir, dataDir(i));
            matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;  %Adjust directory
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (suitGroupRSADir, fileName(i,:))';  %%Select files from matrix
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
    case 'cerebellum_group_randomeffects_LDA' %random effects - LDA
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
        %(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+16+i17+i18+i19+i20+i21+i22+i23+i24)/24
    case 'cerebellum_group_estimate_RSA'      %estimate RSA group
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        contrastN = length(conds);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile (suitGroupRSADir, ['RSA_' conds{i}], 'SPM.mat');  %Adjust directory
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
    case 'cerebelllum_group_estimate_LDA'
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
        %         dataDir = {'MVA_comb_mov', 'MVA_comb_prep'}; %%Save folders for each contrast
        contrastN = length(dataDir);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat = fullfile (suitGroupDir, dataDir(i), 'SPM.mat');  %Adjust directory
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
    case 'cerebellum_map_peaks_to_surface' %maps peaks from random effects analysis to SUIT surface
        %Before running, save .nii maps of all significant peaks from
        %SPM random effects analysis. Label them as cond_phase_region
        %(e.g. int_mov_cerebellum_1, rsa_prep_cerebellum_3).
        %Generate these by opening SPM (SPM fmri) and selecting results.
        %Enter thresholds, then select save cluster. Save all to
        %suit_secondlevel/clusters.
        
        cd(fullfile(suitGroupDir, 'clusters'))
        
        filenames = dir('*.nii'); filenames = {filenames.name};
        
        %Edit below according to peaks
        conMov   = suit_map2surf(filenames(1:2)); %2 con peaks during mov
        conPrep  = suit_map2surf(filenames(3:8)); %6 con peaks during prep
        intMov   = suit_map2surf(filenames(9));   %1 int peak during prod
        rsaCross = suit_map2surf(filenames(10:11)); %2 rsa cross peaks
        rsaPrep  = suit_map2surf(filenames(12)); %1 rsa prep preak
        rsaProd  = suit_map2surf(filenames(13:15)); %3 rsa prod peaks
        
        conMov(isnan(conMov)) = 0;
        conPrep(isnan(conPrep)) = 0;
        intMov(isnan(intMov)) = 0;
        rsaCross(isnan(rsaCross)) = 0;
        rsaPrep(isnan(rsaPrep)) = 0;
        rsaProd(isnan(rsaProd)) = 0;
        
        peakSurfaceDir = fullfile(pwd, 'surface');
        if ~isfolder(peakSurfaceDir)
            mkdir(peakSurfaceDir)
        end
        cd(peakSurfaceDir)
        
        C.cdata=conMov; C=gifti(C);
        save(C,'conMovPeakTCB.func.gii'); %1prep 2prod
        clear C
        
        C.cdata=conPrep; C=gifti(C);
        save(C,'conPrepPeakTCB.func.gii'); %1prep 2prod
        clear C
        
        C.cdata=intMov; C=gifti(C);
        save(C,'intMovPeakTCB.func.gii'); %1prep 2prod
        clear C
        
        C.cdata=rsaCross; C=gifti(C);
        save(C,'rsaCrossPeakTCB.func.gii'); %1prep 2prod
        clear C
        
        C.cdata=rsaPrep; C=gifti(C);
        save(C,'rsaPrepPeakTCB.func.gii'); %1prep 2prod
        clear C
        
        C.cdata=rsaProd; C=gifti(C);
        save(C,'rsaProdPeakTCB.func.gii'); %1prep 2prod
        clear C
        
    case 'cerebellum_plot' %plot area distance and decoding results
        regionNames = suitCBRegions(2,:);%names from second row of variable
        regionNames = strrep(regionNames, '_', ' '); %replace _ with space
        
        %Region names
        regnames = {'lLob4', 'rLob4', 'lLob5', 'rLob5', 'lLob6', 'rLob6', 'lCru1', 'rCru1', 'lCru2', 'rCru2'};
        
        %Text
        titleFontSize = 12;
        fontSize      = 12;
        
        %Colour
        decodeBRG = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880], [0.4660 0.6740 0.1880]};
        rsaGB     = {[0.8 0.8 0.8], [0 0 0]};
        
        %Style
        spl = [zeros(4,1); ones(4,1)];
        label = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        ms = 12;
        ls = 20;
        lw = 3;
        
        %y axis scales
        pscY   =   [-0.5,  1.0];
        rsaY   =   [-0.003, 0.0045];
        ldaY   =   [-4.1,  4.5];
        crossY =   [0,      0.11];
        mdsDistY = [0       0.07];
        eigY     = [-0.0003 0.29];
        
        %%%Load
        dataFile = fullfile(rsaDir, 'cerebellum', 'cbRoiDistances_withmds.mat');
        load(dataFile, 'R', 'Rdist', 'Racc', 'Rperc')%if it exists, load it
        
        %%%Collate data
        Tperc = tapply(Rperc,{'SN', 'region', 'phase'},{'perc', 'mean', 'name', 'perc'});
        Tdist = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset', ...
            Rdist.phase < 3 & Rdist.cond == 1);
        Tacc = tapply(Racc,{'SN', 'region', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'});
        Tcross = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, ...
            'subset', Rdist.phase == 3 & Rdist.cond == 1);
        TmdsDist = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, ...
            'subset', Rdist.cond == 2);
        
        %%%Plot
        
        %         figure %%%General percent signal change overview (zoom to regions of interest)
        %         Tperc = tapply(Rperc,{'SN', 'region', 'phase'},{'perc', 'mean', 'name', 'perc'});
        %         colour={[0 0 0], [1 1 1]};
        %         regions = repmat(regnames, 1, 12);
        %         barplot([Tperc.phase, Tperc.region], Tperc.perc, 'split', Tperc.phase, 'facecolor', colour);
        %         ylim(pscY)
        %         xticklabels(regions)
        %         ylabel('% signal change')
        %         title('Overview - activity increases')
        %         %-------------------------------------------------------------------------------%
        %
        %
        %         figure %%% General distance overview (zoom to regions of interest)
        %         Tdist = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset',...
        %             Rdist.phase < 3 & Rdist.cond == 1);
        %         colour=rsaGB;
        %         regions = repmat(regnames, 1, 12);
        %         barplot([Tdist.phase, Tdist.region], Tdist.dist, 'split', Tdist.phase, 'facecolor', colour);
        %         ylim(rsaY)
        %         xticklabels(regions)
        %         ylabel('Crossnobis dissimilarity')
        %         title('Overview - representational similarity analysis')
        %         %-------------------------------------------------------------------------------%
        %
        %
        %         figure %%% General decoding overview (zoom to regions of interest)
        %         Tacc = tapply(Racc,{'SN', 'region', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'});
        %         colour=decodeBRG;
        %         regions = repmat(regnames, 1, 12);
        %         barplot([Tacc.phase, Tacc.cond, Tacc.region], Tacc.acc, 'split', Tacc.cond, 'facecolor', colour);
        %         drawline(0, 'dir', 'horz', 'linestyle', '- -')
        %         ylabel('Decoding accuracy')
        %         xticklabels(regions)
        %         title('Overview - Linear decoding accuracy')
        %         %-------------------------------------------------------------------------------%
        %
        %
        %         figure %%% General cross-phase distance overview
        %         Tcross = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset',...
        %             Rdist.phase == 3 & Rdist.cond == 1);
        %         colour={[0 0 0]};
        %         regions = repmat(regnames, 1, 12);
        %         barplot(Tcross.region, Tcross.dist, 'facecolor', colour);
        %         ylim(crossY)
        %         xticklabels(regions)
        %         ylabel('Crossnobis dissimilarity')
        %         title('Overview - cross-phase RSA')
        %         %-------------------------------------------------------------------------------%
        %
        %
        %         figure %%% General MDS (PC2 & PC3) distance overview
        %         TmdsDist = tapply(Rdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, ...
        %             'subset', Rdist.cond == 2);
        %         colour={[0 0 0]};
        %         regions = repmat(regnames, 1, 12);
        %         barplot(TmdsDist.region, TmdsDist.dist, 'facecolor', colour);
        %         ylim(mdsDistY)
        %         xticklabels(regions)
        %         ylabel('Euclidean distance')
        %         title('Overview - MDS PC2 & PC3 RSA')
        %         %-------------------------------------------------------------------------------%
        
        
        %         figure %%% PSC Region subplots for prep/prod activity
        %         loopCount = 1;
        %         for i=1:10%plots left hem on the left, right hem on the right
        %             subplot(5,2,loopCount)
        %             loopCount = loopCount + 1;
        %             T = tapply(Rperc,{'SN', 'cond', 'phase'},{'perc', 'mean', 'name', 'perc'}, 'subset',...
        %                 Rperc.region == i);
        %             colour={[0 0 0]}; %black %{[0 0.545 0.545], [1 0.647 0]}; %blue & orange
        %             lineplot([T.phase], T.perc, ...
        %                 'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
        %                 'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
        %             ylim(pscY)
        %             drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
        %             if i==1
        %                 ylabel('% signal change')
        %                 set(gca,'xticklabel',{'Prep', 'Prod'})
        %             else
        %                 ylabel('')
        %                 yticklabels({''})
        %                 set(gca,'xticklabel',{''})
        %             end
        %
        %             set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
        %             title(regionNames{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        %         end%for subcort region
        %         %-------------------------------------------------------------------------------%
        
        
        figure %%% PSC Region subplots for prep/prod activity - box plots
        loopCount = 1;
        for i=1:10%plots left hem on the left, right hem on the right
            subplot(5,2,loopCount)
            loopCount = loopCount + 1;
            T = tapply(Rperc,{'SN', 'cond', 'phase'},{'perc', 'mean', 'name', 'perc'}, 'subset', ...
                Rperc.region == i);
            colour={[0 0 0]}; %black %{[0 0.545 0.545], [1 0.647 0]}; %blue & orange
            myboxplot([T.phase], T.perc, ...
                'linewidth', 2, 'markersize', 5, 'markertype', 'o')%, ...
            %'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
            %'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
            ylim(pscY)
            drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
            if i==1
                ylabel('% signal change')
                set(gca,'xticklabel',{'Prep', 'Prod'})
            else
                ylabel('')
                yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(regionNames{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        end%for subcort region
        %-------------------------------------------------------------------------------%
        
        
%         figure %%% RSA Region subplots for prep/prod distances - box plots
%         loopCount = 1;
%         for i=1:10%plots left hem on the left, right hem on the right
%             subplot(5,2,loopCount)
%             loopCount = loopCount + 1;
%             T = tapply(Rdist,{'SN', 'cond', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset',...
%                 Rdist.region == i & Rdist.phase < 3 & Rdist.cond == 1);
%             colour={[0 0 0]}; %black %{[0 0.545 0.545], [1 0.647 0]}; %blue & orange
%             myboxplot([T.phase], T.dist)%, ...
%             %'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
%             %'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
%             ylim(rsaY)
%             drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
%             if i==1
%                 ylabel('Crossnobis Distance')
%                 set(gca,'xticklabel',{'Prep', 'Prod'})
%             else
%                 ylabel('')
%                 yticklabels({''})
%                 set(gca,'xticklabel',{''})
%             end
%             
%             set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
%             title(regionNames{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
%         end%for subcort region
%         %-------------------------------------------------------------------------------%
        
        figure %%% LDA Overall classifier Region subplots for prep/prod activity - box plots
        loopCount = 1;
        for i=1:10%plots left hem on the left, right hem on the right
            subplot(5,2,loopCount)
            loopCount = loopCount + 1;
            T = tapply(Racc,{'SN', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'}, 'subset', ...
                Racc.region == i & Racc.cond == 1);
            colour={[0 0 0]}; %black %{[0 0.545 0.545], [1 0.647 0]}; %blue & orange
            myboxplot([T.phase], T.acc, ...
                'linewidth', 2, 'markersize', 5, 'markertype', 'o')%, ...
            %'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
            %'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
            ylim([-4.1, 5.7])
            drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
            if i==1
                ylabel('Decoding accuracy (Z)')
                set(gca,'xticklabel',{'Prep', 'Prod'})
            else
                ylabel('')
                yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(regionNames{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        end%for subcort region
        %-------------------------------------------------------------------------------%
        
        
        figure %%% LDA Region subplots for prep/prod order/timing/integrated
        loopCount = 1;
        for i=1:10%plots left hem on the left, right hem on the right
            subplot(5,2,loopCount)
            loopCount = loopCount + 1;
            T = tapply(Racc,{'SN', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'},...
                'subset', Racc.region == i & Racc.cond ~= 5 & Racc.cond ~=1);
            colour=decodeBRG;
            lineplot([T.phase], T.acc, 'split', T.cond, ...
                'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
                'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
            ylim(ldaY)
            drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
            if i==1
                ylabel('Decoding accuracy (Z)')
                set(gca,'xticklabel',{'Prep', 'Prod'})
            else
                ylabel('')
                yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(regionNames{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        end%for subcort region
        %-------------------------------------------------------------------------------%
        
        
        figure %%% LDA Region subplots for prep/prod order/timing/integrated - box plots
        loopCount = 1;
        for i=1:10%plots left hem on the left, right hem on the right
            subplot(5,2,loopCount)
            loopCount = loopCount + 1;
            T = tapply(Racc,{'SN', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'}, ...
                'subset', Racc.region == i & Racc.cond ~= 5 & Racc.cond ~=1);
            colour=decodeBRG;
            myboxplot([T.phase], T.acc, 'split', T.cond, 'fillcolor', colour)%...
            %'markertype', 'o', 'markercolor', colour, 'markerfill', colou%, 'markersize', 5, ...
            %'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
            ylim(ldaY)
            drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
            if i==1
                ylabel('Decoding accuracy (Z)')
                set(gca,'xticklabel',{'Prep', '', '', 'Prod', '', ''})
            else
                ylabel('')
                %yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(regionNames{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        end%for subcort region
        %-------------------------------------------------------------------------------%
        
        
        %         figure %%% LDA Region subplots for prep/prod order/timing/integrated - box plots including new and old integrated
        %         loopCount = 1;
        %         for i=[1 5 2 6 3 7 4 8]%plots left hem on the left, right hem on the right
        %             subplot(4,2,loopCount)
        %             loopCount = loopCount + 1;
        %             T = tapply(Racc,{'SN', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'}, 'subset', Racc.region == i);
        %             colour=decodeBRG;
        %             myboxplot([T.phase], T.acc, 'split', T.cond, 'fillcolor', colour)%...
        %             %'markertype', 'o', 'markercolor', colour, 'markerfill', colou%, 'markersize', 5, ...
        %             %'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour);
        %             ylim(ldaY)
        %             drawline(0, 'dir', 'horz', 'linestyle', '-', 'linewidth', 1)
        %             if i==1
        %                 ylabel('Decoding accuracy (Z)')
        %                 set(gca,'xticklabel',{'Prep', '', 'Old', 'New', 'Prod', '', 'Old', 'New'})
        %             else
        %                 ylabel('')
        %                 %yticklabels({''})
        %                 set(gca,'xticklabel',{''})
        %             end
        %             set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
        %             title(regionNames{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        %         end%for subcort region
        %         %-------------------------------------------------------------------------------%
        
        
        %%% RDMs and multi-dimensional scaling plots (for visualisation)
        labels = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        
        %Extract data
        for s=1:length(R.G)%for subj * region
            G(:,:,s)= R.G{s}; %extract representational dissimilarity matrix into 3D matrix
        end%for subj * region
        
        figure %Plot RDMs
        loopCount = 1;
        for i=1:10%plots left hem on the left, right hem on the right
            
            GRegion = G(:,:,R.region == loopCount); %extract region variance/covariance matrices
            GmRegion = mean(GRegion, 3); %mean across subjs
            
            subplot(5, 2, i) %%%RDM for each subcortical region
            ind=indicatorMatrix('allpairs',1:8); %matrix for all pairwise distances (k*(k-1))
            imagesc(rsa.rdm.squareRDM(diag(ind*GmRegion*ind'))); %display
            %multiplying variance/covariance by indicator matrix results in
            %dissimilarity values (crossnobis)
            title([regionNames{loopCount} ' RDM (crossnobis)'])
            colorbar
            axis image
            loopCount = loopCount + 1;
        end
        %----------------------------------------------------------------------------------------------%
        
        
        
        loopCount = 1; %Multi-dimensional scaling
        for i=1:10%plots left hem on the left, right hem on the right
            
            figure %plot MDS
            
            GRegion = G(:,:,R.region == loopCount); %extract region variance/covariance matrices
            GmRegion = mean(GRegion, 3); %mean across subjs
            
            [COORD,~]=pcm_classicalMDS(GmRegion);
            
            %3D scatter plot
            scatterplot3(COORD(1:end,2),COORD(1:end,3),COORD(1:end,1),'split',spl, ... %here we plot PC 2&3 first because
                'markersize',ms, 'markercolor',rsaGB, 'markerfill',rsaGB, 'label',label);%PC 1 is just activity differences
            
            %%%Draw coloured lines between distinct order & timing conditions
            %Prep
            colors = decodeBRG{1}; %blue for order
            indx=[1 3]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[2 4]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            colors = decodeBRG{2}; %red for timing
            indx=[1 2]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[3 4]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            
            %Prod
            colors = decodeBRG{1}; %blue for order
            indx=[5 7]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[6 8]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            colors = decodeBRG{2}; %red for timing
            indx=[5 6]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[7 8]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            
            grid off
            hold on; plot3(0,0,0,'+','MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',ms+3, 'LineWidth',lw);
            hold off; xlabel('PC 2'); ylabel('PC 3'); zlabel('PC 1'); set(gca,'fontsize',12);
            %axis equal
            title([regionNames{loopCount} ' multi-dimensional scaling'])
            %xlim([-0.035, 0.035])
            %ylim([-0.027, 0.029])
            %zlim([-0.02, 0.3])
            view(30,30)
            loopCount = loopCount + 1;
        end
        %----------------------------------------------------------------------------------------------%
        
        
        %%%Plot Eigenvalues for PC1, PC2, and PC3
        eigen.vals   = [];%preallocate struct vars
        eigen.pcs    = [];
        eigen.region = [];
        eigen.sn     = [];
        
        figure
        for i = 1:10%plots left hem on the left, right hem on the right %unique(Rdist.region)'
            eigValsAll = [Rdist.eigVals{Rdist.region == i & Rdist.cond == 2}]';
            eigVals    = reshape(eigValsAll(:,1:3), 72,1);%takes the first 3 PCs from eigValsAll and makes them one column
            eigPcs     = [ones(24,1); ones(24,1)*2; ones(24,1)*3];%variable to track which PC
            eigRegion  = ones(72,1)*i;%track region num
            eigSn      = repmat((1:24)', 3, 1);
            
            eigen.vals   = [eigen.vals; eigVals];%concat struct
            eigen.pcs    = [eigen.pcs; eigPcs];
            eigen.region = [eigen.region; eigRegion];
            eigen.sn     = [eigen.sn; eigSn];
            
            subplot(5,2,i)
            myboxplot(eigen.pcs, eigen.vals, 'subset', eigen.region == i, 'xtickoff')
            ylim(eigY)
            if i==1
                ylabel('Eigenvalues (a.u.)')
                set(gca,'xticklabel',{'PC1', 'PC2', 'PC3'})
            else
                ylabel('')
                %yticklabels({''})
                set(gca,'xticklabel',{''})
            end
            set(gca,'FontSize',fontSize, 'FontName', 'Calibri')
            title(regionNames{i}, 'FontSize', titleFontSize, 'FontName', 'Calibri')
        end
        %Save eigen values to excel
        varnames = {...
            'subj', ...
            'PC1_LLob4', 'PC2_LLob4', 'PC3_LLob4', ...
            'PC1_RLob4', 'PC2_RLob4', 'PC3_RLob4', ...
            'PC1_LLob5', 'PC2_LLob5', 'PC3_LLob5', ...
            'PC1_RLob5', 'PC2_RLob5', 'PC3_RLob5', ...
            'PC1_LLob6', 'PC2_LLob6', 'PC3_LLob6', ...
            'PC1_RLob6', 'PC2_RLob6', 'PC3_RLob6', ...
            'PC1_LCru1', 'PC2_LCru1', 'PC3_LCru1', ...
            'PC1_RCru1', 'PC2_RCru1', 'PC3_RCru1', ...
            'PC1_LCru2', 'PC2_LCru2', 'PC3_LCru2', ...
            'PC1_RCru2', 'PC2_RCru2', 'PC3_RCru2', ...
            };
        
        eigenData = [];
        for i = unique(eigen.region)'
            for j = unique(eigen.pcs)'
                eigenData(:,end+1) = eigen.vals(eigen.region == i & eigen.pcs == j);
            end
        end
        
        eigenTable = array2table([anaSubj' eigenData], 'VariableNames', varnames);
        writetable(eigenTable, fullfile(roiCbDir, 'eigenValsROI.xlsx'))
        
        %%%Save all results to respective, formatted excel files for SPSS stats
        %psc
        varnames = {...
            'subj', ...
            'prepLLob4', 'prepRLob4', 'prepLLob5', 'prepRLob5', 'prepLLob6', 'prepRLob6', 'prepLCru1', 'prepRCru1', 'prepLCru2', 'prepRCru2' ...
            'prodLLob4', 'prodRLob4', 'prodLLob5', 'prodRLob5', 'prodLLob6', 'prodRLob6', 'prodLCru1', 'prodRCru1', 'prodLCru2', 'prodRCru2' ...
            };
        percData = [];
        for i=unique(Tperc.phase)'
            for j=unique(Tperc.region)'
                percData(:,end+1) = Tperc.perc(Tperc.phase == i & Tperc.region == j);
            end
        end
        percTable = array2table([anaSubj' percData], 'VariableNames', varnames);
        writetable(percTable, fullfile(roiCbDir, 'percROIspss.xlsx'))
        
        %dist
        varnames = {...
            'subj', ...
            'prepLLob4', 'prepRLob4', 'prepLLob5', 'prepRLob5', 'prepLLob6', 'prepRLob6', 'prepLCru1', 'prepRCru1', 'prepLCru2', 'prepRCru2' ...
            'prodLLob4', 'prodRLob4', 'prodLLob5', 'prodRLob5', 'prodLLob6', 'prodRLob6', 'prodLCru1', 'prodRCru1', 'prodLCru2', 'prodRCru2' ...
            };
        distData = [];
        for i=unique(Tdist.phase)'
            for j=unique(Tdist.region)'
                distData(:,end+1) = Tdist.dist(Tdist.phase == i & Tdist.region == j);
            end
        end
        distTable = array2table([anaSubj' distData], 'VariableNames', varnames);
        writetable(distTable, fullfile(roiCbDir, 'distROIspss.xlsx'))
        
        %acc
        varnames = {...
            'subj', ...
            'ovrPrepLLob4', 'ovrPrepRLob4', 'ovrPrepLLob5', 'ovrPrepRLob5', 'ovrPrepLLob6', 'ovrPrepRLob6', 'ovrPrepLCru1', 'ovrPrepRCru1', 'ovrPrepLCru2', 'ovrPrepRCru2' ...
            'ovrProdLLob4', 'ovrProdRLob4', 'ovrProdLLob5', 'ovrProdRLob5', 'ovrProdLLob6', 'ovrProdRLob6', 'ovrProdLCru1', 'ovrProdRCru1', 'ovrProdLCru2', 'ovrProdRCru2' ...
            'ordPrepLLob4', 'ordPrepRLob4', 'ordPrepLLob5', 'ordPrepRLob5', 'ordPrepLLob6', 'ordPrepRLob6', 'ordPrepLCru1', 'ordPrepRCru1', 'ordPrepLCru2', 'ordPrepRCru2' ...
            'ordProdLLob4', 'ordProdRLob4', 'ordProdLLob5', 'ordProdRLob5', 'ordProdLLob6', 'ordProdRLob6', 'ordProdLCru1', 'ordProdRCru1', 'ordProdLCru2', 'ordProdRCru2' ...
            'tempPrepLLob4', 'tempPrepRLob4', 'tempPrepLLob5', 'tempPrepRLob5', 'tempPrepLLob6', 'tempPrepRLob6', 'tempPrepLCru1', 'tempPrepRCru1', 'tempPrepLCru2', 'tempPrepRCru2' ...
            'tempProdLLob4', 'tempProdRLob4', 'tempProdLLob5', 'tempProdRLob5', 'tempProdLLob6', 'tempProdRLob6', 'tempProdLCru1', 'tempProdRCru1', 'tempProdLCru2', 'tempProdRCru2' ...
            'intPrepLLob4', 'intPrepRLob4', 'intPrepLLob5', 'intPrepRLob5', 'intPrepLLob6', 'intPrepRLob6', 'intPrepLCru1', 'intPrepRCru1', 'intPrepLCru2', 'intPrepRCru2' ...
            'intProdLLob4', 'intProdRLob4', 'intProdLLob5', 'intProdRLob5', 'intProdLLob6', 'intProdRLob6', 'intProdLCru1', 'intProdRCru1', 'intProdLCru2', 'intProdRCru2' ...
            'intPrepOldLLob4', 'intPrepOldRLob4', 'intPrepOldLLob5', 'intPrepOldRLob5', 'intPrepOldLLob6', 'intPrepOldRLob6', 'intPrepOldLCru1', 'intPrepOldRCru1', 'intPrepOldLCru2', 'intPrepOldRCru2' ...
            'intProdOldLLob4', 'intProdOldRLob4', 'intProdOldLLob5', 'intProdOldRLob5', 'intProdOldLLob6', 'intProdOldRLob6', 'intProdOldLCru1', 'intProdOldRCru1', 'intProdOldLCru2', 'intProdOldRCru2' ...
            };
        accData = [];
        for i=unique(Tacc.cond)'
            for j=unique(Tacc.phase)'
                for k=unique(Tacc.region)'
                    accData(:,end+1) = Tacc.acc(Tacc.cond == i & Tacc.phase == j & Tacc.region == k);
                end
            end
        end
        accTable = array2table([anaSubj' accData], 'VariableNames', varnames);
        writetable(accTable, fullfile(roiCbDir, 'accROIspss.xlsx'))
        
        %cross
        varnames = {...
            'subj', ...
            'LLob4', 'RLob4', 'LLob5', 'RLob5', 'LLob6', 'RLob6', 'LCru1', 'RCru1', 'LCru2', 'RCru2' ...
            };
        crossData = [];
        for i=unique(Tcross.region)'
            crossData(:,end+1) = Tcross.dist(Tcross.region == i);
        end
        crossTable = array2table([anaSubj' crossData], 'VariableNames', varnames);
        writetable(crossTable, fullfile(roiCbDir, 'crossROIspss.xlsx'))
        
        %MDS Dist
        varnames = {...
            'subj', ...
            'LLob4', 'RLob4', 'LLob5', 'RLob5', 'LLob6', 'RLob6', 'LCru1', 'RCru1', 'LCru2', 'RCru2' ...
            };
        mdsDistData = [];
        for i=unique(TmdsDist.region)'
            mdsDistData(:,end+1) = TmdsDist.dist(Tcross.region == i);
        end
        mdsDistTable = array2table([anaSubj' mdsDistData], 'VariableNames', varnames);
        writetable(mdsDistTable, fullfile(roiCbDir, 'mdsDistROIspss.xlsx'))
    case 'cerebellum_plot_flatmap_avg' %plot average maps onto cerebellar flat map
        %%%Due to compatibility issues, the flatmap functions tend to error
        %%%when the RSA toolbox is on the path. Make sure you remove the
        %%%RSA functions (specifically gifti) from the path before running.
        
        thresholdLDA = 0.5; %minimum threshold for decoding
        
        
        
        %%%Lobule atlas
        G = gifti(fullfile(suitDir, 'atlasSUIT', 'Lobules.label.gii'));
        %Extract the color map, ignoring the first entry
        CMAP = G.labels.rgba(2:end,1:3);
        %Plot the flatmap. Because the data are discrete category labels, 'type' is set to 'label'.
        figure
        suit_plotflatmap(G.cdata,'type','label','cmap',CMAP,'border',[]);
        title('Lobule atlas')
        
        %%%Buckner 17 Networks atlas
        G = gifti(fullfile(suitDir, 'atlasSUIT', 'Buckner_17Networks.label.gii'));
        %Extract the color map, ignoring the first entry
        CMAP = G.labels.rgba(2:end,1:3);
        %Plot the flatmap. Because the data are discrete category labels, 'type' is set to 'label'.
        figure
        suit_plotflatmap(G.cdata,'type','label','cmap',CMAP,'border',[]);
        title('Buckner 17 Networks atlas')
        
        
        %%%PSC
        cd(fullfile(suitGroupDir, 'average'))
        conds = {...
            'perc_prep',    'perc_mov', ...
            };
        fileNames = strcat('avg_', conds, '.nii');
        surfMap = suit_map2surf(fileNames);
        
        C.cdata=surfMap;
        C=gifti(C);
        save(C,'percAvgCB.func.gii'); %1prep 2prod
        clear C
        
        %--- Prep ---%
        figure
        suit_plotflatmap(surfMap(:,1), 'threshold', 0, 'cscale', [0, 0.2])
        title('% Signal change - preparation')
        
        %--- Prod ---%
        figure
        suit_plotflatmap(surfMap(:,2), 'threshold', 0, 'cscale', [0.0, 0.2])
        title('% Signal change - production')
        
        
        %%%RSA
        cd(fullfile(suitGroupRSADir, 'average'))
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        fileNames = strcat('avg_', conds, '_LDC.nii');
        surfMap = suit_map2surf(fileNames);
        
        C.cdata=surfMap;
        C=gifti(C);
        save(C,'rsaAvgCB.func.gii');
        clear C
        
        %--- Prep ---%
        figure
        suit_plotflatmap(surfMap(:,1), 'threshold', 0.0001, 'cscale', [0.0001, 0.003])
        title('RSA - Overall during preparation')
        
        %--- Prod ---%
        figure
        suit_plotflatmap(surfMap(:,2), 'threshold', 0.0001, 'cscale', [0.0001, 0.003])
        title('RSA - Overall during production')
        
        %--- Switch ---%
        figure
        suit_plotflatmap(surfMap(:,3), 'threshold', 0.0001, 'cscale', [0.0001, 0.003])
        title('RSA - Within sequences across phases')
        % --------------------------------------------------------------- %
        clear fileName
        
        
        %%%LDA
        cd(fullfile(suitGroupDir, 'average'))
        fileName = {...
            'avg_spat_prep_LDA.nii'...
            'avg_temp_prep_LDA.nii'...
            'avg_int_prep_LDA.nii'...
            };
        for i=1:length(fileName)
            surfMapPrep(:,i) = suit_map2surf(fileName{i});
        end
        
        %%%Threshold to z accuracy > 0.5
        surfMapPrepThresh = [surfMapPrep(:,2), surfMapPrep(:,3), surfMapPrep(:,1)]; %R(temp) G(int) B(ord)
        %surfMapPrepThresh(surfMapPrep < thresholdLDA) = NaN;
        
        figure
        suit_plotflatmap(surfMapPrepThresh, 'type', 'rgb', 'alpha', 1)%, 'threshold', 0.01, 'cscale', [0.5, 1])
        % --------------------------------------------------------------- %
        
        %%%Production LDA
        fileName = {...
            'avg_spat_mov_LDA.nii'...
            'avg_temp_mov_LDA.nii'...
            'avg_int_mov_LDA.nii'...
            };
        for i=1:length(fileName)
            surfMapProd(:,i) = suit_map2surf(fileName{i});
        end
        
        C.cdata=[surfMapPrep surfMapProd];
        C=gifti(C);
        save(C,'ldaAvgCB.func.gii');
        
        %%%Threshold to z accuracy > 0.5
        surfMapProdThresh = [surfMapProd(:,2), surfMapProd(:,3), surfMapProd(:,1)]; %R(temp) G(int) B(ord)
        surfMapProdThresh(surfMapProd < thresholdLDA) = NaN;
        
        figure
        suit_plotflatmap(surfMapProdThresh, 'type', 'rgb', 'alpha', 1)%, 'threshold', 0.01, 'cscale', [0.5, 1])
    case 'cerebellum_plot_flatmap_T'   %plot T-maps onto cerebellar flat map
        
        cd(suitGroupDir)
        
        %         conds = {...
        %             'overall_prep',    'overall_prod',    'overall_cross', ...
        %             };
        %         folderNames = strcat('RSA_', conds);
        %
        %         fileName = fullfile(folderNames{4}, 'spmT_0001.nii');
        %         surfMap = suit_map2surf(fileName);
        %         figure
        %         suit_plotflatmap(surfMap, 'threshold', 3.48, 'cscale', [0, 3.5])
        %
        %         clear fileName
        
        fileName{1} = fullfile('MVA_temp_prep', 'spmT_0001.nii');
        fileName{2} = fullfile('MVA_int_prep', 'spmT_0001.nii');
        fileName{3} = fullfile('MVA_spat_prep', 'spmT_0001.nii');
        surfMap(:,1) = suit_map2surf(fileName{1});
        surfMap(:,2) = suit_map2surf(fileName{2});
        surfMap(:,3) = suit_map2surf(fileName{3});
        
        figure
        suit_plotflatmap(surfMap, 'type', 'rgb', 'threshold', 10, 'cscale', [0, 3.5])
        
        fileName{1} = fullfile('MVA_temp_mov', 'spmT_0001.nii');
        fileName{2} = fullfile('MVA_int_mov', 'spmT_0001.nii');
        fileName{3} = fullfile('MVA_spat_mov', 'spmT_0001.nii');
        surfMap(:,1) = suit_map2surf(fileName{1});
        surfMap(:,2) = suit_map2surf(fileName{2});
        surfMap(:,3) = suit_map2surf(fileName{3});
        
        figure
        suit_plotflatmap(surfMap, 'type', 'rgb', 'threshold', 10, 'cscale', [0, 3.5])
        
    case 'subAndCB_plot' %plot subcortical and CB effector regions together - distances only (cross & MDS)
        %Also requires case subcortical&cerebellum_mds_distances to run
        
        %%%Design variables
        %Text
        titleFontSize = 12;
        fontSize      = 12;
        
        %Colour
        decodeBRG = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
        rsaGB     = {[0.8 0.8 0.8], [0 0 0]};
        
        %Style
        spl = [zeros(4,1); ones(4,1)];
        label = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        ms = 12;
        ls = 20;
        lw = 3;
        
        %y axis scales
        crossY   = [-0.04,      0.2];
        mdsDistY = [0       0.075];
        
        %region names
        subcortNames = {'lTha', 'lCau', 'lPut', 'lHip', 'rTha', 'rCau', 'rPut', 'rHip'};
        CBnames      = {'lLob4', 'rLob4', 'lLob5', 'rLob5', 'lLob6', 'rLob6', 'lCru1', 'rCru1', 'lCru2', 'rCru2'};
        
        targetSubNames = {'lTha', 'lCau', 'lPut', 'lHip', 'rHip'};
        targetSubNums  = [ 1       2       3       4       8    ];
        
        targetCBNames  = {'rLob4', 'rLob5'};
        targetCBNums   = [ 2        4     ];
        
        %%%Subcortical ----------------------------------------------------
        subcortStructs = subcortStructs(2, :);%names from second row
        subcortStructs = strrep(subcortStructs, '_', ' '); %replace _ with space
        
        %%%Load
        dataFile = fullfile(rsaDir, 'subcortical', 'subRoiDistances_withmds.mat');
        load(dataFile, 'Rdist');%load it
        
        mask = ismember(Rdist.region, targetSubNums); %remove non-target regions
        names = fieldnames(Rdist);
        for i=1:numel(names)
            Rdist.(names{i})(~mask) = [];
        end
        Rdist.region(Rdist.region == 8) = 5; %change rHip to region 5 for plotting
        subcortdist = Rdist;
        
        
        %%%Cerebellum -----------------------------------------------------
        regionNames = suitCBRegions(2,:);%names from second row of variable
        regionNames = strrep(regionNames, '_', ' '); %replace _ with space
        
        %%%Load
        dataFile = fullfile(rsaDir, 'cerebellum', 'cbRoiDistances_withmds.mat');
        load(dataFile, 'Rdist');%load it
        CBdist = Rdist;
        
        mask = ismember(Rdist.region, targetCBNums); %remove non-target regions
        names = fieldnames(Rdist);
        for i=1:numel(names)
            Rdist.(names{i})(~mask) = [];
        end
        Rdist.region(Rdist.region == 2) = 6; %change rLob4 to region 6 for plotting
        Rdist.region(Rdist.region == 4) = 7; %change rLob5 to region 7 for plotting
        CBdist = Rdist;
        
        %%%New struct with only target regions
        tgtdist = addstruct(subcortdist, CBdist);
        
        
        
        %%%Plot
        figure %%% General cross-phase distance overview
        Tcross = tapply(tgtdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, ...
            'subset', tgtdist.phase == 3 & tgtdist.cond == 1);
        colour={[0 0 0]};
        regions = [targetSubNames, targetCBNames];%repmat({'lTha', 'lCau', 'lPut', 'lHip', 'rTha', 'rCau', 'rPut', 'rHip'}, 1, 12);
        myboxplot(Tcross.region, Tcross.dist, 'whiskerwidth', 2.5, 'xtickoff');
        ylim(crossY)
        xticklabels(regions)
        ylabel('Crossnobis dissimilarity')
        title('Cross-phase RSA')
        %-------------------------------------------------------------------------------%
        
        
        figure %%% General MDS (PC2 & PC3) distance overview
        TmdsDist = tapply(tgtdist,{'SN', 'region', 'phase'},{'dist', 'mean', 'name', 'dist'}, ...
            'subset', tgtdist.cond == 2);
        colour={[0 0 0]};
        regions = [targetSubNames, targetCBNames];
        myboxplot(TmdsDist.region, TmdsDist.dist, 'whiskerwidth', 2.5, 'xtickoff');
        ylim(mdsDistY)
        xticklabels(regions)
        ylabel('Euclidean distance')
        title('MDS PC2 & PC3 RSA')
        %-------------------------------------------------------------------------------%
        
    case 'cortical_make_search'
        
        s=varargin{1};
        
        radius=16;
        numVox=160;
        
        cd(glmDir);
        V=spm_vol(fullfile(subj_name{s}, 'maskbrain.nii')); %preceded by case suit_make_mask
        X=spm_read_vols(V);
        [i,j,k]=ind2sub(size(X),find(X~=0));
        vox=[i j k];
        
        L = [];
        [L.LI,L.voxmin,L.voxmax,L.n]=lmva_voxelselection(vox(:,:)',vox',[radius numVox],V.mat,V.dim,[], fullfile(subj_name{s}, 'mva160_numvoxRSA.nii'));
        L.voxel = vox;
        save(fullfile(subj_name{s}, 'volsearch160RSA.mat'), 'L')
    case 'cortical_run_search'
        
        s=varargin{1}; blueBear=varargin{2};
        nrruns = length(run);
        
        rsaSubjFolder = fullfile(rsaDir, 'cortical', subj_name{s});
        if ~isfolder(rsaSubjFolder)
            mkdir(rsaSubjFolder)
        end
        cd(rsaSubjFolder)
        
        %%% Searchlight file
        load(fullfile(glmDir,subj_name{s},'volsearch160RSA.mat'), 'L');
        
        %%% SPM file
        spmDir = fullfile(glmDir, subj_name{s});
        load(fullfile(spmDir, 'SPM'), 'SPM');
        
        %replace fnames for bluebear compatibility
        if blueBear == 1
            for i=1:length(SPM.xY.VY)
                SPM.xY.VY(i).fname = strrep(SPM.xY.VY(i).fname,'\','/');
                SPM.xY.VY(i).fname = strrep(SPM.xY.VY(i).fname,'Z:','/rds/projects/k/kornyshk-kornyshevalab');
            end
        else
        end
        
        %%% prepare condition and partition (run) vectors
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 0 0]; %prep
        prep=[repmat(prep,1,nrruns) runBSL];%1 x nBeta, 1 2 3 4 = prep sequences 1:4
        
        prod      =[5 0 0 0 6 0 0 0 7 0 0 0 8 0 0 0 0 0 0 0]; %prod
        prod=[repmat(prod,1,nrruns) runBSL];%1 x nBeta, 5 6 7 8 = prod sequences 1:4
        
        condVec = prep + prod; condVec = condVec'; % conditions, including no interest regressors as 0
        
        %%% Run searchlight function on whole CB
        rsa.runSearchlightLDC_RY(L, SPM, 'spmDir', spmDir, 'conditionVec', condVec, ...
            'analysisName', [subj_name{s}, 'RSA_All'], 'outDir', rsaSubjFolder)
    case 'cortical_smooth'
        
        s=varargin{1};
        rsaCorticalDir = fullfile(rsaDir, 'cortical');
        
        comb=fullfile(rsaCorticalDir, subj_name{s},[subj_name{s} 'RSA_All_LDC.nii']); %%MVPA smoother
        scomb=fullfile(rsaCorticalDir, subj_name{s},[subj_name{s} 'RSA_ALL_sLDC.nii']);
        spm_smooth(comb,scomb,[2 2 2]); %smooth with 2mm kernel
        
        %comb=fullfile(rsaCorticalDir, subj_name{s},[subj_name{s} 'RSA_All_integrated_LDC.nii']); %%MVPA smoother
        %scomb=fullfile(rsaCorticalDir, subj_name{s},[subj_name{s} 'RSA_ALL_integrated_sLDC.nii']);
        %spm_smooth(comb,scomb,[2 2 2]); %smooth with 2mm kernel
    case 'cortical_calc_dissimilarity_maps' %extract each condition (overall, order, timing, etc) from RSA_ALL_sLDC.nii
        %each volume of the searchlight corresponds to a pairwise
        %dissimilarity measure between sequences. So here we average within
        %our conditions to give us dissimilarity measures for:
        %overall prep, overall prod, overall cross, order prep, order prod,
        %order cross, timing prep, timing prod, timing cross
        
        s=varargin{1};
        cd(fullfile(rsaCorticalDir, subj_name{s}))
        
        %%% Identify and extract values for overall, order, and timing
        %%% within preparation, production, and cross-phase
        %ordDiff        = [0 1 1 0 0 1 1 1 1 0 0 1 1 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0]';   %pairwise contrast index (1s are where orders are different)
        %timDiff        = [1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1]';   %^ index (1s are where timings are different)
        prepCols       = [1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';   %^ index (1s are contrasts within preparation)
        prodCols       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1]';   %^ index (1s are contrasts within production)
        crossPhaseCols = [0 0 0 1 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0]';   %^ index (1s are contrasts across phases)
        
        conds = {...
            'overall_prep', 'overall_prod', 'overall_cross', ...
            prepCols == 1,  prodCols == 1,  crossPhaseCols == 1 ...
            };
        %'order_prep',                                'order_prod',                                'order_cross', ...
        %'timing_prep',                               'timing_prod',                               'timing_cross';...
        %prepCols == 1 & ordDiff == 1 & timDiff == 0, prodCols == 1 & ordDiff == 1 & timDiff == 0, crossPhaseCols == 1 & ordDiff == 1 & timDiff == 0, ...
        %prepCols == 1 & timDiff == 1 & ordDiff == 0, prodCols == 1 & timDiff == 1 & ordDiff == 0, crossPhaseCols == 1 & timDiff == 1 & ordDiff == 0 ...
        for j=1:length(conds)
            vol = spm_vol([subj_name{s} 'RSA_ALL_sLDC.nii']);
            
            Vi = vol(conds{2,j}, :);
            Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
            Vo.fname = [subj_name{s} '_LDC_' conds{1,j} '.nii'];
            Vo.n = [1 1];
            express = 'mean(X)';
            flags.dmtx = 1;
            
            spm_imcalc(Vi, Vo, express, flags)
        end
    case 'cortical_normalise'
        sn=varargin{1};
        cd([rsaDir '/cortical/']);
        
        if isfolder(rsaCorticalGroupDir) == 0
            mkdir(rsaCorticalGroupDir)
        end
        disp(['normalising ' subj_name{sn}])
        
        %MVPA accuracy maps
        images= {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            }; % please add other images as required, e.g. spmT_...
        
        %'order_prep',      'order_prod',      'order_cross', ...
        %'timing_prep',     'timing_prod',     'timing_cross'...
        %'integrated_prep', 'integrated_prod', 'integrated_cross'...
        
        
        s=varargin{1};
        defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
        for j=1:numel(images)
            [~,name,~]=spm_fileparts(images{j});
            in_images{1}= fullfile(rsaCorticalDir,subj_name{s},[subj_name{s} '_LDC_' images{j} '.nii']);
            out_images{1}= fullfile(rsaCorticalGroupDir,[name '_' subj_name{s} '.nii']);
            spmj_normalization_write(defor, in_images,'outimages',out_images); %Trilinear interpolation
        end
    case 'cortical_group_avg'
        
        if ~isfolder(fullfile(rsaCorticalGroupDir, 'average'))
            mkdir(fullfile(rsaCorticalGroupDir, 'average'))
        end
        cd(fullfile(rsaCorticalGroupDir, 'average'))
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        
        %'order_prep',      'order_prod',      'order_cross', ...
        %'timing_prep',     'timing_prod',     'timing_cross'...
        %'integrated_prep', 'integrated_prod', 'integrated_cross', ...
        
        for i=1:length(conds)
            loopCount = 1;
            for s = anaSubj
                Vi(loopCount) = spm_vol(fullfile(rsaCorticalGroupDir, [conds{i} '_' subj_name{s} '.nii']));
                loopCount = loopCount + 1;
            end
            Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
            Vo.fname = ['avg_' conds{i} '_LDC.nii'];
            Vo.n = [1 1];
            express = 'mean(X)';
            flags.dmtx = 1;
            
            spm_imcalc(Vi, Vo, express, flags)
        end
    case 'cortical_group_randomeffects'
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        
        %'order_prep',      'order_prod',      'order_cross', ...
        %'timing_prep',     'timing_prod',     'timing_cross'...
        %'integrated_prep', 'integrated_prod', 'integrated_cross', ...
        
        dataDir = strcat('RSA_', conds);
        
        contrastN = length(dataDir);
        images = repmat(conds', 1, length(anaSubj));
        subNii = repmat (subj_name(anaSubj), length(dataDir), 1);
        fileName = strcat (images, '_', subNii, '.nii');  %%Concatenate contrast files and subject names
        
        for i=1:contrastN  %%Loop across contrasts
            glmscndDir = fullfile(rsaCorticalGroupDir, dataDir(i));
            matlabbatch{1}.spm.stats.factorial_design.dir = glmscndDir;  %Adjust directory
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fullfile (rsaCorticalGroupDir, fileName(i,:))';  %%Select files from matrix
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
    case 'cortical_group_estimate'
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            };
        
        %'order_prep',      'order_prod',      'order_cross', ...
        %'timing_prep',     'timing_prod',     'timing_cross'...
        %'integrated_prep', 'integrated_prod', 'integrated_cross', ...
        
        contrastN = length(conds);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile (rsaCorticalGroupDir, ['RSA_' conds{i}], 'SPM.mat');  %Adjust directory
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
        
    case 'simulations_RSA_LDA' %simulates data according to inputs and analyses
        
        %%%Prepare variables
        prepProd = 1;
        vord     = [0.3 0.1];
        vtemp    = [0.4 0.7];
        vinter   = [0.6 0.8];
        vnoise   = 0.5;
        nIter    = 500;%iterations of simulation
        
        %Text
        titleFontSize = 12;
        fontSize      = 12;
        
        %Colour
        decodeBRG = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
        rsaGB     = {[0.8 0.8 0.8], [0 0 0]};
        
        %Style
        spl = [zeros(4,1); ones(4,1)];
        label = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        ms = 12;
        ls = 20;
        lw = 3;
        
        %y axis scales
        pscY     = [-0.12,  0.7];
        rsaY     = [-0.003, 0.0045];
        ldaY     = [-1.52,  2];
        crossY   = [0,      0.11];
        mdsDistY = [0       0.07];
        
        vararginoptions(varargin,{'prepProd','prodscaling','prewh','vord','vtemp','vinter','vnoise','nIter'});
        A=[]; %data struct
        B=[]; %output struct - stores distances, LDA, G matrix
        D=[]; %raw data store
        
        %%%Prepare condition and partition (run) vectors
        nrruns = 6; nCond = 8;
        runC      = repmat(1:nrruns,nCond,1); %labelling runs 1:6 (imaging runs)
        runC      = reshape(runC,1,[])';
        condC     = repmat((1:nCond)', nrruns, 1); %labelling conditions 1:4 (sequences)
        
        %%%Generate data
        for i=1:nIter %for 'subject' n
            %make the data based on input specifications
            [Y, Y_hat, ~, snr]=prepProdSimu_makedataPP('prepProd',prepProd,'prodscaling',prodscaling,'vord',vord,'vtemp',vtemp,'vinter',vinter,'vnoise',vnoise);
            
            A.data{i}  = Y;
            AP.data{i} = Y_hat;
            A.snr{i} = snr;
            
            %%%Calculate second moment matrix
            [A.G{i},~]    = pcm_estGCrossval(A.data{i},runC,condC);
            [AP.G{i},~]   = pcm_estGCrossval(AP.data{i},runC,condC);
            
            %%%Generate distances vector
            [A.d{i},     A.Sig{i}]     = rsa.distanceLDC(A.data{i}, runC, condC);
            [AP.d{i},     AP.Sig{i}]     = rsa.distanceLDC(AP.data{i}, runC, condC);
            
            %%%Run linear decoding
            A.prepData = [A.data{i}(1:4,:); A.data{i}(9:12,:);  A.data{i}(17:20,:); A.data{i}(25:28,:); A.data{i}(33:36,:); A.data{i}(41:44,:)];
            A.prodData = [A.data{i}(5:8,:); A.data{i}(13:16,:); A.data{i}(21:24,:); A.data{i}(29:32,:); A.data{i}(37:40,:); A.data{i}(45:48,:)];
            
            %               accuracy columns
            %1: order,          2: timing,         3:integrated
            [A.prepAccuracy(i,1), A.prepAccuracy(i,2), A.prepAccuracy(i,3)] = prepProdSimu_classify(A.prepData);
            [A.prodAccuracy(i,1), A.prodAccuracy(i,2), A.prodAccuracy(i,3)] = prepProdSimu_classify(A.prodData);
        end%for subject n
        
        
        phaseVector = [...%index (1s are where phase is different)
            1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...%prep
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1; ...%prod
            0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0  ...%cross
            ];
        
        %%%Loop over subjs, conditions, and phases, and store distances as
        %%%struct to plot later
        loopCounter = 1;
        
        %%% Distance
        for s=1:nIter%for subj
            for i=1:size(phaseVector, 1)%for phase
                B.dist(loopCounter,1)  = mean(A.d{s}(phaseVector(i,:) == 1));
                B.phase(loopCounter,1) = i;
                B.cond(loopCounter,1)  = 1;
                B.sn(loopCounter,1)    = s;
                BP.dist(loopCounter,1)  = mean(AP.d{s}(phaseVector(i,:) == 1));
                BP.phase(loopCounter,1) = i;
                BP.cond(loopCounter,1)  = 1;
                BP.sn(loopCounter,1)    = s;
                loopCounter = loopCounter + 1;
            end%for phase
        end%for subj
        
        loopCounter = 1;
        
        %%% Decoding
        for s=1:nIter%for subj
            for i=1:size([A.prepAccuracy; A.prodAccuracy], 2)
                C.acc(loopCounter,1)   = A.prepAccuracy(s, i);
                C.phase(loopCounter,1) = 1;%prep
                C.cond(loopCounter,1)  = i;
                C.sn(loopCounter,1)    = s;
                loopCounter = loopCounter + 1;
                
                C.acc(loopCounter,1)   = A.prodAccuracy(s, i);
                C.phase(loopCounter,1) = 2;%prod
                C.cond(loopCounter,1)  = i;
                C.sn(loopCounter,1)    = s;
                loopCounter = loopCounter + 1;
            end
        end%for subj
        
        figure %%% Distances
        T = tapply(B,{'sn', 'phase'},{'dist', 'mean', 'name', 'dist'});
        colour={[0 0 0], [1 1 1], [0.6 0.6 0.6]};
        barplot(T.phase, T.dist, 'split', T.phase, 'facecolor', colour)
        ylabel('Crossnobis dissimilarity')
        xticklabels({'Prep', 'Prod', 'Cross'})
        title(['prepProd =', num2str(prepProd), '  iter=', num2str(nIter), '  ORD=',num2str(vord), '  TEMP=',num2str(vtemp), '  INTER=',num2str(vinter), '  NOISE=',num2str(vnoise)]);
        
        figure %%% Distances prewhitened
        T = tapply(BP,{'sn', 'phase'},{'dist', 'mean', 'name', 'dist'});
        colour={[0 0 0], [1 1 1], [0.6 0.6 0.6]};
        barplot(T.phase, T.dist, 'split', T.phase, 'facecolor', colour)
        ylabel('Crossnobis dissimilarity pwh')
        xticklabels({'Prep', 'Prod', 'Cross'})
        title(['prepProd =', num2str(prepProd), '  iter=', num2str(nIter), '  ORD=',num2str(vord), '  TEMP=',num2str(vtemp), '  INTER=',num2str(vinter), '  NOISE=',num2str(vnoise)]);
        %-------------------------------------------------------------------------------%
        
        figure %%% Decoding
        T = tapply(C,{'sn', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'});
        colour={[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
        barplot([T.phase, T.cond], T.acc, 'split', T.cond, 'facecolor', colour)
        ylabel('Decoding accuracy')
        xticklabels({'Ord prep', 'Tim prep', 'Int prep', 'Ord prod', 'Tim prod', 'Int prod'})
        title(['prepProd =', num2str(prepProd), '  iter=', num2str(nIter), '  ORD=',num2str(vord), '  TEMP=',num2str(vtemp), '  INTER=',num2str(vinter), '  NOISE=',num2str(vnoise)]);
        %-------------------------------------------------------------------------------%
        
        labels = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        
        %%%Second Moment Matrix
        for s=1:length(A.data)%for iter
            G_hat(:,:,s)= A.G{s}; %extract second moment matrix
            G_hatP(:,:,s)= AP.G{s}; %extract second moment matrix
        end%forIter
        Gm = mean(G_hat,3); % Mean estimate
        GmP = mean(G_hatP,3); % Mean estimate
        
        figure%%%RDM
        %subplot(1,2,1);
        ind=indicatorMatrix('allpairs', 1:length(Gm)); %indicator matrix for all pairwise distances (k*(k-1))
        imagesc(rsa.rdm.squareRDM(diag(ind*Gm*ind'))); %display
        %multiplying variance/covariance (second moment matrix) by
        %indicator matrix results in dissimilarity values (crossnobis)
        title(['prepProd =', num2str(prepProd), '  iter=', num2str(nIter), '  ORD=',num2str(vord), '  TEMP=',num2str(vtemp), '  INTER=',num2str(vinter), '  NOISE=',num2str(vnoise)]);
        xticks(1:8)
        xticklabels(labels);%p=prep, P=prod
        yticks(1:8)
        yticklabels(labels)
        
        figure%%%RDM
        %subplot(1,2,1);
        ind=indicatorMatrix('allpairs', 1:length(GmP)); %indicator matrix for all pairwise distances (k*(k-1))
        imagesc(rsa.rdm.squareRDM(diag(ind*GmP*ind'))); %display
        %multiplying variance/covariance (second moment matrix) by
        %indicator matrix results in dissimilarity values (crossnobis)
        title(['Prewhitened: prepProd =', num2str(prepProd), '  iter=', num2str(nIter), '  ORD=',num2str(vord), '  TEMP=',num2str(vtemp), '  INTER=',num2str(vinter), '  NOISE=',num2str(vnoise)]);
        xticks(1:8)
        xticklabels(labels);%p=prep, P=prod
        yticks(1:8)
        yticklabels(labels)
        %-------------------------------------------------------------------------------%
        
        
        figure%%%MDS
        GRegion = cat(3,A.G{:}); %extract region variance/covariance matrices
        GRegionP = cat(3,AP.G{:}); %extract region variance/covariance matrices
        lc=1;
        for j=1:length(GRegion) %get MDS for each subj in the region
            [SubjCOORD(:,:,lc),~] = pcm_classicalMDS(GRegion(:,:,lc));
            [SubjCOORDP(:,:,lc),~] = pcm_classicalMDS(GRegionP(:,:,lc));
            lc = lc+1;
        end
        COORD = mean(SubjCOORD, 3);
        COORDP = mean(SubjCOORDP, 3);
        xCoord = COORD(1:end,2);%PC2
        yCoord = COORD(1:end,3);%PC3
        zCoord = COORD(1:end,1);%PC1
        %3D scatter plot
        scatterplot3(xCoord,yCoord,zCoord,'split',spl, ... %here we plot PC 2&3 first because
            'markersize',ms, 'markercolor',rsaGB, 'markerfill',rsaGB, 'label',label);%PC 1 is just activity differences
        %%%Draw coloured lines between distinct order & timing conditions
        %Prep
        colors = decodeBRG{1}; %blue for order
        indx=[1 3]';
        line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
        indx=[2 4]';
        line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
        colors = decodeBRG{2}; %red for timing
        indx=[1 2]';
        line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
        indx=[3 4]';
        line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
        %Prod
        colors = decodeBRG{1}; %blue for order
        indx=[5 7]';
        line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
        indx=[6 8]';
        line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
        colors = decodeBRG{2}; %red for timing
        indx=[5 6]';
        line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
        indx=[7 8]';
        line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
        
        grid off
        hold on; plot3(0,0,0,'+','MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',ms+3, 'LineWidth',lw);
        % now draw the vertical errorbar for each point
        eCoord = std(SubjCOORD,[],3);%stdev
        ex = eCoord(1:end,2)/sqrt(length(eCoord(1:end,2)));%standard error to plot
        ey = eCoord(1:end,3)/sqrt(length(eCoord(1:end,3)));
        ez = eCoord(1:end,1)/sqrt(length(eCoord(1:end,1)));
        for i=1:length(COORD(1:end,2))%for length of PC2 vector
            xV = [xCoord(i); xCoord(i)];
            yV = [yCoord(i); yCoord(i)];
            zV = [zCoord(i); zCoord(i)];
            
            xMin = xCoord(i) + ex(i);
            xMax = xCoord(i) - ex(i);
            yMin = yCoord(i) + ey(i);
            yMax = yCoord(i) - ey(i);
            zMin = zCoord(i) + ez(i);
            zMax = zCoord(i) - ez(i);
            
            xB = [xMin, xMax];
            yB = [yMin, yMax];
            zB = [zMin, zMax];
            
            % draw error bars
            h=plot3(xV, yV, zB, '-k');
            set(h, 'LineWidth', 2);
            h=plot3(xB, yV, zV, '-k');
            set(h, 'LineWidth', 2);
            h=plot3(xV, yB, zV, '-k');
            set(h, 'LineWidth', 2);
        end
        
        hold off; xlabel('PC 2'); ylabel('PC 3'); zlabel('PC 1'); set(gca,'fontsize',12);
        %axis equal
        title('multi-dimensional scaling')
        %xlim([-0.035, 0.035])
        %ylim([-0.027, 0.029])
        %zlim([-0.02, 0.3])
        view(30,30)
        
        
        %3D scatter plot - prewhitened
        figure
        xCoordP = COORDP(1:end,2);%PC2
        yCoordP = COORDP(1:end,3);%PC3
        zCoordP = COORDP(1:end,1);%PC1
        %3D scatter plot
        scatterplot3(xCoordP,yCoordP,zCoordP,'split',spl, ... %here we plot PC 2&3 first because
            'markersize',ms, 'markercolor',rsaGB, 'markerfill',rsaGB, 'label',label);%PC 1 is just activity differences
        %%%Draw coloured lines between distinct order & timing conditions
        %Prep
        colors = decodeBRG{1}; %blue for order
        indx=[1 3]';
        line(COORDP(indx,2),COORDP(indx,3),COORDP(indx,1),'color',colors, 'linewidth',lw);
        indx=[2 4]';
        line(COORDP(indx,2),COORDP(indx,3),COORDP(indx,1),'color',colors, 'linewidth',lw);
        colors = decodeBRG{2}; %red for timing
        indx=[1 2]';
        line(COORDP(indx,2),COORDP(indx,3),COORDP(indx,1),'color',colors, 'linewidth',lw);
        indx=[3 4]';
        line(COORDP(indx,2),COORDP(indx,3),COORDP(indx,1),'color',colors, 'linewidth',lw);
        %Prod
        colors = decodeBRG{1}; %blue for order
        indx=[5 7]';
        line(COORDP(indx,2),COORDP(indx,3),COORDP(indx,1),'color',colors, 'linewidth',lw);
        indx=[6 8]';
        line(COORDP(indx,2),COORDP(indx,3),COORDP(indx,1),'color',colors, 'linewidth',lw);
        colors = decodeBRG{2}; %red for timing
        indx=[5 6]';
        line(COORDP(indx,2),COORDP(indx,3),COORDP(indx,1),'color',colors, 'linewidth',lw);
        indx=[7 8]';
        line(COORDP(indx,2),COORDP(indx,3),COORDP(indx,1),'color',colors, 'linewidth',lw);
        
        grid off
        hold on; plot3(0,0,0,'+','MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',ms+3, 'LineWidth',lw);
        % now draw the vertical errorbar for each point
        %         eCoordP = std(SubjCOORDP,[],3);%stdev
        %         ex = eCoordP(1:end,2)/sqrt(length(eCoordP(1:end,2)));%standard error to plot
        %         ey = eCoordP(1:end,3)/sqrt(length(eCoordP(1:end,3)));
        %         ez = eCoordP(1:end,1)/sqrt(length(eCoordP(1:end,1)));
        %         for i=1:length(COORD(1:end,2))%for length of PC2 vector
        %             xV = [xCoordP(i); xCoordP(i)];
        %             yV = [yCoordP(i); yCoordP(i)];
        %             zV = [zCoordP(i); zCoordP(i)];
        %
        %             xMin = xCoordP(i) + ex(i);
        %             xMax = xCoordP(i) - ex(i);
        %             yMin = yCoordP(i) + ey(i);
        %             yMax = yCoordP(i) - ey(i);
        %             zMin = zCoordP(i) + ez(i);
        %             zMax = zCoordP(i) - ez(i);
        %
        %             xB = [xMin, xMax];
        %             yB = [yMin, yMax];
        %             zB = [zMin, zMax];
        %
        %             % draw error bars
        %             h=plot3(xV, yV, zB, '-k');
        %             set(h, 'LineWidth', 2);
        %             h=plot3(xB, yV, zV, '-k');
        %             set(h, 'LineWidth', 2);
        %             h=plot3(xV, yB, zV, '-k');
        %             set(h, 'LineWidth', 2);
        %         end
        
        hold off; xlabel('PC 2'); ylabel('PC 3'); zlabel('PC 1'); set(gca,'fontsize',12);
        %axis equal
        title('prewhitened multi-dimensional scaling')
        %xlim([-0.035, 0.035])
        %ylim([-0.027, 0.029])
        %zlim([-0.02, 0.3])
        view(30,30)
        
    case 'subAndCB_snr_ratio' %calculate signal to noise ratio for each region, to match later simulations
        
        regNames = {'lTha', 'lCau', 'lPut', 'lHip', 'rHip', 'rLob4', 'rLob5'};
        
        if ~exist(fullfile(simuDir, 'snr_regions.mat'),'file')
            %%%Regions to load
            %lTha   %lCau   %lPut   %lHip   %rHip
            targetSubNums  = [ 1       2       3       4       8    ];
            %rLob4   %rLob5
            targetCBNums   = [ 2        4     ];
            
            
            %%%Subcortical ----------------------------------------------------
            %%%Load
            dataFile = fullfile(roiSubDir, 'preWhitened_betas.mat');
            Sub = load(dataFile, 'snr', 'betaW', 'region', 'SN');%load it
            
            mask = ismember(Sub.region, targetSubNums); %remove non-target regions
            names = fieldnames(Sub);
            for i=1:numel(names)
                Sub.(names{i})(~mask) = [];
            end
            Sub.region(Sub.region == 8) = 5; %change rHip to region 5 for plotting
            
            %%%Cerebellum -----------------------------------------------------
            %%%Load
            dataFile = fullfile(roiCbDir, 'preWhitened_betas.mat');
            CB = load(dataFile, 'snr', 'betaW', 'region', 'SN');%load it
            
            mask = ismember(CB.region, targetCBNums); %remove non-target regions
            names = fieldnames(CB);
            for i=1:numel(names)
                CB.(names{i})(~mask) = [];
            end
            CB.region(CB.region == 2) = 6; %change rLob4 to region 6 for plotting
            CB.region(CB.region == 4) = 7; %change rLob5 to region 7 for plotting
            
            R = addstruct(Sub, CB); %concatenate subcortical & CB
            clear('Sub','CB')
            
            R.nVox = cell2mat(cellfun(@length,R.betaW,'uni',false));
            R = rmfield(R, 'betaW');
            %save
            save(fullfile(simuDir, 'snr_regions.mat'), 'R')
        else
            load(fullfile(simuDir, 'snr_regions.mat'),'R')
        end
        
        figure %plot signal to noise ratio for all regions
        subplot(1,2,1)
        myboxplot(R.region, R.snr, 'xtickoff');
        ylabel('Signal to noise ratio')
        xticklabels(regNames)
        
        subplot(1,2,2)
        myboxplot(R.region, R.nVox, 'xtickoff');
        ylabel('Number of voxels')
        xticklabels(regNames)
    case 'simulations_regional' %simulates and analyses activity patterns with snr matched to target subcortical regions - slow!
        
        %%%Runs 1000 iterations by default, so takes some time. Change
        %%%varargin 'nIter' to a lower number if needed.
        
        %%%Load signal to noise ratios for target regions
        %%%Generated by case 'subAndCB_snr_ratio'
        load(fullfile(simuDir, 'snr_regions.mat'), 'R')
        R=tapply(R,{'region', 'SN', 'nVox', 'snr'},{'snr', 'mean', 'name', 'snr'});%just to reorder to region first
        Rmean = tapply(R,{'region'},{'snr', 'mean', 'name', 'snr'},{'nVox','mean','name','nVox'});
        Rmean.nVox = round(Rmean.nVox);
        
        %%%Load percent signal change to use as prodscaling variable
        %                  lTha    lCau    lPut    lHip    rHip
        targetSubNums  = [ 1       2       3       4       8    ];
        %                  rLob4    rLob5
        targetCBNums   = [ 2        4     ];
        
        %subcortical
        load(fullfile(rsaDir, 'subcortical', 'subRoiDistances.mat'), 'Rperc')
        subRperc = tapply(Rperc,{'region','SN','phase'},{'perc','mean','name','perc'});
        mask = ismember(subRperc.region, targetSubNums); %remove non-target regions
        names = fieldnames(subRperc);
        for i=1:numel(names)
            subRperc.(names{i})(~mask) = [];
        end
        subRperc.region(subRperc.region == 8) = 5; %change rHip to region 5 for plotting
        clear Rperc
        
        %Cerebellum
        load(fullfile(rsaDir, 'cerebellum', 'cbRoiDistances.mat'), 'Rperc')
        cbRperc = tapply(Rperc,{'region','SN','phase'},{'perc','mean','name','perc'});
        mask = ismember(cbRperc.region, targetCBNums); %remove non-target regions
        names = fieldnames(cbRperc);
        for i=1:numel(names)
            cbRperc.(names{i})(~mask) = [];
        end
        cbRperc.region(cbRperc.region == 2) = 6; %change rLob4 to region 6 for plotting
        cbRperc.region(cbRperc.region == 4) = 7; %change rLob5 to region 7 for plotting
        
        Rperc = addstruct(subRperc, cbRperc); %concatenate subcortical & CB
        clear('subRperc','cbRperc')
        
        %Calculate scaling factor to increase/decrease production by
        Rperctemp   = tapply(Rperc,{'region','phase'},{'perc','mean','name','perc'});
        prodscaling = Rperctemp.perc(Rperctemp.phase == 2) - Rperctemp.perc(Rperctemp.phase == 1); %prep-prod diff
        
        regNames = {'lTha', 'lCau', 'lPut', 'lHip', 'rHip', 'rLob4', 'rLob5'};%change accordingly
        
        %Other input options that can alternatively be specified by user
        prepProd    = 0;         %0: prep = prod, 1: prep ~= prod
        vord        = [0.0005 0.0005];    %distance between two orders - provide two values if prepProd == 1
        vtemp       = [0.0005 0.0005];    %distance between two timings - provide two values if prepProd == 1
        vinter      = [0.0005 0.0005];    %distance between four sequences - provide two values if prepProd == 1
        nIter       = 1000;       %iterations of simulation
        vararginoptions(varargin,{'prepProd','prodscaling','vord','vtemp','vinter','nIter'});
        
        A=[]; A.data=[]; AP.data=[]; A.snr=[]; A.tgtSnr=[]; A.vnoise=[]; A.region=[]; AP.region=[]; A.SN=[]; AP.SN=[]; %data struct
        A.G=[]; AP.G=[]; A.d=[]; A.Sig=[]; AP.d=[]; AP.Sig=[]; A.prepAccuracy=[]; A.prodAccuracy=[];
        B=[]; %output struct - stores distances, LDA, G matrix
        D=[]; %raw data store
        
        %%%Prepare condition and partition (run) vectors
        nrruns = 6; nCond = 8;
        runC      = repmat(1:nrruns,nCond,1); %labelling runs 1:6 (imaging runs)
        runC      = reshape(runC,1,[])';
        condC     = repmat((1:nCond)', nrruns, 1); %labelling conditions 1:4 (sequences)
        
        for r=1:length(Rmean.region)%for region
            for i=1:nIter%for 'subj'
                [Y, Y_hat, ~, snr, tgtSnr, noise] = prepProdSimu_makedataPP_region(Rmean.nVox(r), Rmean.snr(r), 'prepProd',prepProd,'prodscaling',prodscaling(r),...
                    'vord',vord,'vtemp',vtemp,'vinter',vinter);
                
                A.data{end+1}      = Y;
                AP.data{end+1}     = Y_hat;
                A.snr(end+1,:)     = snr;
                A.tgtSnr(end+1,:)  = tgtSnr;
                A.vnoise(end+1,:)  = noise;
                A.region(end+1,:)  = r;
                AP.region(end+1,:) = r;
                A.SN(end+1,:)      = i;
                AP.SN(end+1,:)     = i;
                
                %%%Second moment matrix
                [A.G{end+1,1},~]    = pcm_estGCrossval(A.data{end},runC,condC);
                [AP.G{end+1,1},~]   = pcm_estGCrossval(AP.data{end},runC,condC);
                
                %%%Distances
                [A.d{end+1},     A.Sig{end+1}]     = rsa.distanceLDC(A.data{end}, runC, condC);
                [AP.d{end+1},    AP.Sig{end+1}]    = rsa.distanceLDC(AP.data{end}, runC, condC);
                
                %%%Linear decoding
                prepData = [A.data{end}(1:4,:); A.data{end}(9:12,:);  A.data{end}(17:20,:); A.data{end}(25:28,:); A.data{end}(33:36,:); A.data{end}(41:44,:)];
                prodData = [A.data{end}(5:8,:); A.data{end}(13:16,:); A.data{end}(21:24,:); A.data{end}(29:32,:); A.data{end}(37:40,:); A.data{end}(45:48,:)];
                
                %               accuracy columns
                %1: order,            2: timing,           3:integrated
                [A.prepAccuracy(end+1,1), A.prepAccuracy(end+1,2), A.prepAccuracy(end+1,3)] = prepProdSimu_classify(prepData);
                [A.prodAccuracy(end+1,1), A.prodAccuracy(end+1,2), A.prodAccuracy(end+1,3)] = prepProdSimu_classify(prodData);
            end%for 'subj'
        end%for region
        
        figure %plot signal to noise ratio for all regions
        subplot(1,2,1)
        myboxplot(R.region, R.snr, 'xtickoff');
        ylabel('Empirical signal to noise ratio')
        ylim([0 0.7])
        xticklabels(regNames)
        
        subplot(1,2,2)
        myboxplot(A.region, A.snr, 'xtickoff');
        ylabel('Simulated signal to noise ratio')
        ylim([0 0.7])
        xticklabels(regNames)
        %-------------------------------------------------------------------------------%
        
        
        %%%Distances
        phaseVector = [...%index (1s are where phase is different)
            1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...%prep
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1; ...%prod
            0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0  ...%cross
            ];
        %Loop over subjs, conditions, and phases, and store distances as
        %struct to plot later
        loopCounter = 1;
        for r=unique(A.region)'%for region
            Atemp = A.d(A.region == r);
            APtemp = AP.d(A.region == r);
            for s=1:nIter%for subj
                for i=1:size(phaseVector, 1)%for phase
                    B.dist(loopCounter,1)   = mean(Atemp{s}(phaseVector(i,:) == 1));
                    B.phase(loopCounter,1)  = i;
                    B.cond(loopCounter,1)   = 1;
                    B.sn(loopCounter,1)     = s;
                    B.region(loopCounter,1) = r;
                    BP.dist(loopCounter,1)   = mean(APtemp{s}(phaseVector(i,:) == 1));
                    BP.phase(loopCounter,1)  = i;
                    BP.cond(loopCounter,1)   = 1;
                    BP.sn(loopCounter,1)     = s;
                    BP.region(loopCounter,1) = r;
                    loopCounter = loopCounter + 1;
                end%for phase
            end%for subj
        end%for region
        
        
        %%% Decoding
        loopCounter = 1;
        for r=unique(A.region)'%for region
            APrepTemp = A.prepAccuracy(A.region == r,:);
            AProdTemp = A.prodAccuracy(A.region == r,:);
            for s=1:nIter%for subj
                for i=1:size([APrepTemp; AProdTemp], 2)%for cond
                    C.acc(loopCounter,1)    = APrepTemp(s, i);
                    C.phase(loopCounter,1)  = 1;%prep
                    C.cond(loopCounter,1)   = i;
                    C.sn(loopCounter,1)     = s;
                    C.region(loopCounter,1) = r;
                    loopCounter = loopCounter + 1;
                    
                    C.acc(loopCounter,1)    = AProdTemp(s, i);
                    C.phase(loopCounter,1)  = 2;%prod
                    C.cond(loopCounter,1)   = i;
                    C.sn(loopCounter,1)     = s;
                    C.region(loopCounter,1) = r;
                    loopCounter = loopCounter + 1;
                end%for cond
            end%for subj
        end%for region
        
        %%%Plot - Prepare variables
        %Colour      %Blue(ord)         %Red(tim)               %Green(int)
        decodeBRG = {[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
        rsaGB     = {[0.8 0.8 0.8], [0 0 0]};
        
        %Style
        spl = [zeros(4,1); ones(4,1)];
        label = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        ms = 12;
        ls = 20;
        lw = 3;
        
        %y axis scales
        pscY     = [-0.12,  0.7];
        rsaY     = [-0.003, 0.0045];
        ldaY     = [-1.52,  2];
        crossY   = [0,      0.11];
        mdsDistY = [0       0.07];
        
        
        %%%Plot - do plotting
        figure %%% Distances
        for r=unique(B.region')
            subplot(4,2,r)
            T = tapply(B,{'sn', 'phase'},{'dist', 'mean', 'name', 'dist'},'subset', B.region == r);
            colour={[0 0 0], [1 1 1], [0.6 0.6 0.6]};
            barplot(T.phase, T.dist, 'split', T.phase, 'facecolor', colour);
            ylabel('Crossnobis dissimilarity')
            xticklabels({'Prep', 'Prod', 'Cross'})
            if r==1
                xlabel(['pP =', num2str(prepProd), ' iter=', num2str(nIter), ' ORD=',num2str(vord), ' TIM=',num2str(vtemp), ' INT=',num2str(vinter), ' NOISE=',num2str(A.vnoise(1))]);
            end
            title(regNames{r});
        end
        
        figure %%% Distances prewhitened
        for r=unique(BP.region')
            subplot(4,2,r)
            T = tapply(BP,{'sn', 'phase'},{'dist', 'mean', 'name', 'dist'},'subset', B.region == r);
            colour={[0 0 0], [1 1 1], [0.6 0.6 0.6]};
            barplot(T.phase, T.dist, 'split', T.phase, 'facecolor', colour);
            ylabel('Crossnobis dissimilarity prewhitened')
            xticklabels({'Prep', 'Prod', 'Cross'})
            if r==1
                xlabel(['pP =', num2str(prepProd), ' iter=', num2str(nIter), ' ORD=',num2str(vord), ' TIM=',num2str(vtemp), ' INT=',num2str(vinter), ' NOISE=',num2str(A.vnoise(1))]);
            end
            title(regNames{r});
        end
        %-------------------------------------------------------------------------------%
        
        figure %%% Decoding
        for r=unique(C.region')
            subplot(4,2,r)
            T = tapply(C,{'sn', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'},'subset', B.region == r);
            colour={[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
            barplot([T.phase, T.cond], T.acc, 'split', T.cond, 'facecolor', colour);
            ylabel('Decoding accuracy')
            xticklabels({'Ord prep', 'Tim prep', 'Int prep', 'Ord prod', 'Tim prod', 'Int prod'})
            if r==1
                xlabel(['pP =', num2str(prepProd), ' iter=', num2str(nIter), ' ORD=',num2str(vord), ' TIM=',num2str(vtemp), ' INT=',num2str(vinter), ' NOISE=',num2str(A.vnoise(1))]);
            end
            title(regNames{r});
        end
        %-------------------------------------------------------------------------------%
        
        
        %%% RDMs and multi-dimensional scaling
        labels = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        %Extract data
        for s=1:length(A.G)%for subj * region
            G(:,:,s)  = A.G{s}; %extract representational dissimilarity matrix into 3D matrix
            GP(:,:,s) = AP.G{s}; %extract representational dissimilarity matrix into 3D matrix
        end%for subj * region
        
        
        figure%%%RDM
        for i=unique(A.region')%plots left hem on the left, right hem on the right
            
            GRegion = G(:,:,A.region == i); %extract region variance/covariance matrices
            GmRegion = mean(GRegion, 3); %mean across subjs
            
            subplot(4, 2, i) %%%RDM for each subcortical region
            ind=indicatorMatrix('allpairs',1:8); %matrix for all pairwise distances (k*(k-1))
            imagesc(rsa.rdm.squareRDM(diag(ind*GmRegion*ind'))); %display
            %multiplying variance/covariance by indicator matrix results in
            %dissimilarity values (crossnobis)
            title([regNames{i} ' RDM (crossnobis)'])
        end
        figure %%%Prewhitened RDM
        for i=unique(AP.region')%for region
            
            GRegion = GP(:,:,AP.region == i); %extract region variance/covariance matrices
            GmRegion = mean(GRegion, 3); %mean across subjs
            
            subplot(4, 2, i) %%%RDM for each subcortical region
            ind=indicatorMatrix('allpairs',1:8); %matrix for all pairwise distances (k*(k-1))
            imagesc(rsa.rdm.squareRDM(diag(ind*GmRegion*ind'))); %display
            %multiplying variance/covariance by indicator matrix results in
            %dissimilarity values (crossnobis)
            colorbar;
            axis image
            title([regNames{i} ' Prewhitened RDM (crossnobis)'])
        end
        %----------------------------------------------------------------------------------------------%
        
        
        %%%Multi-dimensional scaling
        figure %plot MDS
        for r=unique(A.region')%for region
            %figure %plot MDS
            subplot(4,2,r)
            GRegion = G(:,:,A.region == r); %extract region variance/covariance matrices
            for j=1:length(GRegion) %get MDS for each subj in the region
                [SubjCOORD(:,:,j),~] = pcm_classicalMDS(GRegion(:,:,j));
            end
            COORD = mean(SubjCOORD, 3);
            %COORDP = mean(SubjCOORDP, 3);
            xCoord = COORD(1:end,2);%PC2
            yCoord = COORD(1:end,3);%PC3
            zCoord = COORD(1:end,1);%PC1
            %3D scatter plot
            scatterplot3(COORD(1:end,2),COORD(1:end,3),COORD(1:end,1),'split',spl, ... %here we plot PC 2&3 first because
                'markersize',ms, 'markercolor',rsaGB, 'markerfill',rsaGB, 'label',label);%PC 1 is just activity differences
            
            %%%Draw coloured lines between distinct order & timing conditions
            %Prep
            colors = decodeBRG{1}; %blue for order
            indx=[1 3]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[2 4]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            colors = decodeBRG{2}; %red for timing
            indx=[1 2]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[3 4]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            
            %Prod
            colors = decodeBRG{1}; %blue for order
            indx=[5 7]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[6 8]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            colors = decodeBRG{2}; %red for timing
            indx=[5 6]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[7 8]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            
            grid off
            hold on; plot3(0,0,0,'+','MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',ms+3, 'LineWidth',lw);
            %             % now draw the vertical errorbar for each point
            %             eCoord = std(SubjCOORD,[],3);%stdev
            %             ex = eCoord(1:end,2)/sqrt(length(eCoord(1:end,2)));%standard error to plot
            %             ey = eCoord(1:end,3)/sqrt(length(eCoord(1:end,3)));
            %             ez = eCoord(1:end,1)/sqrt(length(eCoord(1:end,1)));
            %             for i=1:length(COORD(1:end,2))%for length of PC2 vector
            %                 xV = [xCoord(i); xCoord(i)];
            %                 yV = [yCoord(i); yCoord(i)];
            %                 zV = [zCoord(i); zCoord(i)];
            %
            %                 xMin = xCoord(i) + ex(i);
            %                 xMax = xCoord(i) - ex(i);
            %                 yMin = yCoord(i) + ey(i);
            %                 yMax = yCoord(i) - ey(i);
            %                 zMin = zCoord(i) + ez(i);
            %                 zMax = zCoord(i) - ez(i);
            %
            %                 xB = [xMin, xMax];
            %                 yB = [yMin, yMax];
            %                 zB = [zMin, zMax];
            %
            %                 % draw error bars
            %                 h=plot3(xV, yV, zB, '-k');
            %                 set(h, 'LineWidth', 2);
            %                 h=plot3(xB, yV, zV, '-k');
            %                 set(h, 'LineWidth', 2);
            %                 h=plot3(xV, yB, zV, '-k');
            %                 set(h, 'LineWidth', 2);
            %             end%for length of PC2 vector
            hold off; xlabel('PC 2'); ylabel('PC 3'); zlabel('PC 1'); set(gca,'fontsize',12);
            %axis equal
            title([regNames{r} ' MDS'])
            %ylim([-0.03, 0.04])
            %xlim([-0.02, 0.16])
            view(30,30)
        end%for region
        
        %%%Pre-whitened
        %figure %plot MDS
        for r=unique(AP.region')%for region
            %subplot(4,2,r)
            figure %plot MDS
            GRegion = GP(:,:,AP.region == r); %extract region variance/covariance matrices
            for j=1:length(GRegion) %get MDS for each subj in the region
                [SubjCOORD(:,:,j),~] = pcm_classicalMDS(GRegion(:,:,j));
            end
            COORD = mean(SubjCOORD, 3);
            %COORDP = mean(SubjCOORDP, 3);
            xCoord = COORD(1:end,2);%PC2
            yCoord = COORD(1:end,3);%PC3
            zCoord = COORD(1:end,1);%PC1
            %3D scatter plot
            scatterplot3(COORD(1:end,2),COORD(1:end,3),COORD(1:end,1),'split',spl, ... %here we plot PC 2&3 first because
                'markersize',ms, 'markercolor',rsaGB, 'markerfill',rsaGB, 'label',label);%PC 1 is just activity differences
            
            %%%Draw coloured lines between distinct order & timing conditions
            %Prep
            colors = decodeBRG{1}; %blue for order
            indx=[1 3]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[2 4]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            colors = decodeBRG{2}; %red for timing
            indx=[1 2]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[3 4]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            
            %Prod
            colors = decodeBRG{1}; %blue for order
            indx=[5 7]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[6 8]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            colors = decodeBRG{2}; %red for timing
            indx=[5 6]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            indx=[7 8]';
            line(COORD(indx,2),COORD(indx,3),COORD(indx,1),'color',colors, 'linewidth',lw);
            
            grid off
            hold on; plot3(0,0,0,'+','MarkerFaceColor', [0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerSize',ms+3, 'LineWidth',lw);
            %             % now draw the vertical errorbar for each point
            %             eCoord = std(SubjCOORD,[],3);%stdev
            %             ex = eCoord(1:end,2)/sqrt(length(eCoord(1:end,2)));%standard error to plot
            %             ey = eCoord(1:end,3)/sqrt(length(eCoord(1:end,3)));
            %             ez = eCoord(1:end,1)/sqrt(length(eCoord(1:end,1)));
            %             for i=1:length(COORD(1:end,2))%for length of PC2 vector
            %                 xV = [xCoord(i); xCoord(i)];
            %                 yV = [yCoord(i); yCoord(i)];
            %                 zV = [zCoord(i); zCoord(i)];
            %
            %                 xMin = xCoord(i) + ex(i);
            %                 xMax = xCoord(i) - ex(i);
            %                 yMin = yCoord(i) + ey(i);
            %                 yMax = yCoord(i) - ey(i);
            %                 zMin = zCoord(i) + ez(i);
            %                 zMax = zCoord(i) - ez(i);
            %
            %                 xB = [xMin, xMax];
            %                 yB = [yMin, yMax];
            %                 zB = [zMin, zMax];
            %
            %                 % draw error bars
            %                 h=plot3(xV, yV, zB, '-k');
            %                 set(h, 'LineWidth', 2);
            %                 h=plot3(xB, yV, zV, '-k');
            %                 set(h, 'LineWidth', 2);
            %                 h=plot3(xV, yB, zV, '-k');
            %                 set(h, 'LineWidth', 2);
            %             end%for length of PC2 vector
            hold off; xlabel('PC 2'); ylabel('PC 3'); zlabel('PC 1'); set(gca,'fontsize',12);
            %axis equal
            title([regNames{r} ' MDS prewhitened pP ' num2str(prepProd)])
            %ylim([-0.03, 0.04])
            %xlim([-0.02, 0.16])
            view(30,30)
        end%for region
        %-------------------------------------------------------------------------------%
        
        %%%Calculate euclidean distances
        eucPcDist=[]; A.simEucDist = nan(length(A.SN),1); A.simEucScaling = nan(length(A.SN),1);
        for i=unique(A.region)'
            COORD = [];
            GRegion = G(:,:,A.region == i); %extract region variance/covariance matrices
            
            for j=1:length(GRegion) %get MDS for each subj in the region
                [COORD(:,:,j),~] = pcm_classicalMDS(GRegion(:,:,j));
            end
            
            PCs = COORD(:,2:3,:); %PCs 2 and 3
            
            %Euclidean distance between phases within sequences for PC2
            %and PC3
            for j=1:length(PCs)%for subjs
                for k=1:4%for within-sequence prep vs prod
                    eucPcDist(k) = pdist([PCs(k,:,j)' PCs(k+4,:,j)'], 'euclidean');
                end
                
                %Calculate the size of the quadrangle for prep and prod to
                %normalise the distance between phases despite ord temp
                %and int increases or decreases
                
                %pairwise distance measures k(k-1) / 2
                prepDists(1) = pdist([PCs(1,:,j)', PCs(2,:,j)'],'euclidean');
                prepDists(2) = pdist([PCs(1,:,j)', PCs(3,:,j)'],'euclidean');
                prepDists(3) = pdist([PCs(1,:,j)', PCs(4,:,j)'],'euclidean');
                prepDists(4) = pdist([PCs(2,:,j)', PCs(3,:,j)'],'euclidean');
                prepDists(5) = pdist([PCs(2,:,j)', PCs(4,:,j)'],'euclidean');
                prepDists(6) = pdist([PCs(3,:,j)', PCs(4,:,j)'],'euclidean');
                
                prodDists(1) = pdist([PCs(5,:,j)', PCs(6,:,j)'],'euclidean');
                prodDists(2) = pdist([PCs(5,:,j)', PCs(7,:,j)'],'euclidean');
                prodDists(3) = pdist([PCs(5,:,j)', PCs(8,:,j)'],'euclidean');
                prodDists(4) = pdist([PCs(6,:,j)', PCs(7,:,j)'],'euclidean');
                prodDists(5) = pdist([PCs(6,:,j)', PCs(8,:,j)'],'euclidean');
                prodDists(6) = pdist([PCs(7,:,j)', PCs(8,:,j)'],'euclidean');
                %%%FIND A BETTER WAY TO CODE THIS ^
                
                prepSize = mean(prepDists);
                prodSize = mean(prodDists);
                avgSize  = mean([prepSize, prodSize]);
                
                %Collect variables to later add to RDist
                A.simEucDist(A.region == i & A.SN == j,:)    = mean(eucPcDist);
                A.simEucScaling(A.region == i & A.SN == j,:) = avgSize;
            end
        end
        A.simEucDistScaled = (A.simEucDist ./ A.simEucScaling) * 100;
        
        R = A; clear A; %convert to R struct
        R.data = AP.data; %replace regular data with prewhitened data
        R.G    = AP.G;
        R.d    = AP.d;
        R.Sig  = AP.Sig;
        
        %add input specifications for reference
        R.prepProd = prepProd;
        R.vord     = vord;
        R.vtemp    = vtemp;
        R.vinter   = vinter;
        R.nIter    = nIter;
        R.prodscaling = prodscaling;
        
        maintSwitch = {'maint', 'switch'}; %for filename, whether data switches or maints
        sord = num2str(vord); stemp = num2str(vtemp); sinter = num2str(vinter); snIter = num2str(nIter);
        
        %save(fullfile(simuDir, ['simRoiEucDistances_' ...
        %    maintSwitch{prepProd+1} '_O' sord '_T' stemp '_I' sinter '_N' snIter '.mat']), 'R')
    case 'simu_empirical_regional_comparison'
        %region names
        subcortNames = {'lTha', 'lCau', 'lPut', 'lHip', 'rTha', 'rCau', 'rPut', 'rHip'};
        CBnames      = {'lLob4', 'rLob4', 'lLob5', 'rLob5', 'lLob6', 'rLob6', 'lCru1', 'rCru1', 'lCru2', 'rCru2'};
        targetSubNames = {'lTha', 'lCau', 'lPut', 'lHip', 'rHip'};
        targetSubNums  = [ 1       2       3       4       8    ];
        targetCBNames  = {'rLob4', 'rLob5'};
        targetCBNums   = [ 2        4     ];
        
        targetRegNames = [targetSubNames, targetCBNames];
        
        %%%Load simulated data
        simuDataFile = fullfile(simuDir, 'simRoiEucDistances_maint_O0.0005      0.0005_T0.0005      0.0005_I0.0005      0.0005_N1000.mat');
        load(simuDataFile, 'R');
        RsimuMaint = tapply(R,{'region', 'SN', 'simEucScaling','simEucDist'},{'simEucDistScaled', 'mean', 'name', 'distScaled'});
        clear R
        
        simuDataFile = fullfile(simuDir, 'simRoiEucDistances_switch_O0.0005      0.0005_T0.0005      0.0005_I0.0005      0.0005_N1000.mat');
        load(simuDataFile, 'R');
        RsimuSwitch = tapply(R,{'region', 'SN', 'simEucScaling','simEucDist'},{'simEucDistScaled', 'mean', 'name', 'distScaled'});
        clear R
        
        
        %%%Subcortical ----------------------------------------------------
        %%%Load
        dataFile = fullfile(rsaDir, 'subcortical', 'subRoiDistances_withmds.mat');
        load(dataFile, 'Rdist');%load it
        
        mask = ismember(Rdist.region, targetSubNums); %remove non-target regions
        names = fieldnames(Rdist);
        for i=1:numel(names)
            Rdist.(names{i})(~mask) = [];
        end
        Rdist.region(Rdist.region == 8) = 5; %change rHip to region 5 for plotting
        subDist = tapply(Rdist,{'region', 'SN', 'eucScaling'},{'dist', 'mean', 'name', 'dist'},'subset',Rdist.cond == 2);
        clear Rdist
        
        %%%Cerebellum -----------------------------------------------------
        %%%Load
        dataFile = fullfile(rsaDir, 'cerebellum', 'cbRoiDistances_withmds.mat');
        load(dataFile, 'Rdist');%load it
        
        mask = ismember(Rdist.region, targetCBNums); %remove non-target regions
        names = fieldnames(Rdist);
        for i=1:numel(names)
            Rdist.(names{i})(~mask) = [];
        end
        Rdist.region(Rdist.region == 2) = 6; %change rLob4 to region 6 for plotting
        Rdist.region(Rdist.region == 4) = 7; %change rLob5 to region 7 for plotting
        CBdist = tapply(Rdist,{'region', 'SN', 'eucScaling'},{'dist', 'mean', 'name', 'dist'},'subset',Rdist.cond == 2);
        %restricts to euclidean PC distance, excludes crossnobis ^
        
        %concatenate subcortical and CB distance structs
        R = addstruct(subDist, CBdist);
        
        %calculate the scaled euclidean distance based on scaling variable
        R.distScaled = (R.dist ./ R.eucScaling) * 100;
        R            = rmfield(R, {'dist', 'eucScaling'});
        R.distDiff   = nan(length(R.SN), 1);
        
        %%%Use simulations as baseline for one-sample tests
        for i=unique(RsimuMaint.region)'%for region
            baseDist(:,i)   = mean(RsimuMaint.distScaled(RsimuMaint.region == i));%calc maint base
            switchDist(:,i) = mean(RsimuSwitch.distScaled(RsimuSwitch.region == i));%calc switch hypothetical
            R.distDiff(R.region == i) = R.distScaled(R.region == i) - baseDist(i); %calc diff between scaled dist and baseline, for plot
        end%for region
        
        
        figure%Plot empirical region distances, with lines showing simulated levels
        myboxplot(R.region, R.distScaled, 'xtickoff');
        ylabel('Cross-phase Euclidean distance')
        ylim([80 130])
        xticklabels(targetRegNames)
        
        %draw horizontal lines at simulated maint and switch distances for each region
        lineAxes = [0.7 1.3]; maintColour = ([64,224,208] / 255); switchColour = ([242,140,40] / 255);
        for i=unique(RsimuMaint.region)'%for region
            drawline(baseDist(i), 'dir', 'horz', ... %maint lines
                'lim', lineAxes + (i-1), 'color', maintColour, 'linestyle', '-', 'linewidth', 2)
            
            drawline(switchDist(i), 'dir', 'horz', ... %switch lines
                'lim', lineAxes + (i-1), 'color', switchColour, 'linestyle', '-', 'linewidth', 2)
        end%for region
        
        %test significance of distance elevations
        t = NaN(length(unique(R.region)') - 1,1); p = NaN(length(unique(R.region)') - 1,1);
        for i = unique(R.region)'%for region
            [h(i),p(i),ci{i},stats{i}] = ttest_RY(R.distScaled(R.region == i),baseDist(i));
        end%for region
        
        astXLoc = 1:7; %x coordinates for the significance *s on plot
        astYLoc = repelem(125, 7);
        %display asterisks between sig results
        text(astXLoc(p < .05),astYLoc(p < .05),'*','fontSize',20);
        %disp(p)
        %-----------------------------------------------------------------%
        
        
        figure%Plot empirical region distances, with 0 being simulated level
        myboxplot(R.region, R.distDiff, 'xtickoff');
        ylabel('Cross-phase Euclidean distance')
        ylim([-10 30])
        xticklabels(targetRegNames)
        drawline(0, 'dir', 'horz', 'linestyle', '- -')
        
        switchDistToBase = switchDist - baseDist;
        for i=unique(RsimuMaint.region)'%for region
            drawline(switchDistToBase(i), 'dir', 'horz', ... %switch lines
                'lim', lineAxes + (i-1), 'color', switchColour, 'linestyle', '-', 'linewidth', 2)
        end%for region
        
        %test significance of distance elevations
        p = NaN(length(unique(R.region)') - 1,1);
        for i = unique(R.region)'%for region
            [h(i),p(i),ci{i},stats{i}] = ttest_RY(R.distDiff(R.region == i), 0, ...
                'tail', 'right');
            stats{i}.mean = mean(R.distDiff(R.region == i)); %mean of data
            stats{i}.d    = stats{i}.mean / stats{i}.sd; %cohen's d for one sample
            
            p(i) = p(i) * 7; %bonferroni corrections for region
            output = ['t(' num2str(stats{i}.df) ') = ' num2str(stats{i}.tstat) ', p = ' num2str(p(i)) ', d = ' num2str(stats{i}.d)];
            disp(output)%display in commmand window
        end%for region
        
        
        
        astXLoc = 1:7; %x coordinates for the significance *s on plot
        astYLoc = repelem(28, 7);
        %display asterisks between sig results
        text(astXLoc(p < .05),astYLoc(p < .05),'*','fontSize',20);
        %disp(p)
        
        
        
        
end%switch
end%function