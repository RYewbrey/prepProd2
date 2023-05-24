function prepProd2_imana_RY(what,varargin)


%PC paths Load relevant Matlab directories
% addpath(genpath('Z:/rhys/prepProd2/matlab')); %Adjust! loaded with subdirectories (genpath command)
% addpath('Z:/toolboxes/spm12');
% addpath(genpath('Z:/toolboxes/tools')); %joern's extensions for spm
% addpath(genpath('Z:/toolboxes/userfun')); %joern's util tools (open source)
% addpath(genpath('Z:/toolboxes/region-master')); %joern's region toolbox for spm
% addpath(genpath('Z:/toolboxes/spm12/toolbox/suit')); %SUIT Cerebellum
% addpath(genpath('Z:/toolboxes/spm12/toolbox/DARTEL')); %DARTEL deformation for suit reslice
% addpath(genpath('Z:/toolboxes/spm12/toolbox/OldSeg')); %for suit reslice
% addpath(genpath('Z:/toolboxes/permutest')); %permutest for crossSection analysis
% addpath(genpath('Z:/toolboxes/rsatoolbox_matlab')); %RSA toolbox
% addpath(genpath('Z:/toolboxes/pcm_toolbox')); %PCM toolbox
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
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/rsatoolbox_matlab')); %RSA toolbox
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/pcm_toolbox')); %PCM toolbox
% addpath(genpath('/rds/projects/k/kornyshk-kornyshevalab/toolboxes/surfing')); %surfing, used in RSA toolbox

%%%definition and variables

%% Data paths:
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
cerebellumDir=[baseDir filesep 'imaging' filesep 'cerebellum'];
crossSectionDir=[baseDir filesep 'imaging' filesep 'crossSection'];
modelDir=[baseDir filesep 'imaging' filesep 'simulations' filesep 'models'];
rsaDir=[baseDir filesep 'imaging' filesep 'rsa'];
rsaCorticalDir=[rsaDir filesep 'cortical'];
rsaCorticalGroupDir=[baseDir filesep 'imaging' filesep 'rsa_secondlevel' filesep 'cortical'];

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
    1, 1, 1, 1; ... %s01 placeholder
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

allSubcortStructs = ... %subcortical structures of interest for later analysis, defined from freesurfer aseg
    {'left_thalamus','left_caudate','left_putamen','left_pallidum','brain_stem','left_hippocampus',...
    'left_amygdala','left_accumbens_area','left_ventral_dc','left_vessel','left_choroid_plexus',...
    'right_thalamus','right_caudate','right_putamen','right_pallidum','right_hippocampus',...
    'right_amygdala','right_accumbens_area','right_ventral_dc','right_vessel','right_choroid_plexus'
    };

subcortStructs = {... %subcortical structures of interest for later analysis, defined from freesurfer aseg
    10,              11,             12,             13,              17,...
    29,              50,             51,             52,              53;...
    'left_thalamus', 'left_caudate', 'left_putamen', 'left_pallidum', 'left_hippocampus',...
    'right_thalamus','right_caudate','right_putamen','right_pallidum','right_hippocampus',...
    };

suitCBRegions = {...
    3,               4,                5,               7,                8,             10,             11,            13;
    'left_lobule_5', 'right_lobule_5', 'left_lobule_6', 'right_lobule_6', 'left_crus_1', 'right_crus_1', 'left_crus_2', 'right_crus_2'};

corticalRegions = ...
    {'LM1', 'LPMd', 'LSMA', 'LSPC', 'RM1', 'RPMd', 'RSMA', 'RSPC'};

%%%
switch(what)
    
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
        
        %% *1. Coregister meanEPI to anatomical (no reslicing into anat space, dimension of EPI preserved!): %%%
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
                
            end
            
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
        
        
    case 'suit_segment' %%%%%%%%%%%%%%%%%%% BEGINNING OF SUIT DECODING ANALYSIS (CEREBELLUM) %%%%%%%%%%%%%%%%%%%%%%%
        %Isolate the cerebellum of each participant - produces 'c_<source>_pcereb' (cerebellar mask) and '<source>_seg1/2' (grey and white matter respectively) images
        
        sn=varargin{1};
        cd([baseDir '/imaging/anatomicals/' subj_name{sn}]);
        disp(['suit_segmenting ' subj_name{sn}])
        
        anatomical = {[subj_name{sn} '_anatomical.nii']};
        
        suit_isolate_seg(anatomical, 'maskp', 0.2) %change maskp for probability value. Higher = tighter mask. Hand-correct using MRIcron if necessary
    case 'suit_make_mask' %restrict area of analysis to grey matter - produces 'maskbrainSUIT.nii'
        
        s=varargin{1};
        
        if isfolder([baseDir '/imaging/suit/' subj_name{s}]) == 0
            mkdir([baseDir '/imaging/suit/' subj_name{s}])
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
        cd([baseDir '/imaging/anatomicals/' subj_name{sn}]);
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
        cd([baseDir '/imaging/suit/']);
        
        if isfolder(suitGroupDir) == 0
            mkdir(suitGroupDir)
        end
        disp(['suit_reslicing ' subj_name{sn}])
        
        inDir = [suitDir '/' subj_name{sn} '/']; %path to where data is stored (to be normalised)
        outDir = suitGroupDir;
        filenames = {'szacc_Comb_160_Prep', 'szacc_Comb_160_Mov', 'szacc_Spat_160_Prep', 'szacc_Spat_160_Mov', 'szacc_Temp_160_Prep', 'szacc_Temp_160_Mov', 'szacc_Int_160_Prep', 'szacc_Int_160_Mov'};
        
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
        cd([baseDir '/imaging/suit/' subj_name{sn}]);
        
        disp(['suit_reslicing_contrast ' subj_name{sn}])
        
        inDir = groupDir; %path to where data is stored (to be normalised)
        
        contrasts = {'scon_0001', 'scon_0002', 'scon_0003', 'scon_0004', 'scon_0005', 'scon_0006'};
        
        outDir = suitGroupDir;
        
        % prepare files for input
        affine = {[anatDir '/' subj_name{sn} '/' 'Affine_' subj_name{sn} '_anatomical_seg1.mat']};
        flowfield = {[anatDir '/' subj_name{sn} '/' 'u_a_' subj_name{sn} '_anatomical_seg1.nii']};
        
        dataFiles = cell(length(contrasts),1);
        for i=1:length(contrasts)
            dataFiles{i} = [inDir, '/', contrasts{i}, '_', subj_name{sn} '.nii'];
        end
        
        mask = {[anatDir '/' subj_name{sn} '/' 'c_' subj_name{sn} '_anatomical_pcereb.nii']};
        
        outFiles = cell(length(contrasts),1);
        for i=1:length(contrasts)
            outFiles{i} = [outDir, '/', contrasts{i}, '_', subj_name{sn}, '.nii'];
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
        cd([baseDir '/imaging/suit/']);
        
        %         if isdir([baseDir '/imaging/suit/' subj_name{sn}]) == 0
        %             mkdir([baseDir '/imaging/suit/' subj_name{sn}])
        %         end
        disp(['suit_reslicing ' subj_name{sn}])
        
        inDir = [glmDir '/' subj_name{sn} '/']; %path to where data is stored (to be normalised)
        outDir = [baseDir '/imaging/suit/'];
        filenames = {'_szacc_Comb_160_Prep.nii', '_szacc_Comb_160_Mov.nii', '_szacc_Spat_160_Prep.nii', '_szacc_Spat_160_Mov.nii', '_szacc_Temp_160_Prep.nii', '_szacc_Temp_160_Mov.nii', '_szacc_Int_160_Prep.nii', '_szacc_Int_160_Mov.nii'};
        filenamesT = {'spmT_0001.nii', 'spmT_0002.nii', 'spmT_0003.nii', 'spmT_0004.nii', 'spmT_0005.nii', 'spmT_0006.nii'};
        
        %% prepare files for input
        affine = {[anatDir '/' subj_name{sn} '/' 'Affine_' subj_name{sn} '_anatomical_seg1.mat']};
        flowfield = {[anatDir '/' subj_name{sn} '/' 'u_a_' subj_name{sn} '_anatomical_seg1.nii']};
        
        dataFiles = cell(length(filenames),1); %separate loops to fill filenames for classifiers and T contrast maps
        for i=1:length(filenames)
            dataFiles{i} = [inDir, subj_name{sn}, filenames{i}];
        end
        dataFilesT = cell(length(filenamesT),1);
        for i=1:length(filenamesT)
            dataFilesT{i} = [inDir, filenamesT{i}];
        end
        dataFiles = [dataFiles; dataFilesT];
        
        mask = {[anatDir '/' subj_name{sn} '/' 'c_' subj_name{sn} '_anatomical_pcereb.nii']};
        
        outFiles = cell(length(filenames),1);
        for i=1:length(filenames)
            outFiles{i} = [outDir, '/', subj_name{sn}, filenames{i}];
        end
        outFilesT = cell(length(filenamesT),1);
        for i=1:length(filenamesT)
            outFilesT{i} = [outDir, '/', subj_name{sn}, '_' filenamesT{i}];
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
        dataFolder = {'data/Mov','data/Prep','data/PrepProd','data/ProdPrep','data/Error','data/Rest',... %folders containing data to be projected
            'MVA_comb_mov','MVA_comb_prep','MVA_spat_mov','MVA_spat_prep','MVA_temp_mov','MVA_temp_prep','MVA_int_mov','MVA_int_prep'};
        
        for i =1:length(dataFolder)
            cd([suitDir, '/', dataFolder{i}]) %go into folder location for data and to save into.
            
            flatmapVector = suit_map2surf('spmT_0001.nii'); %produces vector information regarding flatmap
            flatmapGifti = gifti(flatmapVector); %convert to gifti format (compatibility with connectome workbench)
            save(flatmapGifti,'spmT_0001.func.gii')
            save('spmT_0001_flatmap.mat','flatmapVector')
        end
        
        %% for within-matlab visualisation of flatmap, load respective 'spmT_0001_flatmap.mat'
        % then use suit_plotflatmap(flatmapVector, 'cmap', hot, 'cscale', [0 7.00], 'threshold', 3.48) replacing numbers with desired values.
    case 'suit_roi_whole'
        
        dataFolder = {'data/Mov','data/Prep','data/PrepProd','data/ProdPrep','data/Error','data/Rest',... %folders containing data to be projected
            'MVA_comb_mov','MVA_comb_prep','MVA_spat_mov','MVA_spat_prep','MVA_temp_mov','MVA_temp_prep','MVA_int_mov','MVA_int_prep'};
        
        dataNames = {'Mov','Prep','PrepProd','ProdPrep','Error','Rest',... %names for saving
            'MVA_comb_mov','MVA_comb_prep','MVA_spat_mov','MVA_spat_prep','MVA_temp_mov','MVA_temp_prep','MVA_int_mov','MVA_int_prep'};
        
        for i =1:length(dataFolder)
            
            cd([suitDir, '/', dataFolder{i}]) %go into folder location for data.
            saveDir = [suitDir, '/', 'RoI', '/', dataNames{i}, '.txt'];
            filename = {'spmT_0001.nii'};
            
            suit_ROI_summarize(filename, 'outfilename', saveDir) %For definitions of RoIs using standard atlas visit http://www.diedrichsenlab.org/imaging/mdtb.htm
            
        end
        
    case 'RSA' %----------------- BEGINNING OF RSA -----------------%
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
            
            cd([subcorticalAreaDir, '/', subj_name{sn}])
            
            matlabbatch{1}.spm.util.imcalc.input = {[subcorticalAreaDir, '/' subj_name{sn} '/' subj_name{sn}, '_aseg.nii,1']};
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
    case 'subcortical_make_mask'  %masks subcortical niftis using grey matter mask
        
        s=varargin{1};
        subcortStructs = subcortStructs(2,:);%names from second row
        
        for i=1:length(subcortStructs)
            funMask=fullfile(glmDir, subj_name{s},'maskbrain.nii');
            omask=fullfile(subcorticalDir, subj_name{s},['mr' subj_name{s} '_' subcortStructs{i} '.nii']); %output mask to be used in the future
            subcort = fullfile(subcorticalDir, subj_name{s},['r' subj_name{s} '_' subcortStructs{i} '.nii']);
            
            spm_imcalc({funMask,subcort},omask,'i1 >0.01 & i2 > 0.1',{}); %recorded activity in brain (grey + white matter)
        end
    case 'subcortical_makeROIs' %uses region toolbox to define R struct for prewhitening
        
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
    case 'subcortical_preWhiten' %pre-whiten data from subcortical ROIs ready for PCM
        
        T = []; blueBear = varargin{1}; %s=varargin{1};
        
        for s=anaSubj
            fprintf('%d.',s); fprintf('/n')
            load([glmDir, '/', subj_name{s}, '/', 'SPM.mat'],            'SPM')
            load([roiSubDir, '/', subj_name{s}, '_subcortical_roi.mat'], 'R')
            cd(fullfile(subcorticalDir, subj_name{s}))
            
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
                
                [betaW,resMS,~,beta] = rsa.spm.noiseNormalizeBeta(Y,SPM);
                
                S.betaW = {betaW};
                S.betaUW = {bsxfun(@rdivide,beta,sqrt(resMS))};
                S.betaRAW = {beta};
                S.resMS = {resMS};
                S.SN = s;
                S.region = r;
                
                T = addstruct(T,S);
                fprintf('%d.',r)
            end
            fprintf('\n');
        end
        
        save(fullfile(roiSubDir, 'preWhitened_betas.mat'),'-struct','T');
    case 'subcortical_calculateRDMs'
        
        cd(roiSubDir)
        R = load('preWhitened_betas.mat');
        nrruns = length(run); nCond = 8;
        saveDir = fullfile(rsaDir, 'subcortical');
        
        %%% prepare condition and partition (run) vectors
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
        
        %%% loop through prewhitened data, calculate crossnobis dissimilarities
        %pre-allocate output variables
        R.d   = NaN(length(R.SN), nCond * (nCond - 1) / 2); %pairwise distance measures
        R.Sig = cell(length(R.SN), 1); %covariance matrix of the beta estimates across different imaging runs.
        R.G   = cell(length(R.SN), 1); %second moment matrix
        R.matrix = cell(length(R.SN), 1); % Pairwise contrast matrix
        
        for s = anaSubj
            load(fullfile(glmDir, subj_name{s}, 'SPM'), 'SPM') %load SPM design matrix for distance function
            for r = unique(R.region)'
                B = R.betaW{R.region == r & R.SN == s};
                
                [d, Sig] = rsa.distanceLDC(B, partVec, condVec, SPM.xX.X);
                [G,~]    = pcm_estGCrossval(B,partVec,condVec, 'X', SPM.xX.X);
                matrix   = indicatorMatrix('allpairs',1:nCond); % Pairwise contrast matrix
                
                R.d     (R.region == r & R.SN == s, :) = d;
                R.Sig   {R.region == r & R.SN == s}    = Sig;
                R.G     {R.region == r & R.SN == s}    = G;
                R.matrix{R.region == r & R.SN == s}    = matrix;
            end
        end
        
        %%% Identify and extract values for overall, order, and timing
        %%% within preparation, production, and cross-phase
        orderDiff      = [0 1 1 0 0 1 1 1 1 0 0 1 1 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0];   %pairwise contrast index (1s are where orders are different)
        timingDiff     = [1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1];   %^ index (1s are where timings are different)
        prepCols       = [1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];   %^ index (1s are contrasts within preparation)
        prodCols       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1];   %^ index (1s are contrasts within production)
        crossPhaseCols = [0 0 0 1 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0];   %^ index (1s are contrasts across phases)
        
        for i=1:length(R.d)
            
            R.dOrderPrep(i,:)  = R.d(i, orderDiff == 1 & timingDiff == 0 & prepCols == 1);
            R.dOrderProd(i,:)  = R.d(i, orderDiff == 1 & timingDiff == 0 & prodCols == 1);
            R.dOrderCross(i,:) = R.d(i, orderDiff == 1 & timingDiff == 0 & crossPhaseCols == 1);
            R.dOrderAll(i,:)   = R.d(i, orderDiff == 1 & timingDiff == 0);
            
            R.dTimingPrep(i,:)  = R.d(i, timingDiff == 1 & orderDiff == 0 & prepCols == 1);
            R.dTimingProd(i,:)  = R.d(i, timingDiff == 1 & orderDiff == 0 & prodCols == 1);
            R.dTimingCross(i,:) = R.d(i, timingDiff == 1 & orderDiff == 0 & crossPhaseCols == 1);
            R.dTimingAll(i,:)   = R.d(i, timingDiff == 1 & orderDiff == 0);
            
            R.dOverallPrep(i,:)  = R.d(i, prepCols == 1);
            R.dOverallProd(i,:)  = R.d(i, prodCols == 1);
            R.dOverallCross(i,:) = R.d(i, crossPhaseCols == 1);
        end
        
        if ~isfolder(saveDir)
            mkdir(saveDir)
        end
        
        save(fullfile(saveDir, 'subRoiDistances.mat'), 'R')
    case 'subcortical_calculateRDMs_integrated'
        %subtracts averaged order and timing patterns within runs
        %subtraction happens during distance function
        %(distanceLDC_cor4main_RY)
        
        cd(roiSubDir)
        R = load('preWhitened_betas.mat');
        nrruns = length(run); nCond = 8;
        saveDir = fullfile(rsaDir, 'subcortical');
        
        %%% prepare condition and partition (run) vectors
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
        
        %%% loop through prewhitened data, calculate crossnobis dissimilarities
        %pre-allocate output variables
        R.d   = NaN(length(R.SN), nCond * (nCond - 1) / 2); %pairwise distance measures
        R.Sig = cell(length(R.SN), 1); %covariance matrix of the beta estimates across different imaging runs.
        R.G   = cell(length(R.SN), 1); %second moment matrix
        R.matrix = cell(length(R.SN), 1); % Pairwise contrast matrix
        
        for s = anaSubj
            load(fullfile(glmDir, subj_name{s}, 'SPM'), 'SPM') %load SPM design matrix for distance function
            for r = unique(R.region)'
                B = R.betaW{R.region == r & R.SN == s};
                
                [d, Sig] = rsa.distanceLDC_cor4main_RY(B, partVec, condVec, SPM.xX.X); %integrated function
                [G,~]    = pcm_estGCrossval(B,partVec,condVec, 'X', SPM.xX.X);
                matrix   = indicatorMatrix('allpairs',1:nCond); % Pairwise contrast matrix
                
                R.d     (R.region == r & R.SN == s, :) = d;
                R.Sig   {R.region == r & R.SN == s}    = Sig;
                R.G     {R.region == r & R.SN == s}    = G;
                R.matrix{R.region == r & R.SN == s}    = matrix;
            end
        end
        
        %%% Identify and extract values for overall, order, and timing
        %%% within preparation, production, and cross-phase
        orderDiff      = [0 1 1 0 0 1 1 1 1 0 0 1 1 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0];   %pairwise contrast index (1s are where orders are different)
        timingDiff     = [1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1];   %^ index (1s are where timings are different)
        prepCols       = [1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];   %^ index (1s are contrasts within preparation)
        prodCols       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1];   %^ index (1s are contrasts within production)
        crossPhaseCols = [0 0 0 1 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0];   %^ index (1s are contrasts across phases)
        
        for i=1:length(R.d)
            
            R.dOrderPrepInt(i,:)  = R.d(i, orderDiff == 1 & timingDiff == 0 & prepCols == 1);
            R.dOrderProdInt(i,:)  = R.d(i, orderDiff == 1 & timingDiff == 0 & prodCols == 1);
            R.dOrderCrossInt(i,:) = R.d(i, orderDiff == 1 & timingDiff == 0 & crossPhaseCols == 1);
            R.dOrderAllInt(i,:)   = R.d(i, orderDiff == 1 & timingDiff == 0);
            
            R.dTimingPrepInt(i,:)  = R.d(i, timingDiff == 1 & orderDiff == 0 & prepCols == 1);
            R.dTimingProdInt(i,:)  = R.d(i, timingDiff == 1 & orderDiff == 0 & prodCols == 1);
            R.dTimingCrossInt(i,:) = R.d(i, timingDiff == 1 & orderDiff == 0 & crossPhaseCols == 1);
            R.dTimingAllInt(i,:)   = R.d(i, timingDiff == 1 & orderDiff == 0);
            
            R.dIntegratedPrep(i,:)  = R.d(i, prepCols == 1);
            R.dIntegratedProd(i,:)  = R.d(i, prodCols == 1);
            R.dIntegratedCross(i,:) = R.d(i, crossPhaseCols == 1);
        end
        
        if ~isfolder(saveDir)
            mkdir(saveDir)
        end
        
        save(fullfile(saveDir, 'subRoiDistances_integrated.mat'), 'R')
    case 'subcortical_plotRDMs'
        subcortStructs = subcortStructs(2,:);%names from second row
        subcortStructs = strrep(subcortStructs, '_', ' '); %replace _ with space
        concatFile = fullfile(rsaDir, 'subcortical', 'subRoiDistancesAll.mat'); %file with all region distances
        
        if ~exist(concatFile, 'file')
            load(fullfile(rsaDir, 'subcortical', 'subRoiDistances.mat'), 'R');
            Rovr = R;
            load(fullfile(rsaDir, 'subcortical', 'subRoiDistances_integrated.mat'), 'R');
            Rint = R;
            
            conds = {...
                'dOverallPrep',    'dOverallProd',    'dOverallCross', ...%variable name
                'dOrderPrep',      'dOrderProd',      'dOrderCross', ...
                'dTimingPrep',     'dTimingProd',     'dTimingCross';...
                1,                 1,                 1, ...%overall/order/timing - cond
                2,                 2,                 2, ...
                3,                 3,                 3; ...
                1,                 2,                 3, ...%prep/prod/cross - phase
                1,                 2,                 3, ...
                1,                 2,                 3, ...
                };
            condsInt = {...
                'dIntegratedPrep', 'dIntegratedProd', 'dIntegratedCross';...
                4,                 4,                 4; ...%integrated - cond
                1,                 2,                 3, ...%prep/prod/cross - phase
                };
            nConds = length(conds);
            nCondsInt = length(condsInt);
            
            A = [];
            
            for i=1:nConds%for factorial conds
                Atemp.SN     = Rovr.SN;
                Atemp.region = Rovr.region;
                Atemp.dist   = mean(Rovr.(conds{1, i}), 2);%mean of relevant pairwise distances
                Atemp.cond   = ones(length(Atemp.dist), 1) * conds{2, i};
                Atemp.phase  = ones(length(Atemp.dist), 1) * conds{3, i};
                
                A = addstruct(A, Atemp);
                clear Atemp
            end%for factorial conds
            for i=1:nCondsInt%for int conds
                Atemp.SN     = Rint.SN;
                Atemp.region = Rint.region;
                Atemp.dist   = mean(Rint.(condsInt{1, i}), 2);%mean of relevant pairwise distances
                Atemp.cond   = ones(length(Atemp.dist), 1) * condsInt{2, i};
                Atemp.phase  = ones(length(Atemp.dist), 1) * condsInt{3, i};
                
                A = addstruct(A, Atemp);
                clear Atemp
            end%for int conds
            
            save(fullfile(rsaDir, 'subcortical', 'subRoiDistancesAll.mat'), 'A');
        else
            load(concatFile, 'A')%if it exists, load it
        end
        
        figure %%% General overview (zoom to regions of interest)
        T = tapply(A,{'SN', 'region', 'cond', 'phase'},{'dist', 'mean', 'name', 'dist'});
        colour={[0 0 0], [0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
        regions = repmat({'l thal', 'l caud', 'l put', 'l pal', 'l hip', 'r thal', 'r caud', 'r put', 'r pal', 'r hip'}, 1, 12);
        barplot([T.phase, T.cond, T.region], T.dist, 'split', T.cond, 'facecolor', colour)
        ylim([-0.0017 0.022])
        set(gca,'xticklabel',regions)
        ylabel('Crossnobis dissimilarity')
        title('Overview')
        %-------------------------------------------------------------------------------%
        
        for i=1:length(subcortStructs)
            figure %%% Region plots for prep/prod order/timing/integrated
            T = tapply(A,{'SN', 'cond', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset', A.region == i & A.cond > 1 & A.phase < 3);
            colour={[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
            lineplot([T.phase], T.dist, 'split', T.cond, ...
                'markertype', 'o', 'markercolor', colour, 'markerfill', colour, 'markersize', 5, ...
                'linecolor', colour, 'linewidth', 3, 'errorwidth', 2, 'errorcolor', colour)
            ylim([-0.0018 0.0055])
            drawline(0, 'dir', 'horz', 'linestyle', '- -')
            ylabel('Crossnobis dissimilarity')
            set(gca,'xticklabel',{'Prep', 'Prod'})
            title(subcortStructs{i})
        end
        %-------------------------------------------------------------------------------%
        
        figure %%% Investigate cross-phase distances
        T = tapply(A,{'SN', 'region', 'cond', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset', A.cond == 1 & A.phase == 3);
        %colour={[0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
        regions = repmat({'l thal', 'l caud', 'l put', 'l pal', 'l hip', 'r thal', 'r caud', 'r put', 'r pal', 'r hip'}, 1, 12);
        barplot([T.phase, T.cond, T.region], T.dist, 'split', T.cond) %'facecolor', colour)
        ylim([0 0.022])
        set(gca,'xticklabel',regions)
        ylabel('Crossnobis dissimilarity')
        title('Cross-phase')
        
        
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
        
    case 'cerebellum_make_nii' %----------------- BEGINNING OF SUIT RSA (CEREBELLUM) -----------------%
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
    case 'cerebellum_make_mask' %restrict area of analysis to grey matter - produces 'maskbrainSUIT.nii'
        
        s=varargin{1};
        
        if isfolder([baseDir '/imaging/suit/' subj_name{s}]) == 0
            mkdir([baseDir '/imaging/suit/' subj_name{s}])
        end
        
        mask=fullfile(glmDir, subj_name{s},'mask.nii');
        %suit=fullfile(anatDir, subj_name{s},['c_', subj_name{s},'_anatomical_pcereb.nii']); %pcereb holds all cerebellum-related regions to a value of 1...
        suit=fullfile(anatDir, subj_name{s},[subj_name{s}, '_anatomical_seg1.nii']); %whereas _seg1 is only grey matter and sets extra-cerebellar regions (e.g. pons) to values other than 1...
        omask=fullfile(suitDir, subj_name{s},'maskbrainSUIT.nii');
        
        spm_imcalc_ui({mask,suit},omask,'i1>0 & i2>0.999',{}); %so including a mask of 0.999 makes sure we only include cerebellar regions.
    case 'cerebellum_make_search'
        
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
    case 'cerebellum_run_search'
        
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
    case 'cerebellum_run_search_integrated' %runs the searchlight, but subtracts averaged order and timing patterns within runs
        %subtraction happens during distance function
        %(distanceLDC_cor4main_RY) within searchlight function
        
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
        end
        
        %%% prepare condition and partition (run) vectors
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 0 0]; %prep
        prep=[repmat(prep,1,nrruns) runBSL];%1 x nBeta, 1 2 3 4 = prep sequences 1:4
        
        prod      =[5 0 0 0 6 0 0 0 7 0 0 0 8 0 0 0 0 0 0 0]; %prod
        prod=[repmat(prod,1,nrruns) runBSL];%1 x nBeta, 5 6 7 8 = prod sequences 1:4
        
        condVec = prep + prod; condVec = condVec'; % conditions, including no interest regressors as 0
        
        %%% Run searchlight function on whole CB
        rsa.runSearchlightLDC_cor4main_RY(L, SPM, 'spmDir', spmDir, 'conditionVec', condVec, ...
            'analysisName', [subj_name{s}, 'RSA_All_integrated'], 'outDir', fullfile(suitDir, subj_name{s}))
    case 'cerebellum_smooth'
        
        s=varargin{1};
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} 'RSA_All_LDC.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} 'RSA_ALL_sLDC.nii']);
        spm_smooth(comb,scomb,[2 2 2]); %smooth with 2mm kernel
        
        comb=fullfile(suitDir, subj_name{s},[subj_name{s} 'RSA_All_integrated_LDC.nii']); %%MVPA smoother
        scomb=fullfile(suitDir, subj_name{s},[subj_name{s} 'RSA_ALL_integrated_sLDC.nii']);
        spm_smooth(comb,scomb,[2 2 2]); %smooth with 2mm kernel
    case 'cerebellum_calc_dissimilarity_maps' %extract each condition (overall, order, timing, etc) from RSA_ALL_sLDC.nii
        %each volume of the searchlight corresponds to a pairwise
        %dissimilarity measure between sequences. So here we average within
        %our conditions to give us dissimilarity measures for:
        %overall prep, overall prod, overall cross, order prep, order prod,
        %order cross, timing prep, timing prod, timing cross
        
        s=varargin{1};
        cd(fullfile(suitDir, subj_name{s}))
        
        %%% Identify and extract values for overall, order, and timing
        %%% within preparation, production, and cross-phase
        ordDiff        = [0 1 1 0 0 1 1 1 1 0 0 1 1 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0]';   %pairwise contrast index (1s are where orders are different)
        timDiff        = [1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1]';   %^ index (1s are where timings are different)
        prepCols       = [1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';   %^ index (1s are contrasts within preparation)
        prodCols       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1]';   %^ index (1s are contrasts within production)
        crossPhaseCols = [0 0 0 1 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0]';   %^ index (1s are contrasts across phases)
        
        conds = {...
            'overall_prep', 'overall_prod', 'overall_cross', ...
            'order_prep',                                'order_prod',                                'order_cross', ...
            'timing_prep',                               'timing_prod',                               'timing_cross';...
            prepCols == 1,  prodCols == 1,  crossPhaseCols == 1, ...
            prepCols == 1 & ordDiff == 1 & timDiff == 0, prodCols == 1 & ordDiff == 1 & timDiff == 0, crossPhaseCols == 1 & ordDiff == 1 & timDiff == 0, ...
            prepCols == 1 & timDiff == 1 & ordDiff == 0, prodCols == 1 & timDiff == 1 & ordDiff == 0, crossPhaseCols == 1 & timDiff == 1 & ordDiff == 0 ...
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
    case 'cerebellum_calc_dissimilarity_maps_integrated' %extract integrated prep, prod, cross from RSA_ALL_integrated_sLDC.nii
        %each volume of the searchlight corresponds to a pairwise
        %dissimilarity measure between sequences. So here we average within
        %our conditions to give us dissimilarity measures for:
        %Int prep, Int prod, int cross, order prep, order prod,
        %order cross, timing prep, timing prod, timing cross
        
        s=varargin{1};
        cd(fullfile(suitDir, subj_name{s}))
        
        %%% Identify and extract values for overall, order, and timing
        %%% within preparation, production, and cross-phase
        ordDiff        = [0 1 1 0 0 1 1 1 1 0 0 1 1 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0]';   %pairwise contrast index (1s are where orders are different)
        timDiff        = [1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1]';   %^ index (1s are where timings are different)
        prepCols       = [1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';   %^ index (1s are contrasts within preparation)
        prodCols       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1]';   %^ index (1s are contrasts within production)
        crossPhaseCols = [0 0 0 1 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0]';   %^ index (1s are contrasts across phases)
        
        conds = {...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
            'orderInt_prep',                                'orderInt_prod',                                'orderInt_cross', ...
            'timingInt_prep',                               'timingInt_prod',                               'timingInt_cross';...
            prepCols == 1,     prodCols == 1,     crossPhaseCols == 1, ...
            prepCols == 1 & ordDiff == 1 & timDiff == 0,    prodCols == 1 & ordDiff == 1 & timDiff == 0,    crossPhaseCols == 1 & ordDiff == 1 & timDiff == 0, ...
            prepCols == 1 & timDiff == 1 & ordDiff == 0,    prodCols == 1 & timDiff == 1 & ordDiff == 0,    crossPhaseCols == 1 & timDiff == 1 & ordDiff == 0 ...
            };
        
        for j=1:length(conds)
            vol = spm_vol([subj_name{s} 'RSA_ALL_integrated_sLDC.nii']);
            
            Vi = vol(conds{2,j}, :);
            Vo = Vi(1); Vo = rmfield(Vo, 'pinfo');
            Vo.fname = [subj_name{s} '_LDC_' conds{1,j} '.nii'];
            Vo.n = [1 1];
            express = 'mean(X)';
            flags.dmtx = 1;
            
            spm_imcalc(Vi, Vo, express, flags)
        end
    case 'cerebellum_normalise'
        
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
            'order_prep',      'order_prod',      'order_cross', ...
            'timing_prep',     'timing_prod',     'timing_cross'...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
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
        
        %% prepare struct for function
        job.subj.affineTr = affine; %fill job.subj. struct with respective items
        job.subj.flowfield = flowfield;
        job.subj.resample = dataFiles;
        job.subj.mask = mask;
        job.subj.outname = outFiles;
        
        %function
        suit_reslice_dartel(job)
    case 'cerebellum_group_avg'
        
        if ~isfolder(fullfile(suitGroupRSADir, 'average'))
            mkdir(fullfile(suitGroupRSADir, 'average'))
        end
        cd(fullfile(suitGroupRSADir, 'average'))
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            'order_prep',      'order_prod',      'order_cross', ...
            'timing_prep',     'timing_prod',     'timing_cross'...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
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
    case 'cerebellum_group_randomeffects'
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            'order_prep',      'order_prod',      'order_cross', ...
            'timing_prep',     'timing_prod',     'timing_cross'...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
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
    case 'cerebellum_group_estimate'
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            'order_prep',      'order_prod',      'order_cross', ...
            'timing_prep',     'timing_prod',     'timing_cross'...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
            };
        contrastN = length(conds);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile (suitGroupRSADir, ['RSA_' conds{i}], 'SPM.mat');  %Adjust directory
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
    case 'cerebellum_plot_flatmap'
        
        cd(suitGroupRSADir)
        
        conds = {...
            'overall_prep',    'overall_prod',    'overall_cross', ...
            'order_prep',      'order_prod',      'order_cross', ...
            'timing_prep',     'timing_prod',     'timing_cross'...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
            };
        folderNames = strcat('RSA_', conds);
        
        fileName = fullfile(folderNames{4}, 'spmT_0001.nii');
        surfMap = suit_map2surf(fileName);
        figure
        suit_plotflatmap(surfMap, 'threshold', 3.48, 'cscale', [0, 3.5])
        
        clear fileName
        
        fileName{1} = fullfile('RSA_timing_prep', 'spmT_0001.nii');
        fileName{2} = fullfile('RSA_integrated_prep', 'spmT_0001.nii');
        fileName{3} = fullfile('RSA_order_prep', 'spmT_0001.nii');
        surfMap(:,1) = suit_map2surf(fileName{1});
        surfMap(:,2) = suit_map2surf(fileName{2});
        surfMap(:,3) = suit_map2surf(fileName{3});
        
        figure
        suit_plotflatmap(surfMap, 'type', 'rgb', 'threshold', 10, 'cscale', [0, 3.5])
        
        fileName{1} = fullfile('RSA_timing_prod', 'spmT_0001.nii');
        fileName{2} = fullfile('RSA_integrated_prod', 'spmT_0001.nii');
        fileName{3} = fullfile('RSA_order_prod', 'spmT_0001.nii');
        surfMap(:,1) = suit_map2surf(fileName{1});
        surfMap(:,2) = suit_map2surf(fileName{2});
        surfMap(:,3) = suit_map2surf(fileName{3});
        
        figure
        suit_plotflatmap(surfMap, 'type', 'rgb', 'threshold', 10, 'cscale', [0, 3.5])
        
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
    case 'cortical_run_search_integrated' %runs the searchlight, but subtracts averaged order and timing patterns within runs
        %subtraction happens during distance function
        %(distanceLDC_cor4main_RY) within searchlight function
        
        s=varargin{1}; blueBear=varargin{2};
        nrruns = length(run);
        
        cd(fullfile(rsaDir, 'cortical', subj_name{s}))
        
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
        end
        
        %%% prepare condition and partition (run) vectors
        runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
        prep      =[0 0 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 0 0]; %prep
        prep=[repmat(prep,1,nrruns) runBSL];%1 x nBeta, 1 2 3 4 = prep sequences 1:4
        
        prod      =[5 0 0 0 6 0 0 0 7 0 0 0 8 0 0 0 0 0 0 0]; %prod
        prod=[repmat(prod,1,nrruns) runBSL];%1 x nBeta, 5 6 7 8 = prod sequences 1:4
        
        condVec = prep + prod; condVec = condVec'; % conditions, including no interest regressors as 0
        
        %%% Run searchlight function on whole CB
        rsa.runSearchlightLDC_cor4main_RY(L, SPM, 'spmDir', spmDir, 'conditionVec', condVec, ...
            'analysisName', [subj_name{s}, 'RSA_All_integrated'], 'outDir', fullfile(rsaDir, 'cortical', subj_name{s}))
    case 'cortical_smooth'
        
        s=varargin{1};
        rsaCorticalDir = fullfile(rsaDir, 'cortical');
        
        comb=fullfile(rsaCorticalDir, subj_name{s},[subj_name{s} 'RSA_All_LDC.nii']); %%MVPA smoother
        scomb=fullfile(rsaCorticalDir, subj_name{s},[subj_name{s} 'RSA_ALL_sLDC.nii']);
        spm_smooth(comb,scomb,[2 2 2]); %smooth with 2mm kernel
        
        comb=fullfile(rsaCorticalDir, subj_name{s},[subj_name{s} 'RSA_All_integrated_LDC.nii']); %%MVPA smoother
        scomb=fullfile(rsaCorticalDir, subj_name{s},[subj_name{s} 'RSA_ALL_integrated_sLDC.nii']);
        spm_smooth(comb,scomb,[2 2 2]); %smooth with 2mm kernel
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
        ordDiff        = [0 1 1 0 0 1 1 1 1 0 0 1 1 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0]';   %pairwise contrast index (1s are where orders are different)
        timDiff        = [1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1]';   %^ index (1s are where timings are different)
        prepCols       = [1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';   %^ index (1s are contrasts within preparation)
        prodCols       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1]';   %^ index (1s are contrasts within production)
        crossPhaseCols = [0 0 0 1 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0]';   %^ index (1s are contrasts across phases)
        
        conds = {...
            'overall_prep', 'overall_prod', 'overall_cross', ...
            'order_prep',                                'order_prod',                                'order_cross', ...
            'timing_prep',                               'timing_prod',                               'timing_cross';...
            prepCols == 1,  prodCols == 1,  crossPhaseCols == 1, ...
            prepCols == 1 & ordDiff == 1 & timDiff == 0, prodCols == 1 & ordDiff == 1 & timDiff == 0, crossPhaseCols == 1 & ordDiff == 1 & timDiff == 0, ...
            prepCols == 1 & timDiff == 1 & ordDiff == 0, prodCols == 1 & timDiff == 1 & ordDiff == 0, crossPhaseCols == 1 & timDiff == 1 & ordDiff == 0 ...
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
    case 'cortical_calc_dissimilarity_maps_integrated' %extract integrated prep, prod, cross from RSA_ALL_integrated_sLDC.nii
        %each volume of the searchlight corresponds to a pairwise
        %dissimilarity measure between sequences. So here we average within
        %our conditions to give us dissimilarity measures for:
        %Int prep, Int prod, int cross, order prep, order prod,
        %order cross, timing prep, timing prod, timing cross
        
        s=varargin{1};
        cd(fullfile(rsaCorticalDir, subj_name{s}))
        
        %%% Identify and extract values for overall, order, and timing
        %%% within preparation, production, and cross-phase
        ordDiff        = [0 1 1 0 0 1 1 1 1 0 0 1 1 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0]';   %pairwise contrast index (1s are where orders are different)
        timDiff        = [1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1]';   %^ index (1s are where timings are different)
        prepCols       = [1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';   %^ index (1s are contrasts within preparation)
        prodCols       = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1]';   %^ index (1s are contrasts within production)
        crossPhaseCols = [0 0 0 1 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0]';   %^ index (1s are contrasts across phases)
        
        conds = {...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
            'orderInt_prep',                                'orderInt_prod',                                'orderInt_cross', ...
            'timingInt_prep',                               'timingInt_prod',                               'timingInt_cross';...
            prepCols == 1,     prodCols == 1,     crossPhaseCols == 1, ...
            prepCols == 1 & ordDiff == 1 & timDiff == 0,    prodCols == 1 & ordDiff == 1 & timDiff == 0,    crossPhaseCols == 1 & ordDiff == 1 & timDiff == 0, ...
            prepCols == 1 & timDiff == 1 & ordDiff == 0,    prodCols == 1 & timDiff == 1 & ordDiff == 0,    crossPhaseCols == 1 & timDiff == 1 & ordDiff == 0 ...
            };
        
        for j=1:length(conds)
            vol = spm_vol([subj_name{s} 'RSA_ALL_integrated_sLDC.nii']);
            
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
            'order_prep',      'order_prod',      'order_cross', ...
            'timing_prep',     'timing_prod',     'timing_cross'...
            'integrated_prep', 'integrated_prod', 'integrated_cross'...
            }; % please add other images as required, e.g. spmT_...
        
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
            'order_prep',      'order_prod',      'order_cross', ...
            'timing_prep',     'timing_prod',     'timing_cross'...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
            };
        
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
            'order_prep',      'order_prod',      'order_cross', ...
            'timing_prep',     'timing_prod',     'timing_cross'...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
            };
        
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
            'order_prep',      'order_prod',      'order_cross', ...
            'timing_prep',     'timing_prod',     'timing_cross'...
            'integrated_prep', 'integrated_prod', 'integrated_cross', ...
            };
        contrastN = length(conds);
        
        for i=1:contrastN
            matlabbatch{1}.spm.stats.fmri_est.spmmat{1} = fullfile (rsaCorticalGroupDir, ['RSA_' conds{i}], 'SPM.mat');  %Adjust directory
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
        end
        
    case 'simulations_RSA_LDA'
        
        %%%Prepare variables
        prepProd = 1;
        vord     = [0.3 0.1];
        vtemp    = [0.4 0.7];
        vinter   = [0.6 0.8];
        vnoise   = 0.5;
        nIter    = 500;%acts as subj number
        
        vararginoptions(varargin,{'prepProd','vord','vtemp','vinter','vnoise','nIter','sn'});
        A=[]; %data struct
        B=[]; %output struct - stores distances, LDA, G matrix
        D=[]; %raw data store
        
        %%%Prepare condition and partition (run) vectors
        nrruns = 6; nCond = 8;
        runC      = repmat(1:nrruns,nCond,1); %labelling runs 1:6 (imaging runs)
        runC      = reshape(runC,1,[])';
        condC     = repmat((1:nCond)', nrruns, 1); %labelling conditions 1:4 (sequences)
        condCOrd  = repmat([1 1 2 2 3 3 4 4]', 6, 1);
        condCTemp = repmat([1 2 1 2 3 4 3 4]', 6, 1);
        
        %%%Generate data
        for i=1:nIter %for 'subject' n
            %make the data based on input specifications
            Y=prepProdSimu_makedataPP('prepProd',prepProd,'vord',vord,'vtemp',vtemp,'vinter',vinter,'vnoise',vnoise);
            
            A.data{i} = Y;
            
            %%%Generate distances vector
            [A.d{i},     A.Sig{i}]     = rsa.distanceLDC(A.data{i}, runC, condC);
            [A.dOrd{i},  A.SigOrd{i}]  = rsa.distanceLDC_oneOut_RY(A.data{i}, runC, condC, condCOrd);
            [A.dTemp{i}, A.SigTemp{i}] = rsa.distanceLDC_oneOut_RY(A.data{i}, runC, condC, condCTemp);
            [A.dInt{i},  A.SigInt{i}]  = rsa.distanceLDC_cor4main_RY(A.data{i}, runC, condC);
            
            %%%Run linear decoding
            prepData = [A.data{i}(1:4,:); A.data{i}(9:12,:);  A.data{i}(17:20,:); A.data{i}(25:28,:); A.data{i}(33:36,:); A.data{i}(41:44,:)];
            prodData = [A.data{i}(5:8,:); A.data{i}(13:16,:); A.data{i}(21:24,:); A.data{i}(29:32,:); A.data{i}(37:40,:); A.data{i}(45:48,:)];
            
            %               accuracy columns
            %1: overall,        2: order,          3: timing,         4:integrated(Z)    5: integrated(acc) for debugging
            [prepAccuracy(i,1), prepAccuracy(i,2), prepAccuracy(i,3), prepAccuracy(i,4), prepAccuracy(i,5)]=prepProdSimu_classify(prepData);
            [prodAccuracy(i,1), prodAccuracy(i,2), prodAccuracy(i,3), prodAccuracy(i,4), prodAccuracy(i,5)]=prepProdSimu_classify(prodData);
        end%for subject n
        
        % Identify and extract values for overall, order, and timing
        % within preparation, production, and cross-phase
        diffVectorOvr  = [... %index (1s are where condition is different)
            1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1; ...%overall
            ];
        %0 1 1 0 0 1 1 1 1 0 0 1 1 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0; ...%order
        %1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 1 0 1; ...%timing
        %0 0 1 0 0 0 1 1 0 0 0 1 0 0 0 1 0 0 1 0 0 0 0 0 1 1 0 0; ...%integrated
        
        %diffVector  = [...
        %diffVector(1,:); ...%overall stays the same
        %diffVector(2,:) == 1 & diffVector(3,:) == 0; ...%ONLY sequences with diff order & SAME timing
        %diffVector(3,:) == 1 & diffVector(2,:) == 0; ...%vice versa as above
        %diffVector(4,:); ...%integrated stays the same
        %];
        %diffVectorInt = [...
            %1 1 0 0 1 1 0 0 1 1 0 0 1 1 1 0 0 1 0 1 1 0 1 1 0 0 1 1; ...%integrated
            %];
        phaseVector = [...%index (1s are where phase is different)
            1 1 1 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...%prep
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1; ...%prod
            0 0 0 1 1 1 1 0 0 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0  ...%cross
            ];
        phaseVectorOneOut = [...
            1 0 0 0 0 0; ...%prep
            0 0 0 0 0 1; ...%prod
            0 1 1 1 1 0; ...%cross
            ];
        
        %%%Loop over subjs, conditions, and phases, and store distances as
        %%%struct to plot later
        loopCounter = 1;
        
        %%% Overall
        for s=1:nIter%for subj
            for i=1:size(phaseVector, 1)%for phase
                for j=1:size(diffVectorOvr, 1)%for conds
                    B.dist(loopCounter,1)  = mean(A.d{s}(phaseVector(i,:) & diffVectorOvr(j,:)));
                    B.phase(loopCounter,1) = i;
                    B.cond(loopCounter,1)  = j;
                    B.sn(loopCounter,1)    = s;
                    loopCounter = loopCounter + 1;
                end%for conds
            end%for phase
        end%for subj
        
        %%% Order and timing
        for s=1:nIter%for subj
            for i=1:size(phaseVectorOneOut, 1)%for phase
                    B.dist(loopCounter,1)  = mean(A.dOrd{s}(phaseVectorOneOut(i,:) == 1));
                    B.phase(loopCounter,1) = i;
                    B.cond(loopCounter,1)  = 2;
                    B.sn(loopCounter,1)    = s;
                    loopCounter = loopCounter + 1;
                    B.dist(loopCounter,1)  = mean(A.dTemp{s}(phaseVectorOneOut(i,:) == 1));
                    B.phase(loopCounter,1) = i;
                    B.cond(loopCounter,1)  = 3;
                    B.sn(loopCounter,1)    = s;
                    loopCounter = loopCounter + 1;
            end%for phase
        end%for subj
        
        %%% Integrated 
        for s=1:nIter%for subj
            for i=1:size(phaseVector, 1)%for phase
                B.dist(loopCounter,1)  = mean(A.dInt{s}(phaseVector(i,:) & diffVectorOvr(1,:)));
                B.phase(loopCounter,1) = i;
                B.cond(loopCounter,1)  = 4;
                B.sn(loopCounter,1)    = s;
                loopCounter = loopCounter + 1;
            end%for phase
        end%for subj
        
        loopCounter = 1;
        for s=1:nIter%for subj
            for i=1:size(prepAccuracy, 2)
                C.acc(loopCounter,1)   = prepAccuracy(s, i);
                C.phase(loopCounter,1) = 1;
                C.cond(loopCounter,1)  = i;
                C.sn(loopCounter,1)    = s;
                loopCounter = loopCounter + 1;
            end
        end%for subj
        
        for s=1:nIter%for subj
            for i=1:size(prodAccuracy, 2)
                C.acc(loopCounter,1)   = prodAccuracy(s, i);
                C.phase(loopCounter,1) = 2;
                C.cond(loopCounter,1)  = i;
                C.sn(loopCounter,1)    = s;
                loopCounter = loopCounter + 1;
            end
        end%for subj
        
        
        
        figure %%% General overview (zoom to regions of interest)
        T = tapply(B,{'sn', 'cond', 'phase'},{'dist', 'mean', 'name', 'dist'}, 'subset', B.phase < 3);
        colour={[0 0 0], [0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
        barplot([T.phase, T.cond], T.dist, 'split', T.cond, 'facecolor', colour)
        ylabel('Crossnobis dissimilarity')
        title(['prepProd =', num2str(prepProd), '  iter=', num2str(nIter), '  ORD=',num2str(vord), '  TEMP=',num2str(vtemp), '  INTER=',num2str(vinter), '  NOISE=',num2str(vnoise)]);
        %-------------------------------------------------------------------------------%
        
        figure %%% General overview (zoom to regions of interest)
        T = tapply(C,{'sn', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'});
        colour={[0 0 0], [0 0.4470 0.7410], [0.6350 0.0780 0.1840], [0.4660 0.6740 0.1880]};
        barplot([T.phase, T.cond], T.acc, 'split', T.cond, 'facecolor', colour)
        ylabel('Decoding accuracy (Z)')
        title(['prepProd =', num2str(prepProd), '  iter=', num2str(nIter), '  ORD=',num2str(vord), '  TEMP=',num2str(vtemp), '  INTER=',num2str(vinter), '  NOISE=',num2str(vnoise)]);
        %-------------------------------------------------------------------------------%
        
        labels = {'O1T1p', 'O1T2p', 'O2T1p', 'O2T2p', 'O1T1P', 'O1T2P', 'O2T1P', 'O2T2P'};%p=prep, P=prod
        
        %%%Multi-dimensional scaling - RDM
        for s=1:length(A.data)%for iter
            d_hat(:,:,s)= squareform(A.d{s}); %calc representational dissimilarity matrix
        end%forIter
        dm = mean(d_hat,3); % Mean estimate
        
        figure
        %subplot(1,2,1);
        H = eye(8)-ones(8)/8;
        imagesc(H*dm*H');
        title(['prepProd =', num2str(prepProd), '  iter=', num2str(nIter), '  ORD=',num2str(vord), '  TEMP=',num2str(vtemp), '  INTER=',num2str(vinter), '  NOISE=',num2str(vnoise)]);
        xticks(1:8)
        xticklabels(labels);%p=prep, P=prod
        yticks(1:8)
        yticklabels(labels)
        
        % %%% Multi-dimensional scaling
        try
            RDMs.RDM   = dm;
            RDMs.name  = 'simulation';
            RDMs.color = [0 0 1];
            
            userOptions = prepProdSimuPP_defineUserOptions(baseDir);
            
            rsa.MDSConditions(RDMs, userOptions);
        catch
            disp('MDS Failed.')
        end
        %----------------------------------------------------------------------------------------------%
        
        %%%Multi-dimensional scaling - RDM
        for s=1:length(A.data)%for iter
            d_hatInt(:,:,s)= squareform(A.dInt{s}); %calc representational dissimilarity matrix
        end%forIter
        dmInt = mean(d_hatInt,3); % Mean estimate
        
        figure
        %subplot(1,2,1);
        H = eye(8)-ones(8)/8;
        imagesc(H*dmInt*H');
        title('Corrected for main (Integrated)');
        xticks(1:8)
        xticklabels(labels);%p=prep, P=prod
        yticks(1:8)
        yticklabels(labels)
        
        % %%% Multi-dimensional scaling
        try
            RDMsInt.RDM   = dmInt;
            RDMsInt.name  = 'simulation_integrated';
            RDMsInt.color = [0 0 1];
            
            userOptions = prepProdSimuPP_defineUserOptions;
            
            rsa.MDSConditions(RDMsInt, userOptions);
        catch
            disp('MDS Failed.')
        end
        
        T = tapply(C,{'sn', 'cond', 'phase'},{'acc', 'mean', 'name', 'acc'}, 'subset', C.phase == 1);
        barplot(T.cond,T.acc,'split',T.cond,'style_bold','leg',{'ovr', 'ord', 'temp' ,'inter'},'leglocation','north');
        title(['prepProd =', num2str(prepProd), '  iter=', num2str(nIter), '  ORD=',num2str(vord), '  TEMP=',num2str(vtemp), '  INTER=',num2str(vinter), '  NOISE=',num2str(vnoise)]);
        ylabel('acc');
        xlabel('factor');
        
end%switch
end%function