function prepProd2_imana_RY(what,varargin)


%% Load relevant Matlab directories
% addpath(genpath('D:\projects\rhys\prepProd2\matlab')); %Ajust! loaded with subdirectories (genpath command)
% addpath(genpath('D:\projects\toolboxes\spm12'));
% addpath(genpath('D:\projects\toolboxes\tools')); %joern's extensions for spm
% addpath(genpath('D:\projects\toolboxes\userfun')); %joern's util tools (open source)
% addpath(genpath('D:\projects\toolboxes\region')); %joern's region toolbox for spm

%%%definition and variables

%% Data paths:
baseDir= 'D:\projects\rhys\prepProd2\data';
baseDirR= 'D:\projects\rhys\prepProd2\data';
rawDir=[baseDirR filesep 'imaging' filesep 'raw']; %original files before conversion
anatDir=[baseDir filesep 'imaging' filesep 'anatomicals'];
epiDir=[baseDir filesep 'imaging' filesep 'epi'];
glmDir=[baseDir filesep 'imaging' filesep 'GLM_firstlevel'];
behDir=[baseDir filesep 'behavioural'];
groupDir=[baseDir filesep 'imaging' filesep 'GLM_secondlevel'];
scndDir=[baseDir filesep 'imaging' filesep 'GLM_secondlevel' filesep 'data'];
regDir=[baseDir filesep 'imaging' filesep 'RegionOfInterest'];




%% Names and IDs
subj=[2,3,4,5,6,7,8,9,10,11];
subj_name={'s02','s03','s04','s05','s06','s07','s08','s09','s10','s11'};
subjID={'SS0094_09','SS0094_10','SS0094_11','SS0094_14','SS0094_13','SS0094_12','SS0094_15','SS0094_16','SS0094_17','SS0094_18'}; %folder

subjIDsession={'SS0094_09_20201113_631887115','SS0094_10_20201119_632400213','SS0094_11_20201119_632410210','SS0094_14_20201126_633022974',...
    'SS0094_13_20201126_633008222','SS0094_12_20201126_632998125','SS0094_15_20201203_631017618','PrepProd2_s09','PrepProd2_s10','PrepProd2_s11'}; %beginning of filename raw files

scanID={'301','401','501','601','701','801'; ...
    '401','501','601','701','801','901'; ...
    '301','401','501','601','701','901'; ...
    '301','401','501','601','701','801'; ...
    '301','401','501','601','701','801'; ...
    '301','501','701','801','901','1001'; ...
    '301','401','501','601','701','801'; ...
    '301','401','501','601','701','801'; ...
    '301','401','501','601','701','801'; ...
    '301','401','501','601','701','801'};

fieldMapID={};
fieldMapID2={};
fieldMapName={'magnitude','phase'};
runID={'R1','R2','R3','R4','R5','R6'};
run={'1','2','3','4','5','6'};
dataext='run';
hemName={'LeftHem','RightHem'};
hem={'lh','rh'};
atlasA={'i','x'};
atlas= 2;
atlasname={'fsaverage','fsaverage_sym'};


%% Parameters
delay=1; %SPM starts with TR 0; designfiles start counting with 1; -> has to be substracted out CHECK!!!
TR=2;
nrVolumes=230; %number of volumes
nrSlices=60; %different across subjects
radius=16;
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
originAC=[127,-131,-72 1; ...
    123,-129,-76 1; ...
    119,-131,-71 1; ...
    123,-135,-75 1; ...
    122,-131,-76 1; ...
    121,-129,-72 1; ...
    122,-138,-83 1; ...
    125,-141,-72 1; ...
    123,-130,-75 1; ...
    130,-138,-70 1];


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
            source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} '_202_scombi.nii.gz']);
            dest = fullfile(rawDir, subjID{sn});
            gunzip(source,dest); %unzip
            
            source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} '_202_scombi.nii']);
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
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn}  '_' scanID{sn,i} '_FE_EPI_2_iso_MB2' '.nii.gz']);
                dest = fullfile(rawDir, subjID{sn});
                gunzip(source,dest); %unzip
                
                source = fullfile(rawDir, subjID{sn}, [subjIDsession{sn} '_' scanID{sn,i} '_FE_EPI_2_iso_MB2' '.nii']);
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
        for r= 1:length(run);
            for i=1:nrVolumes
                N{i,1} = [fullfile(epiDir, subj_name{sn}, [ subj_name{sn},'_run',run{r},'.nii,',num2str(i)])];
            end;
            J.scans{r} = N;
        end
        matlabbatch{1}.spm.temporal.st.scans = J.scans;
        matlabbatch{1}.spm.temporal.st.nslices = nrSlices;
        matlabbatch{1}.spm.temporal.st.tr = TR;
        matlabbatch{1}.spm.temporal.st.ta = TR-(TR/nrSlices);
        matlabbatch{1}.spm.temporal.st.so = (sliceAcquisition);
        matlabbatch{1}.spm.temporal.st.refslice = 1;
        matlabbatch{1}.spm.temporal.st.prefix = 'a';
        spm_jobman('run',matlabbatch);
        
        
    case 'realign_unwarp' %based on JD's code; extracts motion regressors for glm    http://www.diedrichsenlab.org/imaging/robustWLS.html
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
        
        %%% *1. Coregister meanEPI to anatomical (no reslicing into anat space, dimension of EPI preserved!): %%%
        %NOTE: If original image completely off in terms of alignment to anat - First coregister per hand via 'coregtool' command
        %%then run the spm algorithm:
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
        
                %%% *4. Check whether runs are really aligned to the rmeanepi on
                %%% MRIcron: using anatomical and ua* data overlap
        
%         spm_jobman('run',matlabbatch);
        
    case 'glm_set'
        sn=varargin{1};
        prefix='ua'; %% if bias correction was performed previously
        delay=0;
        dur=nrVolumes;  % TRs per Block
        fMRIblockNr=(47:52); %BN behavioural file
        T=[];
        
        for r=1:numel(run)
            
            %Load behavioural file (structure B)
            fname=fullfile(behDir,subj_name{sn},['exp_BN' num2str(fMRIblockNr(r)) '.mat']);
            load(fname);
            
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
                J.sess(r).cond(c).onset = [B.tZeroCue(idx)/1000];
                J.sess(r).cond(c).duration = max(B.tZero- B.tZeroCue)/1000; %max prep time in sec; %[(B.tZero(idx) - B.tZeroCue(idx))/1000];
                J.sess(r).cond(c).tmod = 0;
                J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
            end;
            
            
            %%% Sequence initiation/1st press
            for c=1:2:8                  % Loop over the conditions (production trials only!!) conf: 1     3     5     7
                idx=find(B.cond==c);
                
                % Make new structure(S) that reflects what each of the
                % regressors in the design matrix -> goes into T
                S.RN(c,1)=r; %run n
                S.SN(c,1)=sn; %subject nr
                S.condition(c,1)=c;
                S.numtrials(c,1)=length(idx);
                
                % Info on condition name, onset and duration
                J.sess(r).cond(c).name = sprintf('%d_%d',r,c);
                J.sess(r).cond(c).onset = [(B.tZero(idx)+ B.timing(idx,1))/1000];
                J.sess(r).cond(c).duration = [(max(B.timing(idx),[],2)-B.timing(idx,1))/1000];
                J.sess(r).cond(c).tmod = 0;
                J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
            end;
            
            
            %%% Error (feedback) %work in progress
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
            J.sess(r).cond(cond).onset = B.tZeroCue(idx)/1000;
            J.sess(r).cond(cond).duration = [(B.tZero(idx) - B.tZeroCue(idx))/1000]; %Time between Sequence and Go cues
            J.sess(r).cond(cond).tmod = 0;
            J.sess(r).cond(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
            %end;
            c=cond;
            
            
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
        rest =[-1 0 -1 0 -1 0 -1 0 -1 0 -1 0 -1 0 -1 0 -1 0  -1 0]; %rest
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
        
        s=varargin{1};
        con=fullfile(glmDir, subj_name{s},'con_0002.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0002.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        con=fullfile(glmDir, subj_name{s},'con_0003.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0003.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        con=fullfile(glmDir, subj_name{s},'con_0004.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0004.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        con=fullfile(glmDir, subj_name{s},'con_0005.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0005.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        con=fullfile(glmDir, subj_name{s},'con_0006.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} '_scon_0006.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        con=fullfile(glmDir, subj_name{s},'spmT_0001.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} 'SspmT_0001.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        s=varargin{1};
        con=fullfile(glmDir, subj_name{s},'spmT_0002.nii'); %%contrast smoother
        scon=fullfile(glmDir, subj_name{s},[subj_name{s} 'SspmT_0002.nii']);
        spm_smooth(con,scon,[4 4 4]); %smooth with 4mm kernel
        
        %smooth other images here as required
        
    case 'glm_contrastGroup'

        %Produces second-level contrasts. Edit contrast folders, image names, and subject nifti files.
        
        dataDir = {'Mov', 'Prep', 'Error', 'PrepProd', 'ProdPrep', 'Rest'}; %%Save folders for each contrast
        images = {'scon_0001';'scon_0002';'scon_0003';'scon_0004';'scon_0005';'scon_0006'};
%         subNii = {'_s01.nii','_s02.nii','_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s08.nii','_s09.nii','_s10.nii'};
        subNii = {'_s05.nii','_s06.nii','_s07.nii'};
        
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
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@prepProd_combinedclass_corrected4Main,'params',{c,run,train,test});
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
        lmva_spm(fullfile(glmDir,subj_name{s},'volsearch160.mat'),Pselect,out,@prepProd_combinedclass_corrected4Main,'params',{c,run,train,test});
        telapsed = toc(tstart)
        
    case 'MVA_group'
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
        images = {'szacc_Comb_160_Mov';'szacc_Comb_160_Prep';'szacc_Int_160_Mov';'szacc_Int_160_Prep';'szacc_Spat_160_Mov';'szacc_Spat_160_Prep';'szacc_Temp_160_Mov';'szacc_Temp_160_Prep'};
        subNii = {'_s01.nii','_s02.nii','_s03.nii','_s05.nii','_s06.nii','_s07.nii','_s08.nii','_s09.nii','_s10.nii'};
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
        
    case 'MVA_estimate'
        dataDir = {'MVA_comb_mov', 'MVA_comb_prep', 'MVA_int_mov', 'MVA_int_prep', 'MVA_spat_mov', 'MVA_spat_prep','MVA_temp_mov','MVA_temp_prep'}; %%Save folders for each contrast
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
        
        outimages={'_zacc_Comb_160_Mov','_zacc_Comb_160_Prep','_zacc_Int_160_Mov','_zacc_Int_160_Prep'};
        
        
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
        
    case 'MVA_smooth' %%%Smoothing before suit reslice %%Subject space
        
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
        %         images= {'szacc_Comb_160_Mov.nii','szacc_Comb_160_Prep.nii','_saccuracy_Comb_160_Mov.nii','_saccuracy_Comb_160_Prep.nii',...
        %             'szacc_Spat_160_Mov','szacc_Spat_160_Prep','szacc_Temp_160_Mov','szacc_Temp_160_Prep','szacc_Int_160_Mov','szacc_Int_160_Prep'...
        %             'scon_0001.nii','scon_0002.nii','scon_0003.nii','scon_0004.nii','scon_0005.nii','scon_0006.nii'}; % please add other images as required, e.g. spmT_...
        
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
                images= {'scon_0001.nii','scon_0002.nii','scon_0003.nii','scon_0004.nii','scon_0005.nii','scon_0006.nii'}; % please add other images as required, e.g. spmT_...
        
                %          images= {'szacc_Comb_160_Prep.nii'}; % please add other images as required, e.g. spmT_...
        
                s=varargin{1};
                defor= fullfile(anatDir, subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
                for j=1:numel(images)
                    [dir,name,ext]=spm_fileparts(images{j});
                    sn_images{j}= fullfile(glmDir,subj_name{s},[images{j}]);
                    out_images{j}= fullfile(groupDir,[name '_' subj_name{s} '.nii']);
                end
                spmj_normalization_write(defor, sn_images,'outimages',out_images); %Trilinear interpolation
        
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
        roiOUT='D:\projects\rhys\prepProd\data\imaging\ROI\roi_Elife2014_bilateral.mat';
        
        
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
        roiIN='D:\projects\rhys\prepProd\data\imaging\ROI\roi_Elife2014_bilateral.mat';
        
        load(roiIN);
        P={};
        %%%% Comb
        cd('D:\projects\rhys\prepProd\data\imaging\ROI');
        MVA_ID={'Comb','Temp','Spat','Int'};%%%TODO
        Phase_ID={'Prep','Mov'};
        ROI=[];
        K=[];
        for phase=1:2 %prep and production
            for m=1:numel(MVA_ID) %loop for MVA type
                for s=1:numel(subj_name)
                    %                P{s}=fullfile(groupDir, ['szacc_' MVA_ID{m} '_160_Mov_' subj_name{s}, '.nii']);
                    P{s}=fullfile(groupDir, ['szacc_' MVA_ID{m} '_160_' Phase_ID{phase} '_' subj_name{s}, '.nii']);
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
        
        
end;