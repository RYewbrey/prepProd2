function varargout=prepProd2_surfImana(what,varargin)


% addpath(genpath('/Users/psu9d3/Documents/prepProd2/matlab'))
% addpath('/Users/psu9d3/Documents/MATLAB/toolboxes/spm12')
% addpath(genpath('/Users/psu9d3/Documents/MATLAB/toolboxes/freesurfer'))
% addpath(genpath('/Users/psu9d3/Documents/MATLAB/toolboxes/userfun'))
% addpath(genpath('/Users/psu9d3/Documents/MATLAB/toolboxes/suit'))
% addpath(genpath('/Users/psu9d3/Documents/MATLAB/toolboxes/caret'))
% addpath(genpath('/Users/psu9d3/Documents/MATLAB/toolboxes/tools'))
% addpath(genpath('/Users/psu9d3/Documents/MATLAB/toolboxes/region'))
% addpath(genpath('/Users/psu9d3/Documents/MATLAB/toolboxes/fmristat'))
% addpath(genpath('/Users/psu9d3/Documents/MATLAB/toolboxes/raincloud'))

% addpath(genpath('/Users/psu9d3/spm8-maint'))

baseDir= '/Users/psu9d3/Documents/prepProd2/data';
anatDir=[baseDir filesep 'imaging' filesep 'anatomicals'];
epiDir=[baseDir filesep 'imaging' filesep 'imaging_data'];
glmDir=[baseDir filesep 'imaging' filesep 'GLM_firstlevel'];
behDir=[baseDir filesep 'behavioural_data'];
dataext='run';
groupDir=[baseDir filesep 'imaging' filesep 'group_GLM_MNI'];
suitDir=[baseDir filesep 'imaging' filesep 'group_GLM_SUIT'];
BGDir=[baseDir filesep 'imaging' filesep 'group_GLM_BG'];
freesurferDir=[baseDir filesep 'imaging' filesep 'surfaceFreesurfer'];
caretDir=[baseDir filesep 'imaging' filesep 'surfaceCaret'];
regDir=[baseDir filesep 'imaging' filesep 'RegionOfInterest'];
simDir=[baseDir filesep 'imaging' filesep 'simuDim'];
regname={'M1','PM','V1'};

subj = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 ...
    37 38 39 40 41 42];
subj_name={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11','s12','s13','s14','s15','s16',...
    's17','s18','19','s20','s21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34',...
    's35','s36','s37','s38','s39','s40','s41','s42'}; %isotropic new data
anaSubj=[3 5 6 7 9 10 13 16 17 18 20 21 22 25 26 31 32 34 36 38 39 40 41 42];

nSubj = length(anaSubj);

atlasA={'i','x'};
atlas= 2;
atlasname={'fsaverage','fsaverage_sym'};

%Parameters
run={'1','2','3','4','5','6'};
%run={'2','3','4','5','6'};
nrVolumes=230;

radius=16;
numVox=160;

hem={'lh','rh'};
hemName={'LeftHem','RightHem'};
dim_name={'1','2','3','4','5','6','7','8','9'};


switch(what)
    %=========================================================================
    % Map Icosahedron
    %=========================================================================
    %     case 'surf_map_ico'
    %         sn=varargin{1};
    %         for i=sn
    %             freesurfer_mapicosahedron(subj_name{i},freesurferDir,'smoothing',1);
    %         end;
    
    
    case 'surf_xhemireg'             % STEP S1: cross-Register surfaces left / right hem
        sn=varargin{1};
        for i=sn
            freesurfer_registerXhem({subj_name{i}},freesurferDir,'hemisphere',[1:2]);
        end
        
    case 'surf_map_ico'              % STEP S2: Align to the new atlas surface
        
        sn=varargin{1};
        atlas=varargin{2};
        if (atlas==2)
            for i=sn
                freesurfer_mapicosahedron_xhem(subj_name{i},freesurferDir,'smoothing',1,'hemisphere',[1:2]);
            end
        else
            for i=sn
                freesurfer_mapicosahedron(subj_name{i},freesurferDir,'smoothing',1);
            end
        end
        
        %=========================================================================
        % Transfer data into Caret-format
        %=========================================================================
        
    case 'surf_make_caret'
        %         sn=varargin{1};
        %         for i=sn
        %             caret_importfreesurfer(['i' subj_name{i}],freesurferDir,caretDir);
        %         end;
        
        sn=varargin{1};
        atlas=varargin{2};
        for i=sn
            caret_importfreesurfer([atlasA{atlas} subj_name{i}],freesurferDir,caretDir);
        end;
        %=========================================================================
        % Map functional volumes on surface over caret
        %=========================================================================
        
    case 'surf_define_search' % Define Search light on the surface representation
        sn=varargin{1};
        hem=varargin{2};
        for s=sn
            for h=hem
                if h== 1
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
                    coord_pial= caret_load(fullfile(caret_subjDIR, 'lh.PIAL.coord'));
                    coord_white= caret_load(fullfile(caret_subjDIR, 'lh.WHITE.coord'));
                    topo= caret_load(fullfile(caret_subjDIR, 'lh.CLOSED.topo'));
                elseif h==2
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'RightHem');
                    coord_pial= caret_load(fullfile(caret_subjDIR, 'rh.PIAL.coord'));
                    coord_white= caret_load(fullfile(caret_subjDIR, 'rh.WHITE.coord'));
                    topo= caret_load(fullfile(caret_subjDIR, 'rh.CLOSED.topo'));
                else
                    disp('wrong hemisphere! Should be lh or rh ');
                end
                %load data get the roi for the lda
                refV= spm_vol(fullfile(glmDir, subj_name{s},'mask.nii'));
                epiInfo.mat= refV.mat;
                epiInfo.dim= refV.dim;
                epiInfo.mask=spm_read_vols(refV);
                node_range=1:size(coord_white.data,1);
                [LI,voxmin,voxmax,vORr]= surfing_voxelselection(coord_white.data',coord_pial.data',topo.data', [0 160],epiInfo, node_range ,[5,0,1]);
                LI=LI';
                save([caret_subjDIR '/surface_roi_160vox.mat'], 'LI','voxmin','voxmax','vORr');
                M=caret_struct('metric','data',vORr);
                caret_save(fullfile(caret_subjDIR,'radius160.metric'),M);
            end;
        end;
        
        %%%% Volume based searchlight on surface
        
    case 'MVA_search_surf'        % Volume-based search light for surface
        sn=varargin{1};
        atlas=varargin{2};
        refDir= glmDir;
        
        for s=sn
            for h=1:2
                caret_subjDIR = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                coord_pial= caret_load(fullfile(caret_subjDIR, [hem{h} '.PIAL.coord']));
                coord_white= caret_load(fullfile(caret_subjDIR, [hem{h} '.WHITE.coord']));
                topo= caret_load(fullfile(caret_subjDIR, [hem{h} '.CLOSED.topo']));
                surf(h).c1=coord_white.data';
                surf(h).c2=coord_pial.data';
                surf(h).f=topo.data';
            end;
            clear coord_pial coord_white coord_caret;
            volDef= spm_vol(fullfile(refDir, subj_name{s},'mask.nii'));
            volDef.mask=spm_read_vols(volDef);
            [LI,voxmin,voxmax,voxel,node,surfindx,depth,radvox]= lmva_voxelselection_surf(surf, [6 160],volDef);
            LI=LI';
            save(fullfile(refDir,subj_name{s}, 'vol_roi_160vox.mat'), 'LI','voxmin','voxmax','voxel','radvox');
            save(fullfile(refDir,subj_name{s}, 'vol_surf.mat'), 'voxel','node','surfindx','depth');
            
            V=volDef;
            X=zeros(V.dim);
            depth=depth+1;
            depth(isnan(depth))=1;
            X(voxel)=depth;
            V=volDef;
            V.dt=[4 0];
            V.pinfo=[3/100 0 0]';
            V.fname=fullfile(refDir,subj_name{s},'vol_roi_depth.nii');
            spm_write_vol(V,X);
        end;
        
        %         V.fname=fullfile(refDir,subj_name{s},'vol_roi_test.nii');
        %         X=zeros(V.dim);
        %         for i=1:length(LI);
        %             indx(i,1)=~isempty(LI{i});
        %         end;
        %         indx=find(indx);
        %
        %         X(LI{indx(1)})=1;
        %         spm_write_vol(V,X);
        %         keyboard;
        %
        
    case 'MVA_do_overall_Mov_surf'                 % Conduct the classification analysis 4 sequences
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
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
        
        for h=1:2 %hemisphere
            if h== 1
                caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
                metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Comb_160_Mov.metric',subj_name{s}));
                out={[metric ',accMov_L']};
            elseif h==2
                caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'RightHem');
                metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Comb_160_Mov.metric',subj_name{s}));
                out={[metric ',accMov_R']};
            end
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            tstart = tic
            lmva_spm(fullfile(caret_subjDIR,'surface_roi_160vox.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart)
            
        end
        
    case 'MVA_do_overall_Prep_surf'                 % Conduct the classification analysis 4 sequences
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
        end
        
        for h=1:2 %hemisphere
            if h== 1
                caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
                metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Comb_160_Prep.metric',subj_name{s}));
                out={[metric ',accPrep_L']};
            elseif h==2
                caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'RightHem');
                metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Comb_160_Prep.metric',subj_name{s}));
                out={[metric ',accPrep_R']};
            end
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end
            
            tstart = tic
            lmva_spm(fullfile(caret_subjDIR,'surface_roi_160vox.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
            telapsed = toc(tstart)
            
        end
        
    case 'MVA_do_spatOneout_Mov_surf'
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
            
            for h=1:2 %hemisphere
                if h== 1
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
                    metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Spat_160_Mov.metric',subj_name{s}));
                    out={[metric ',accSpatMov_L']};
                elseif h==2
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'RightHem');
                    metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Spat_160_Mov.metric',subj_name{s}));
                    out={[metric ',accSpatMov_R']};
                end
                
                [row,col] = find(prod>0);
                P={SPM.Vbeta(1:126).fname}';
                Pselect=[];
                for i=1:length(col)
                    Pselect{i,1}= P{col(i)};
                end
                
                tstart = tic
                lmva_spm(fullfile(caret_subjDIR,'surface_roi_160vox.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
                telapsed = toc(tstart);
            end
        end
        
    case 'MVA_do_spatOneout_Prep_surf'
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
            
            for h=1:2 %hemisphere
                if h== 1
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
                    metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Spat_160_Prep.metric',subj_name{s}));
                    out={[metric ',accSpatPrep_L']};
                elseif h==2
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'RightHem');
                    metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Spat_160_Prep.metric',subj_name{s}));
                    out={[metric ',accSpatPrep_R']};
                end
                
                
                [row,col] = find(prep>0);
                P={SPM.Vbeta(1:126).fname}';
                Pselect=[];
                for i=1:length(col)
                    Pselect{i,1}= P{col(i)};
                end
                
                tstart = tic
                lmva_spm(fullfile(caret_subjDIR,'surface_roi_160vox.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
                telapsed = toc(tstart);
            end
        end
        
    case 'MVA_do_tempOneout_Mov_surf'
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
            
            for h=1:2 %hemisphere
                if h== 1
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
                    metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Temp_160_Mov.metric',subj_name{s}));
                    out={[metric ',accTempMov_L']};
                elseif h==2
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'RightHem');
                    metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Temp_160_Mov.metric',subj_name{s}));
                    out={[metric ',accTempMov_R']};
                end
                
                [row,col] = find(prod>0);
                P={SPM.Vbeta(1:126).fname}';
                Pselect=[];
                for i=1:length(col)
                    Pselect{i,1}= P{col(i)};
                end
                
                
                
                tstart = tic
                lmva_spm(fullfile(caret_subjDIR,'surface_roi_160vox.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
                telapsed = toc(tstart);
            end
        end
        
    case 'MVA_do_tempOneout_Prep_surf'
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
            
            for h=1:2 %hemisphere
                if h== 1
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
                    metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Temp_160_Prep.metric',subj_name{s}));
                    out={[metric ',accTempPrep_L']};
                elseif h==2
                    caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'RightHem');
                    metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Temp_160_Prep.metric',subj_name{s}));
                    out={[metric ',accTempPrep_R']};
                end;
                
                [row,col] = find(prep>0);
                P={SPM.Vbeta(1:126).fname}';
                Pselect=[];
                for i=1:length(col)
                    Pselect{i,1}= P{col(i)};
                end;
                
                tstart = tic
                lmva_spm(fullfile(caret_subjDIR,'surface_roi_160vox.mat'),Pselect,out,@combinedclass,'params',{c,run,train,test});
                telapsed = toc(tstart);
            end
        end
        
    case 'MVA_do_Int_Mov_surf'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
        s=varargin{1};
        
        cd(fullfile(glmDir, subj_name{s}));
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
        
        for h=1:2 %hemisphere
            if h== 1
                caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
                metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Int_160_Mov.metric',subj_name{s}));
                out={[metric ',accIntMov_L']};
            elseif h==2
                caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'RightHem');
                metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Int_160_Mov.metric',subj_name{s}));
                out={[metric ',accIntMov_R']};
            end;
            
            [row,col] = find(prod>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(caret_subjDIR,'surface_roi_160vox.mat'),Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
            telapsed = toc(tstart)
            
        end
        
    case 'MVA_do_Int_Prep_surf'    %'integrated' (subtracts out main effects, i.e. common patterns for T1, T2, S1, S2 and classifies residual)
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
        
        for h=1:2 %hemisphere
            if h== 1
                caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'LeftHem');
                metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Int_160_Prep.metric',subj_name{s}));
                out={[metric ',accIntPrep_L']};
            elseif h==2
                caret_subjDIR = fullfile(caretDir,['x',subj_name{s}],'RightHem');
                metric=fullfile(caret_subjDIR,sprintf('%s_accuracy_Int_160_Prep.metric',subj_name{s}));
                out={[metric ',accIntPrep_R']};
            end;
            
            [row,col] = find(prep>0);
            P={SPM.Vbeta(1:126).fname}';
            Pselect=[];
            for i=1:length(col)
                Pselect{i,1}= P{col(i)};
            end;
            
            tstart = tic
            lmva_spm(fullfile(caret_subjDIR,'surface_roi_160vox.mat'),Pselect,out,@prepProd2_combinedclass_corrected4Main,'params',{c,run,train,test});
            telapsed = toc(tstart)
            
        end
        
    case 'glm_set_surf' %
        sn=varargin{1};
        prefix='ua';
        delay=0;
        dur=nrVolumes; %TRs per block
        
        for s=sn
            T=[]; %%% will save info for design file (SPM_info)
            
            %%% Correct for the number of dummy scans (not in nii file)
            num_dummys=3;
            
            %%%%%%%%%%%%%%%%%% for s01 only:
            %fname=fullfile(behDir,subj_name{s},['to1_' subj_name{s}
            %'_block.ana']);
            %%% for all others:
            %%%%%%%%%%%%%%%%%%
            
            %             %%% Load Beh. file with conditions and onsets
            fname=fullfile(behDir,subj_name{s},['data/to3_' subj_name{s} '.mat']);
            %             D=dload(fname);  % Load the behavioral file (if *.ana)
            load(fname); %if *mat file
            D=T;
            T=[];
            
            %%% Load Beh. file with additional beh. info, i.a. Errors
            fname=fullfile(behDir,subj_name{s},['data/to3_' subj_name{s} '_addRegressors.ana']);
            E=dload(fname);
            
            
            %%% General info for GLM
            J.dir = {fullfile(glmDir, subj_name{s})}; % Working directory of the GLM
            J.timing.units = 'scans';                   % Units (scans or sec for ons)
            J.timing.RT = 2.72;                     % TR
            J.timing.fmri_t = 16;                   % Mircotime resolution
            J.timing.fmri_t0 = 1;
            
            %%% Gather info for each run
            for r=1:numel(run)              % Loop over all the 6 runs
                %cond regressors:
                R=getrow(D,D.BN==r+subj_priorBN(sn)); %depends on subject: +0 in s1, +63 in s2 etc.
                %additional error regressors:
                aR=getrow(E,E.BN==r+subj_priorBN(sn)); %depends on subject: +0 in s1, +63 in s2 etc.
                
                %%% Load all images/TRs (N)
                for i=1:187                 % Loop over all images and put into the data structure
                    N{i} = [fullfile(baseDir, 'imaging_data',subj_name{s}, [prefix subj_name{s},'_run',run{r},'.nii,',num2str(i)])];
                end;
                J.sess(r).scans= N;         % Number of images
                
                %%% Assign Conditions
                for c=1:9                  % Loop over the conditions
                    idx=find(R.conditiontype==c);
                    
                    % Make new structure(S) that reflects what each of the
                    % regressors in the design matrix -> goes into T
                    S.RN(c,1)=r; %run n
                    S.SN(c,1)=s; %subject nr
                    S.condition(c,1)=c;
                    S.numtrials(c,1)=length(idx);
                    
                    % Info on condition name, onset and duration
                    J.sess(r).cond(c).name = sprintf('%d_%d',r,c);
                    J.sess(r).cond(c).onset = [R.startTR(idx)-num_dummys-delay]; %%Correct for # of dummies!
                    J.sess(r).cond(c).duration = dur;
                    J.sess(r).cond(c).tmod = 0;
                    J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                end;
                
                %%% Assign additional regressors:
                if sum(aR.error)==0 %%if no errors have been made in this block
                    c=c+1; %next higher regressor # for errors;
                    J.sess(r).cond(c).name = 'error';
                    J.sess(r).cond(c).onset = -99;
                    J.sess(r).cond(c).duration = 0; %%%error modelled as one TR
                    J.sess(r).cond(c).tmod = 0;
                    J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                else
                    c=c+1; %next higher regressor # for errors;
                    J.sess(r).cond(c).name = 'error';
                    J.sess(r).cond(c).onset = [aR.startTR(aR.error==1)-num_dummys-delay]; %%Correct for # of dummies!
                    J.sess(r).cond(c).duration = 0; %%%error modelled as one TR
                    J.sess(r).cond(c).tmod = 0;
                    J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                end;
                % Make new structure(S) that reflects what each of the
                % regressors in the design matrix -> goes into T
                idx=sum(aR.error); %%% INFO: missing responses are modelled as Errors
                S.RN(c,1)=r; %run n
                S.SN(c,1)=s; %subject nr
                S.condition(c,1)=c;
                S.numtrials(c,1)=idx; %Nr of trials in a run
                
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = 128;
                
                
                T=addstruct(T,S);
            end;
            J.fact = struct('name', {}, 'levels', {});
            
            %%%Derivatives:
            %J.bases.hrf.derivs = [0 0]; %%%used in first analysis
            J.bases.hrf.derivs = [1 0];
            
            J.volt = 1;
            J.global = 'None';
            %%J.mask = {[fullfile(anatDir,  subj_name{sn},
            %%[subj_name{sn},'_mask4glm.img'])],1};
            J.mask = {''};
            J.cvi =  'wls';
            matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec=J;
            spm_jobman('run',matlabbatch);
            dsave(fullfile(glmDir, subj_name{s},'SPM_info.ana'),T);  % save the information fil;e for the SPM
        end;
        
    case 'estimate_glm'
        % sl1_imana('estimate_glm')
        sn=varargin{1};
        for s=sn
            matlabbatch{1}.spm.tools.rwls.fmri_rwls_est.spmmat = { fullfile(glmDir, subj_name{s},'SPM.mat')};
            matlabbatch{1}.spm.tools.rwls.fmri_rwls_est.method.Classical = 1;
            spm_jobman('run',matlabbatch);
        end;
        
    case 'surf_map_acc'  % map the volumetric MVPA accuracy onto the surface, for comparison
        %sl4_imana('surf_map_con', 1:11, 1:2)
        % map volume images to metric file and save them in individual surface folder
        sn=varargin{1};
        hemisphere=[1:2];
        atlas=2;
        
        fileList={'_accuracy_Comb_160Uncorrected_surfvol.nii',...
            '_accuracy_Comb_160Corrected_surfvol.nii',...
            '_accuracy_Comb_160TempOneOut_surfvol.nii',...
            '_accuracy_Comb_160OrdOneOut_surfvol.nii'};
        
        %         fileList={'con_0001',... %for contrasts, use surf_map_con instead
        %             'con_0002',...
        %             'con_0003',...
        %             'con_0004',...
        %             'con_0005',...
        %             'con_0006'};
        
        for s=sn
            for h=hemisphere
                caretSDir = fullfile(caretDir,[atlasA{atlas},subj_name{s}],hemName{h});
                specname=fullfile(caretSDir,[atlasA{atlas},subj_name{s} '.' hem{h}   '.spec']);
                white=fullfile(caretSDir,[hem{h} '.WHITE.coord']);
                pial=fullfile(caretSDir,[hem{h} '.PIAL.coord']);
                
                C1=caret_load(white);
                C2=caret_load(pial);
                
                for f=1:length(fileList)
                    %                     images=fullfile(glmDir, subj_name{s},[subj_name{s}, fileList{f}]);
                    images=fullfile(glmDir, subj_name{s},[fileList{f} '.nii']);
                    M=caret_vol2surf_own(C1.data,C2.data,images,'ignore_zeros',1);
                    caret_save(fullfile(caretSDir,sprintf([subj_name{s}, fileList{f}, '.metric'])),M);
                    
                end;
                fprintf('%d\n',s);
            end;
        end;
        
    case 'surf_map_con' %map the volumetric contrasts onto the surface
%         fileList={'con_0001.nii','con_0002.nii','perc_0001.nii','perc_0002.nii',...
%             'con_0003.nii','con_0004.nii','con_0005.nii','con_0006.nii'};
%         outList={'con_0001','con_0002','perc_0001','perc_0002',...
%             'con_0003','con_0004','con_0005','con_0006'};
        
        fileList={'perc_0001.nii','perc_0002.nii'};
        outList={'perc_0001','perc_0002'};
        
        sn=varargin{1};
        he=1:2; %for both hemispheres
        %         he=varargin{2}; %for hemisphere = second input
        for s=sn
            for h=he
                caretSDir = fullfile(caretDir,['x',subj_name{s}],hemName{h});
                specname=fullfile(caretSDir,['x' subj_name{s} '.' hem{h}   '.spec']);
                for f=1:length(fileList)
                    
                    images{1}=fullfile(glmDir,subj_name{s},fileList{f});
                    images=caret_volprep(images);
                    
                    caret_vol2surf(specname,images,'outname',fullfile(caretSDir,[subj_name{s}, outList{f}, '.metric']));
                end;
                
            end;
        end;
        
    case 'surf_avrgcoord'       % Makes an average FIDICUAL (WHITE) surface
        
        atlas=2;
        j=0;
        avrgDir=[caretDir filesep atlasname{atlas}];
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
            
            if exist(surfaceGroupDir, 'dir') == 0
                mkdir(surfaceGroupDir)
            end
            
            cd(surfaceGroupDir);
            for i=anaSubj %subjname
                j=j+1;
                coord_name{j}=[caretDir filesep 'x' subj_name{i} filesep hemName{h} filesep hem{h} '.WHITE.coord'];
            end;
            topo_name= [hem{h} '.CLOSED.topo'];
            caret_avrgsurface('avrgDir', avrgDir,'topo_name', topo_name,'coord_name', coord_name);
        end
        
    case 'surf_makeGroup' % __________STEP S6:Make the group metric files Make group metric on func and accuracy data. Do contrast and zacc later on these files!
        
        atlas=2;
        %         anaSubj=[1:length(subj_name)]; %subjects used in the group analysis
        
        
        INname={'_accuracy_Comb_160_Mov',...
            '_accuracy_Comb_160_Prep',...
            '_accuracy_Int_160_Mov',...
            '_accuracy_Int_160_Prep',...
            '_accuracy_Spat_160_Mov',...
            '_accuracy_Spat_160_Prep'...
            '_accuracy_Temp_160_Mov',...
            '_accuracy_Temp_160_Prep',...
            'con_0001',...
            'con_0002',...
            'perc_0001',...
            'perc_0002',...
            'con_0003',...
            'con_0004',...
            'con_0005',...
            'con_0006'};
        
        
        OUTname={'_accuracy_Comb_160_Mov',...
            '_accuracy_Comb_160_Prep',...
            '_accuracy_Int_160_Mov'...
            '_accuracy_Int_160_Prep'...
            '_accuracy_Spat_160_Mov'...
            '_accuracy_Spat_160_Prep'...
            '_accuracy_Temp_160_Mov',...
            '_accuracy_Temp_160_Prep',...
            'con_0001'...
            'con_0002'...
            'perc_0001',...
            'perc_0002',...
            'con_0003'...
            'con_0004'...
            'con_0005'...
            'con_0006'};
        
        
        
        %         INname={'con_0001.nii'}; %for debugging, just one contrast
        %         OUTname={'con_0001.nii'};
        
        
        inputcol= ones(length(INname));
        replaceNaN= zeros(length(INname));
        
        vararginoptions(varargin,{'atlas'});
        
        for h=1:2
            
            surfaceGroupDir=[caretDir filesep atlasname{atlas} filesep hemName{h} ];
            cd(surfaceGroupDir);
            for j=1:length(INname) %----loop over each input metric file and make a group metric file
                k=0;
                for i=anaSubj %----define names of subj metric files
                    k=k+1;
                    infilenames{j}{k}=[caretDir filesep atlasA{atlas} subj_name{i} filesep hemName{h} filesep subj_name{i} INname{j} '.metric'];
                end
                
                %--`--name for group metric file in average surface folder
                outfilenames{j}=[surfaceGroupDir filesep hem{h} '.' OUTname{j} '.metric'];
                %----make the group metric
                fprintf('hem: %i  image: %i \n', h,j);
                
                caret_metricpermute(infilenames{j},'outfilenames',outfilenames(j),'inputcol',inputcol(j),'replaceNaNs',replaceNaN(j));
            end
        end
        
    case 'surf_zacc'            % ____STEP S7.1: Make the zacc files (on group level) and remove the NaN
        %sl4_imana('surf_zacc')
        
        numTests=6;
        numCat=4;
        mu=1/numCat; %mu=0.25;
        N=numTests*numCat;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        atlas =2;
        %         INname={'acc_Comb_160','acc_OrdOneout_160','acc_TempOneout_160'};
        %         OUTname={'zacc_Comb_160','zacc_OrdOneout_160','zacc_TempOneout_160'};
        
        %         INname={'acc_Comb_160'};
        %         OUTname={'zacc_Comb_160'};
        
        %         INname={'acc_Comb_160Uncorrected'};
        %         OUTname={'zacc_Comb_160Uncorrected'};
        
        INname={'_accuracy_Comb_160_Mov',...
            '_accuracy_Comb_160_Prep',...
            '_accuracy_Int_160_Mov'...
            '_accuracy_Int_160_Prep'};
        
        OUTname={'_zacc_Comb_160_Mov',...
            '_zacc_Comb_160_Prep',...
            '_zacc_Int_160_Mov'...
            '_zacc_Int_160_Prep'};
        
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas} filesep hemName{h} ];
            cd(surfaceGroupDir);
            for i=1:length(INname)
                iname=[hem{h} '.' INname{i} '.metric'];
                oname=[hem{h} '.' OUTname{i} '.metric'];
                fprintf('hem: %i  image: %i \n', h,i);
                caret_metriccalc({iname},oname,...
                    sprintf(' (X./X).*((X-%d)/%d)',mu,sigma),'replaceNaNs',1)
                %                 '(i1-1/4)/0.0884','replaceNaNs',1);
            end
        end
        
        
        %%%OneOut:
        
        takeOneOutIter=2;
        numTests=6;
        numCat=2;
        mu=1/numCat; %mu=0.5;
        N=numTests*numCat*takeOneOutIter;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        
        INname={'_accuracy_Spat_160_Mov'...
            '_accuracy_Spat_160_Prep'...
            '_accuracy_Temp_160_Mov',...
            '_accuracy_Temp_160_Prep'};
        OUTname={'_zacc_Spat_160_Mov'...
            '_zacc_Spat_160_Prep'...
            '_zacc_Temp_160_Mov',...
            '_zacc_Temp_160_Prep'};
        
        
        
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas} filesep hemName{h} ];
            cd(surfaceGroupDir);
            for i=1:length(INname)
                iname=[hem{h} '.' INname{i} '.metric'];
                oname=[hem{h} '.' OUTname{i} '.metric'];
                fprintf('hem: %i  image: %i \n', h,i);
                caret_metriccalc({iname},oname,...
                    sprintf(' (X./X).*((X-%d)/%d)',mu,sigma), 'replaceNaNs',1)
                %                '(i1-1/2)/0.1021','replaceNaNs',1);
            end;
        end;
        
    case 'surf_groupSmooth' %_________STEP S8: Smooth the group maps
        
        atlas=2;
        %fwhm kernel size for smoothing
%         kernelSize = 2;
%         kernelSize = 6; %same as elife, for comparison. Our default for this experiment.
        kernelSize = 4; %same as volumetric, but surface analysis requires broader smoothing?
%         kernelSize = 8; %elife was lower resolution, so perhaps this is more comparable?
%         kernelSize = 20; %very high kernel size, for troubleshooting
        
        vararginoptions(varargin,{'excludeSubj', 'Smooth_iterations','atlas'});
        
        
        SPMname={'_zacc_Comb_160_Mov'...
            '_zacc_Comb_160_Prep'...
            '_zacc_Int_160_Mov'...
            '_zacc_Int_160_Prep'...
            '_zacc_Spat_160_Mov'...
            '_zacc_Spat_160_Prep'...
            '_zacc_Temp_160_Mov'...
            '_zacc_Temp_160_Prep',...
            'con_0001',...
            'con_0002',...
            'perc_0001',...
            'perc_0002',...
            'con_0003',...
            'con_0004',...
            'con_0005',...
            'con_0006'};
        
        Smooth_iterations=ones(length(SPMname))*15; %smooth 15
        
        setenv('PATH', '/Applications/caret/macosx64_apps/caret_command.app/Contents/MacOS/') %sets terminal variables necessary to run caret_command
        
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h}];
            cd(surfaceGroupDir)
            
            %----define name of coord and topology
            coordfile=[caretDir filesep atlasname{atlas}  filesep hemName{h} filesep hem{h} '.WHITE.coord'];
            topofile=[caretDir filesep atlasname{atlas}  filesep hemName{h} filesep hem{h} '.CLOSED.topo'];
            
            %----get the full directory name of the metric files and the smoothed metric files that we create below
            for i=1:length(SPMname)
                filename=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '.metric']; % unsmoothed
                
             %----smooth the metric files and save them with the prefix 's'
                caret_smooth(filename,'coord',coordfile,'topo',topofile,'algorithm','FWHM','iterations',Smooth_iterations(i), 'fwhm', kernelSize); %FWHM smoothing
%                caret_smooth(filename,'coord',coordfile,'topo',topofile); %nearest neighbour smoothing
            end
        end
        
    case 'surf_groupTtest'
        
        col=1:nSubj; %nr of subjects arranged in columns
        
        atlas=2;
        
        SPMname={...
            '_zacc_Comb_160_Mov'...
            '_zacc_Comb_160_Prep'...
            '_zacc_Int_160_Mov'...
            '_zacc_Int_160_Prep'...
            '_zacc_Spat_160_Mov'...
            '_zacc_Spat_160_Prep'...
            '_zacc_Temp_160_Mov'...
            '_zacc_Temp_160_Prep',...
            'con_0001',...
            'con_0002',...
            'perc_0001',...
            'perc_0002',...
            'con_0003',...
            'con_0004',...
            'con_0005',...
            'con_0006'...
            };
        
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h}];
            cd(surfaceGroupDir)
            for i=1:length(SPMname)
                filename=[surfaceGroupDir filesep 's' hem{h} '.' SPMname{i} '.metric']; %
                %outfilename=[surfaceGroupDir filesep 'cSPM_s' hem{h} '.' SPMname{i} '.metric']; %
                
                %                     cSPM=caret_getcSPM('onesample_t','data', filename, col);
                cSPM=caret_getcSPM('onesample_t','data', filename, col,'maskthreshold', 0);
                
                %save as mat
                outfilename=[surfaceGroupDir filesep 'cSPM_s' hem{h} '.' SPMname{i} '.mat'];
                save(outfilename,'cSPM');
                
                %save as metric
                outfilename=[surfaceGroupDir filesep 'cSPM_s' hem{h} '.' SPMname{i} '.metric'];
                caret_savecSPM(outfilename,cSPM);
                
                % assemble summary structure
                data(:,i)=cSPM.con(1).con; % mean
                data(:,i+length(SPMname))=cSPM.con(1).Z; % T
                column_name{i}=['mean_' SPMname{i}];
                column_name{i+length(SPMname)}=['T_' SPMname{i}];
            end
            C=caret_struct('metric','data',data,'column_name',column_name);
            caret_save(fullfile(surfaceGroupDir,'summary.metric'),C);
        end
        
    case 'surf_groupTtest_unsmoothed'
        
        col=1:nSubj; %nr of subjects arranged in columns
        
        atlas=2;
        
        SPMname={'_zacc_Comb_160_Mov'...
            '_zacc_Comb_160_Prep'...
            '_zacc_Int_160_Mov'...
            '_zacc_Int_160_Prep'...
            '_zacc_Spat_160_Mov'...
            '_zacc_Spat_160_Prep'...
            '_zacc_Temp_160_Mov'...
            '_zacc_Temp_160_Prep',...
            'con_0001',...
            'con_0002',...
            'con_0003',...
            'con_0004',...
            'con_0005',...
            'con_0006'};
        
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h}];
            cd(surfaceGroupDir)
            for i=1:length(SPMname)
                filename=[surfaceGroupDir filesep hem{h} '.' SPMname{i} '.metric']; %
                %outfilename=[surfaceGroupDir filesep 'cSPM_s' hem{h} '.' SPMname{i} '.metric']; %
                
                cSPM=caret_getcSPM('onesample_t','data', filename, col);
                
                %save as mat
                outfilename=[surfaceGroupDir filesep 'cSPM_' hem{h} '.' SPMname{i} '.mat'];
                save(outfilename,'cSPM');
                
                %save as metric
                outfilename=[surfaceGroupDir filesep 'cSPM_' hem{h} '.' SPMname{i} '.metric'];
                caret_savecSPM(outfilename,cSPM);
                
                % assemble summary structure
                data(:,i)=cSPM.con(1).con; % mean
                data(:,i+length(SPMname))=cSPM.con(1).Z; % T
                column_name{i}=['mean_' SPMname{i}];
                column_name{i+length(SPMname)}=['T_' SPMname{i}];
            end;
            C=caret_struct('metric','data',data,'column_name',column_name);
            caret_save(fullfile(surfaceGroupDir,'summaryUnsmoothed.metric'),C);
        end;
        
    case 'percent_signal'
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            Vrunmean=SPM.Vbeta(SPM.xX.iB);
            P={Vrunmean.fname};
            spmj_imcalc_mtx(P,'meanEPI.nii','mean(X)');
            spm_imcalc({'meanEPI.nii','con_0001.nii'},'perc_0001.nii','i2./i1*100');
            spm_imcalc({'meanEPI.nii','con_0002.nii'},'perc_0002.nii','i2./i1*100');
            
        end
        
    case 'surf_ROI_MVA'             % Compute activated / classification area for SPL,M1,S1,SMA together
        
        nodenum=63416; %%L M1
        radius=160;
        white={'lh.FIDUCIAL.coord'};
        pial={'lh.FIDUCIAL.coord'};
        region('surf_circle',nodenum,radius,white,pial,topo,varargin);
        R=region_calcregions(R);
        
    case 'surf_ROI_stat_RY' %Functioning RoI script - uses 'elife_bilateral' RoI file, produces boxplots or bar charts
        
        plotType = varargin{1}; %1 for Bar chart, 2 for box plots, 3 for line graphs (paper figures)
        groupFigs = varargin{2}; %1 to group graphs in 1 figure, 2 to separate graphs into multiple figures
        fontSize = [20, 18]; %set font size for figures, first for title second for axes
        
        ROI_name='ROI.paint'; %load RoIs from RoI caret paint file
        ROI=caret_load(ROI_name);
        
        MVA_ID={'Comb','Temp','Spat','Int'};
        Phase_ID={'Prep','Mov'};
        
        ROIsave = []; %pre-assigning matrices to be fillled in loops
        lhROIData = [];
        rhROIData = [];
        
        for he=1:2 %for each hemisphere, extract data and then plot in bar or box plots
            
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{he} ];
            cd(surfaceGroupDir);
            
            k = 1; %loop counter for data extraction
            
            for i=1:numel(Phase_ID) %do prep & prod separately
                for j = 1:numel(MVA_ID) %look at each classifier within prep/prod
                    for l = unique(ROI.data') %and within each region, extract relevant cSPM.metric file
                        
                        P{k} = fullfile(surfaceGroupDir, ['cSPM_s' hem{he} '._zacc_' MVA_ID{j} '_160_' Phase_ID{i} '.metric']);
                        V = caret_load(P{k});
                        
                        mask = ROI.data == l; %take data from current ROI
                        
                        for sub=1:nSubj
                            data(sub) = mean(V.data(mask,sub));
                        end
                        
                        stats{l + 1} = data';
                        
                        k = k +1;
                    end
                    
                    S=[];
                    
                    for statSize=1:size(stats,2) %loop over ROIs
                        meanROI(:,statSize)=nanmean(stats{statSize},2); %average across voxels
                        D.meanROI=meanROI(:,statSize);
                        D.subj=anaSubj';
                        D.roiID=repmat(statSize,numel(anaSubj),1);
                        D.mvaID=ones(numel(anaSubj),1)*j;
                        D.phase=ones(numel(anaSubj),1)*i;
                        S=addstruct(S,D);%add ROI by ROI
                        D=[];
                    end
                    ROIsave=addstruct(ROIsave,S); %add MVA type by type
                    S=[];
                    fileOut=[hem{he} MVA_ID{j} '_' Phase_ID{i} '_ROI_spss.mat'];
                    save(fileOut,'meanROI');
                end
                
                if he == 1
                    lhROIData= addstruct(lhROIData,ROIsave); %add MVA type by type
                    ROIsave=[];
                else
                    rhROIData=addstruct(rhROIData,ROIsave); %add phase by phase
                    ROIsave=[];
                end
            end
            
            if he==1
                
                %%% Figure
                %black: combined; red: temporal; blue: spatial; green: integrated
                %                 color={[0 0 0],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
                
                %same colours not including combined/overall
                color={[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
                
                %color={[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
                %                 figure('Renderer', 'painters', 'Position', [10 10 400 1000])
                %                 ROIname={'M1','PMd','SMA','SPCa' 'SPCp'};
                ROIname={'L-M1','L-S1','L-PMd','L-PMv','L-SMA','L-SPCa' 'L-SPCp'};
                
                loopCounter = 1;
                
                for i=[3 2 4 5 6 8 9] %loop over ROIs
                    
                    if groupFigs ==1
                        %For all RoIs on one plot
                        subplot(7,1,loopCounter)
                        
                    elseif groupFigs ==2
                        %For RoIs on individual plots
                        figure
                        
                    elseif groupFigs==1 && plotType == 3
                        %For when line plots are shown as one figure
                        subplot(1,7,loopCounter)
                        
                    end
                    
                    
                    %data for plots not including overall/combined classifier
                    index = lhROIData.mvaID > 1;
                    lhROIDataNoComb.meanROI = lhROIData.meanROI(index);
                    lhROIDataNoComb.roiID = lhROIData.roiID(index);
                    lhROIDataNoComb.mvaID = lhROIData.mvaID(index);
                    lhROIDataNoComb.phase = lhROIData.phase(index);
                    %
                    
                    
                    %data for plots with only overall/combined classifier
                    %                     combColor = [0 0 0];
                    index = lhROIData.mvaID == 1;
                    lhROIDataComb.meanROI = lhROIData.meanROI(index);
                    lhROIDataComb.roiID = lhROIData.roiID(index);
                    lhROIDataComb.mvaID = lhROIData.mvaID(index);
                    lhROIDataComb.phase = lhROIData.phase(index);
                    
                    
                    
                    if plotType == 1 %uncomment which type of plot you want
                        
                        %barplot containing all classifiers
%                         barplot([lhROIData.phase lhROIData.mvaID],lhROIData.meanROI,'split',lhROIData.mvaID,'subset',lhROIData.roiID==i&lhROIData.mvaID>0,'facecolor',color,...
%                         'leglocation','northeast') ;
                        
                        %barplot not containing overall/combined classifier
%                         barplot([lhROIDataNoComb.phase lhROIDataNoComb.mvaID],lhROIDataNoComb.meanROI,'split',lhROIDataNoComb.mvaID,'subset',lhROIDataNoComb.roiID==i&lhROIDataNoComb.mvaID>0,'facecolor',color,...
%                         'leglocation','northeast') ;
                        
                        %barplot with only overall/combined classifier
%                         barplot([lhROIDataComb.phase lhROIDataComb.mvaID],lhROIDataComb.meanROI,'split',lhROIDataComb.mvaID,'subset',lhROIDataComb.roiID==i&lhROIDataComb.mvaID>0,'facecolor',combColor,...
%                         'leglocation','northeast') ;
                        
                    elseif plotType == 2 %uncomment which type of plot you want
                        
                        %boxplot containing all classifiers
%                         myboxplot([lhROIData.phase lhROIData.mvaID],lhROIData.meanROI,'split',lhROIData.mvaID,'subset',lhROIData.roiID==i&lhROIData.mvaID>0,'fillcolor',color,...
%                         'leglocation','northeast') ;
                        
                        %boxplot not containing overall/combined classifier
                        myboxplot([lhROIDataNoComb.phase lhROIDataNoComb.mvaID],lhROIDataNoComb.meanROI,'split',lhROIDataNoComb.mvaID,'subset',lhROIDataNoComb.roiID==i&lhROIDataNoComb.mvaID>0,'fillcolor',color,...
                            'leglocation','northeast') ;
                        
                    elseif plotType == 3 %for paper
                        
                        lineplot(lhROIDataNoComb.phase, lhROIDataNoComb.meanROI, 'split', lhROIDataNoComb.mvaID,...
                            'subset',lhROIDataNoComb.roiID==i&lhROIDataNoComb.mvaID>0,'style_thickline', 'markersize', 5, ...
                            'linecolor', color, 'errorcolor', color, 'markerfill', color, 'markercolor',color) ;
                        
                        
                    end
                    
                    
                    title(ROIname{loopCounter}, 'FontSize', fontSize(1), 'FontName', 'Calibri')
                    drawline(0,'dir','horz'); %draw line at chance level
                    
%                     ylim([-2.5 3]); %boxplot axes for writeup
%                     yticks([-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3])
                    
%                     ylim([0 1.5]) %elife comparison bar charts
                    
                    ylim([-0.5 1]) %good ylim for line plot

                    set(gca,'FontSize',fontSize(2), 'FontName', 'Calibri')

                    loopCounter = loopCounter +1;
                end
                
                save('PrepProd_ROI.mat','-struct','lhROIData');
                figure
            elseif he==2
                
                                %%% Figure
                %black: combined; red: temporal; blue: spatial; green: integrated
                %                 color={[0 0 0],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
                
                %same colours not including combined/overall
                color={[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
                
                %color={[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.4660 0.6740 0.1880]};
                %                 figure('Renderer', 'painters', 'Position', [10 10 400 1000])
                %                 ROIname={'M1','PMd','SMA','SPCa' 'SPCp'};
                ROIname={'R-M1','R-S1','R-PMd','R-PMv','R-SMA','R-SPCa' 'R-SPCp'};
                
                loopCounter = 1;
                
                for i=[3 2 4 5 6 8 9] %loop over ROIs
                    
                    if groupFigs ==1
                        %For all RoIs on one plot
                        subplot(7,1,loopCounter)
                        
                    elseif groupFigs ==2
                        %For RoIs on individual plots
                        figure
                    end
                    
                    
                    %data for plots not including overall/combined classifier
                    index = rhROIData.mvaID > 1;
                    rhROIDataNoComb.meanROI = rhROIData.meanROI(index);
                    rhROIDataNoComb.roiID = rhROIData.roiID(index);
                    rhROIDataNoComb.mvaID = rhROIData.mvaID(index);
                    rhROIDataNoComb.phase = rhROIData.phase(index);
                    %
                    
                    
                    %data for plots with only overall/combined classifier
                    %                     combColor = [0 0 0];
                    index = rhROIData.mvaID == 1;
                    rhROIDataComb.meanROI = rhROIData.meanROI(index);
                    rhROIDataComb.roiID = rhROIData.roiID(index);
                    rhROIDataComb.mvaID = rhROIData.mvaID(index);
                    rhROIDataComb.phase = rhROIData.phase(index);
                    
                    
                    
                    if plotType == 1 %uncomment which type of plot you want
                        
                        %barplot containing all classifiers
%                         barplot([rhROIData.phase rhROIData.mvaID],rhROIData.meanROI,'split',rhROIData.mvaID,'subset',rhROIData.roiID==i&rhROIData.mvaID>0,'facecolor',color,...
%                         'leglocation','northeast') ;
                        
                        %barplot not containing overall/combined classifier
                        barplot([rhROIDataNoComb.phase rhROIDataNoComb.mvaID],rhROIDataNoComb.meanROI,'split',rhROIDataNoComb.mvaID,'subset',rhROIDataNoComb.roiID==i&rhROIDataNoComb.mvaID>0,'facecolor',color,...
                        'leglocation','northeast') ;
                        
                        %barplot with only overall/combined classifier
%                         barplot([rhROIDataComb.phase rhROIDataComb.mvaID],rhROIDataComb.meanROI,'split',rhROIDataComb.mvaID,'subset',rhROIDataComb.roiID==i&rhROIDataComb.mvaID>0,'facecolor',combColor,...
%                         'leglocation','northeast') ;
                        
                    elseif plotType == 2 %uncomment which type of plot you want
                        
                        %boxplot containing all classifiers
%                         myboxplot([rhROIData.phase rhROIData.mvaID],rhROIData.meanROI,'split',rhROIData.mvaID,'subset',rhROIData.roiID==i&rhROIData.mvaID>0,'fillcolor',color,...
%                         'leglocation','northeast') ;
                        
                        %boxplot not containing overall/combined classifier
                        myboxplot([rhROIDataNoComb.phase rhROIDataNoComb.mvaID],rhROIDataNoComb.meanROI,'split',rhROIDataNoComb.mvaID,'subset',rhROIDataNoComb.roiID==i&rhROIDataNoComb.mvaID>0,'fillcolor',color...
                        ,'xtickoff') ;
                        
                    elseif plotType == 3 %for paper
                        
                        lineplot(rhROIDataNoComb.phase, rhROIDataNoComb.meanROI, 'split', rhROIDataNoComb.mvaID,...
                            'subset',rhROIDataNoComb.roiID==i&rhROIDataNoComb.mvaID>0,'style_thickline', 'markersize', 5, ...
                            'linecolor', color, 'errorcolor', color, 'markerfill', color, 'markercolor',color) ;                        
                                                
                    end
                    
                    
                    title(ROIname{loopCounter}, 'FontSize', fontSize(1), 'FontName', 'Calibri')
                    drawline(0,'dir','horz'); %draw line at chance level
                    
%                     ylim([-2.5 3]); %boxplot axes for writeup
%                     yticks([-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3])
                    
%                     ylim([0 1.5]); %elife comparison bar charts
                    
                    ylim([-0.5 1]) %good ylim for line plot
                    
                    set(gca,'FontSize',fontSize(2), 'FontName', 'Calibri')

                    loopCounter = loopCounter +1;
                end
                
                save('PrepProd_ROI.mat','-struct','lhROIData');
                
            end
        end
        
    case 'result_crosssectionContrasts_RY' %produces cross-sections of con images and percent signal change
        
        set(0, 'DefaultAxesFontSize', 12, 'DefaultAxesFontName', 'Calibri')
        
        subj=(1:length(anaSubj)); %because we extract data from saved files
        crossSectionType = {'profile.borderproj', 'premotor.borderproj'};
        crossSectionName = {'Profile', 'Premotor'};
        
        hemiName={'lh','rh'};
        
        for h=1:2 %loop through hemispheres, plot for each
            groupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            
            %border files to extract data from - profile (PM to SPC) or premotor (PMv to SMA). Uncomment accordingly
%             border=fullfile(groupDir,'profile.borderproj');
%             border=fullfile(groupDir,'premotor.borderproj');
            border=fullfile(groupDir,crossSectionType(varargin{1}));
            border = char(border);
            
            %             X=caret_crosssection(border,fullfile(groupDir,[hem{h} '.INFLATED.coord']));
            S=caret_crosssection(border,fullfile(groupDir,[hem{h} '.surface_shape']));
            
            T.conMov=caret_crosssection(border,fullfile(groupDir,['s' hemiName{h} '.con_0001.metric']))';
            T.conPrep=caret_crosssection(border,fullfile(groupDir,['s' hemiName{h}  '.con_0002.metric']))'; %
            
            
            T.conMov=T.conMov(subj,:); %caret_crosssection output includes betas, subjN etc. here we just take the data.
            T.conPrep=T.conPrep(subj,:);
            
            T.subj=(1:length(T.conMov))'; %add subj to our struct so we can split by it
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots
            
            %plot preparation data into traceplots
            x=1:size(T.conPrep,2);
            figure;
            subplot(4,1,[1:3]);
            traceplot(x,T.conPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 1 0],'errorfcn','stderr');

            title('Preparation Versus Rest')
            ylim([-1.5, 5])
            drawline(0,'dir','horz');
            
            subplot(4,1,4);
            plot(x,S(:,2));
            ylim([-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            
            %plot movement data into traceplots
            x=1:size(T.conMov,2);
            figure;
            subplot(4,1,[1:3]);
            traceplot(x,T.conMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 1 0],'errorfcn','stderr');

            title('Movement Versus Rest')
            
            hold off;
            ylim([-1.5, 5])
            drawline(0,'dir','horz');
            
            subplot(4,1,4);
            plot(x,S(:,2));
            ylim([-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            
            %----------------------- Percent Signal Change --------------------
            
            T.percMov=caret_crosssection(border,fullfile(groupDir,['s' hemiName{h} '.perc_0001.metric']))';
            T.percPrep=caret_crosssection(border,fullfile(groupDir,['s' hemiName{h}  '.perc_0002.metric']))'; %
            
            
            T.percMov=T.percMov(subj,:); %caret_crosssection output includes betas, subjN etc. here we just take the data.
            T.percPrep=T.percPrep(subj,:);
            
            T.subj=(1:length(T.percMov))'; %add subj to our struct so we can split by it
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots
            
            %plot preparation data into traceplots
            x=1:size(T.percPrep,2);
            figure;
            subplot(4,1,[1:3]);
            traceplot(x,T.percPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 0 0],'errorfcn','stderr','patchcolor',[0.1 0.1 0.1]);
            
            title('Preparation Versus Rest (% Signal Change)')
            ylim([-0.2, 2])
            yticks(0:0.2:2)
            ylabel('% Signal Change')
            drawline(0,'dir','horz');
            
            subplot(4,1,4);
            plot(x,S(:,2));
            ylim([-1.6 1.5]);
            
            set(gcf,'Position',[10 10 700 1000]);
            
            %plot movement data into traceplots
            x=1:size(T.percMov,2);
            figure;
            subplot(4,1,[1:3]);
            traceplot(x,T.percMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 0 0],'errorfcn','stderr','patchcolor',[0.1 0.1 0.1]);
            
            title('Movement Versus Rest (% Signal Change)')
            
            hold off;
            ylim([-0.2, 2])
            yticks(0:0.2:2)
            ylabel('% Signal Change')
            drawline(0,'dir','horz');
            
            subplot(4,1,4);
            plot(x,S(:,2));
            ylim([-1.6 1.5]);
            ylabel('Sulcul Depth (mm)')
            set(gcf,'Position',[10 10 700 1000]);
            
            filename = ['/Users/psu9d3/Documents/prepProd2/docs/crossSections/', 'PrepProd_crossSectionPerc_', hemiName(h), crossSectionName(varargin{1}), '.mat'];
            filename = strcat(filename{1}, filename{2}, filename{3}, filename{4}, filename{5});
            save(filename,'-struct','T')
            
            
        end
        
    case 'result_crosssection_RY'  % Generates Cross-section
        
        subj=(1:length(anaSubj)); %because we extract data from saved files
        
        crossSectionType = {'profile.borderproj', 'premotor.borderproj'};
        crossSectionName = {'Profile', 'Premotor'};
        hemiName={'lh','rh'};
        D=[];
        A=[];
        for h=1:2
            groupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
%             border=fullfile(groupDir,'profile.borderproj');
%             border=fullfile(groupDir,'premotor.borderproj');
            border=fullfile(groupDir,crossSectionType(varargin{1}));
            border = char(border);
%             X=caret_crosssection(border,fullfile(groupDir,[hem{h} '.INFLATED.coord']));
            S=caret_crosssection(border,fullfile(groupDir,[hem{h} '.surface_shape']));
            
            T.combMov=caret_crosssection(border,fullfile(groupDir,['cSPM_s' hemiName{h}  '._zacc_Comb_160_Mov.metric']))'; %
            T.combPrep=caret_crosssection(border,fullfile(groupDir,['cSPM_s' hemiName{h}  '._zacc_Comb_160_Prep.metric']))'; %
            
            T.intMov=caret_crosssection(border,fullfile(groupDir,['cSPM_s' hemiName{h}  '._zacc_Int_160_Mov.metric']))'; %
            T.intPrep=caret_crosssection(border,fullfile(groupDir,['cSPM_s' hemiName{h}  '._zacc_Int_160_Prep.metric']))'; %
            
            T.tempMov=caret_crosssection(border,fullfile(groupDir,['cSPM_s' hemiName{h}  '._zacc_Temp_160_Mov.metric']))'; %
            T.tempPrep=caret_crosssection(border,fullfile(groupDir,['cSPM_s' hemiName{h}  '._zacc_Temp_160_Prep.metric']))'; %
            
            T.spatMov=caret_crosssection(border,fullfile(groupDir,['cSPM_s' hemiName{h}  '._zacc_Spat_160_Mov.metric']))'; %
            T.spatPrep=caret_crosssection(border,fullfile(groupDir,['cSPM_s' hemiName{h} '._zacc_Spat_160_Prep.metric']))'; %
            
            
            T.combMov=T.combMov(subj,:); %caret_crosssection output includes betas, subjN etc. here we just take the data.
            T.combPrep=T.combPrep(subj,:);
            T.intMov=T.intMov(subj,:);
            T.intPrep=T.intPrep(subj,:);
            T.tempMov=T.tempMov(subj,:);
            T.tempPrep=T.tempPrep(subj,:);
            T.spatMov=T.spatMov(subj,:);
            T.spatPrep=T.spatPrep(subj,:);
            
            %             T.hem=ones(size(T.combMov),1)*h;
            T.subj=(1:length(T.combMov))';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plots
            
            %plot preparation data into traceplots
            x=1:size(T.combPrep,2);
            figure;
            subplot(4,1,[1:3]);
            traceplot(x,T.intPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 1 0],'errorfcn','stderr','patchcolor',[0 1 0]);
            hold on;
            traceplot(x,T.tempPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[1 0 0],'errorfcn','stderr','patchcolor',[1 0 0]);
            traceplot(x,T.spatPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 0 1],'errorfcn','stderr','patchcolor',[0 0 1]);
            title('preparation')
            
            hold off;
            set(gca,'XLim',[1 length(x)],'YLim',[-1.1 1.6]);
            drawline(0,'dir','horz');
            %         drawline(49,'dir','vert');
            
            subplot(4,1,4);
            plot(x,S(:,2));
            set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            varargout={T};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            figure;
            subplot(4,1,1:3)
            traceplot(x,T.combPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 0 0],'errorfcn','stderr','patchcolor',[0.1 0.1 0.1]);
            title('combPrep')
            
            set(gca,'XLim',[1 length(x)],'YLim',[-1.1 2]);
            drawline(0,'dir','horz');
            %         drawline(49,'dir','vert');
            
            subplot(4,1,4);
            plot(x,S(:,2));
            set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            varargout={T};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %plot movement data into traceplots
            x=[1:size(T.combMov,2)];
            figure;
            subplot(4,1,[1:3]);
            traceplot(x,T.intMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 1 0],'errorfcn','stderr','patchcolor',[0 1 0]);
            hold on;
            traceplot(x,T.tempMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[1 0 0],'errorfcn','stderr','patchcolor',[1 0 0]);
            traceplot(x,T.spatMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 0 1],'errorfcn','stderr','patchcolor',[0 0 1]);
            title('production')
            
            hold off;
            set(gca,'XLim',[1 length(x)],'YLim',[-1.1 1.6]);
            drawline(0,'dir','horz');
            %         drawline(49,'dir','vert');
            
            subplot(4,1,4);
            plot(x,S(:,2));
            set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            varargout={T};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            figure;
            subplot(4,1,1:3)
            traceplot(x,T.combMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 0 0],'errorfcn','stderr','patchcolor',[0.1 0.1 0.1]);
            title('combMov')
            
            set(gca,'XLim',[1 length(x)],'YLim',[-1.1 2]);
            drawline(0,'dir','horz');
            %         drawline(49,'dir','vert');
            
            subplot(4,1,4);
            plot(x,S(:,2));
            set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            varargout={T};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % All crosssections plotted together
            
            figure
            
            %%%%%%%%%%%%%%%%%%%% Prep group plot
            
            %%%%% comb
            subplot(12,2,[1 2 3 4])
            traceplot(x,T.combPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 0 0],'errorfcn','stderr','patchcolor',[0.1 0.1 0.1]);
            title('Overall Decoding Preparation')
            
            set(gca,'XLim',[1 length(x)],'YLim',[-1.1 2]);
            drawline(0,'dir','horz');
            %         drawline(49,'dir','vert');
            
            subplot(12,2,[5 6]);
            plot(x,S(:,2));
            set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            
            %%%%% features
            x = 1:size(T.combPrep,2);
            subplot(12,2,[7 8 9 10])
            traceplot(x,T.intPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 1 0],'errorfcn','stderr','patchcolor', [0 1 0]);
            hold on;
            traceplot(x,T.tempPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[1 0 0],'errorfcn','stderr','patchcolor', [1 0 0]);
            traceplot(x,T.spatPrep,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 0 1],'errorfcn','stderr','patchcolor', [0 0 1]);
            title('Feature Decoding Preparation')
            
            hold off;
            set(gca,'XLim',[1 length(x)],'YLim',[-1.1 1.6]);
            drawline(0,'dir','horz');
            %         drawline(49,'dir','vert');
            
            subplot(12,2,[11 12]);
            plot(x,S(:,2));
            set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            
            %%%%%%%%%%%%%%%%%%%%% Prod group plot
            
            %%%%% overall
            subplot(12,2,[13 14 15 16])
            traceplot(x,T.combMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 1 0],'errorfcn','stderr', 'patchcolor', [0 1 0]);
            title('Overall Decoding Production')
            
            set(gca,'XLim',[1 length(x)],'YLim',[-1.1 2]);
            drawline(0,'dir','horz');
            %         drawline(49,'dir','vert');
            
            subplot(12,2,[17 18]);
            plot(x,S(:,2));
            set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            varargout={T};
            
            %%%%% features
            x=1:size(T.combMov,2);
            subplot(12,2,[19 20 21 22]);
            traceplot(x,T.intMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 1 0],'errorfcn','stderr', 'patchcolor', [0 1 0]);
            hold on;
            traceplot(x,T.tempMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[1 0 0],'errorfcn','stderr', 'patchcolor', [1 0 0]);
            traceplot(x,T.spatMov,'linecolor','k','linewidth',3,'linestyle','-','linecolor',[0 0 1],'errorfcn','stderr', 'patchcolor', [0 0 1]);
            title('Feature Decoding Production')
            
            hold off;
            set(gca,'XLim',[1 length(x)],'YLim',[-1.1 1.6]);
            drawline(0,'dir','horz');
            %         drawline(49,'dir','vert');
            
            subplot(12,2,[23 24]);
            plot(x,S(:,2));
            set(gca,'Box','off','XLim',[1 length(x)],'YLim',[-1.6 1.5]);
            set(gcf,'Position',[10 10 500 1000]);
            varargout={T};
            
            filename = ['/Users/psu9d3/Documents/prepProd2/docs/crossSections/', 'PrepProd_crossSection_', hemiName(h), crossSectionName(varargin{1}), '.mat'];
            filename = strcat(filename{1}, filename{2}, filename{3}, filename{4}, filename{5});
            save(filename,'-struct','T')
            
        end
            
        
    case 'SUIT_Segment_RY'
        
        cd(anatDir)
        
        for i = 1:length(subj_name)
            
            fileName = sprintf('%s/%s_anatomical.nii', subj_name{i}, subj_name{i});
            
            suit_isolate_seg({fileName})
            
            
        end
        
    case 'surf_ROI_stat'             % Compute activated / classification area for SPL,M1,S1,SMA together
        %sl4_imana('surf_area_stat')
        
        
        atlas=2;
        dim=0; %set dim analysis to 0 (switch on below)
        
        prefix={'right', 'right'};
        
        thres=[2.6 2.6 0 0 0];
        %           thres=[1 1 0 0 0];
        
        %            thres=[5 3 0 0 0];
        %         thres=[5 2.6 0 0 0];
        %          thres=[5 4 0 0 0];
        %          thres=[1 1 0 0 0];
        
        % infilename={'%s.con.metric',...
        %     '%s.zacc_Comb_160Uncorrected.metric',...
        %     '%s.zacc_Comb_160.metric',...
        %     '%s.zacc_OrdOneout_160.metric',...
        %     '%s.zacc_TempOneout_160.metric'};
        
        infilename={'_zacc_Comb_160_Mov'...
            '_zacc_Comb_160_Prep'...
            '_zacc_Int_160_Mov'...
            '_zacc_Int_160_Prep'...
            '_zacc_Spat_160_Mov'...
            '_zacc_Spat_160_Prep'...
            '_zacc_Temp_160_Mov'...
            '_zacc_Temp_160_Prep',...
            'con_0001',...
            'con_0002',...
            'con_0003',...
            'con_0004',...
            'con_0005',...
            'con_0006'};
        
        % varname={'ttest','zacc_combU','zacc_comb','zacc_ord','zacc_temp'};
        varname={'zacc_comb_mov','zacc_comb_prep','zacc_int_mov','zacc_int_prep','zacc_spat_mov','zacc_spat_prep',...
            'zacc_temp_mov','zacc_temp_prep','ttest','ttest','ttest','ttest','ttest','ttest',};
        
        
        
        
        % %         %%Mean Searchlight (uncorrected only)
        %%ROI Dimensionalisty analysis:
        %         dim=1; %if dimensionality analysis
        %
        %         infilename={'%s.con.metric',...
        %             '%s.zacc_Comb_160Uncorrected_Dim1.metric','%s.zacc_Comb_160Uncorrected_Dim2.metric',...
        %             '%s.zacc_Comb_160Uncorrected_Dim3.metric','%s.zacc_Comb_160Uncorrected_Dim4.metric',...
        %             '%s.zacc_Comb_160Uncorrected_Dim5.metric','%s.zacc_Comb_160Uncorrected_Dim6.metric',...
        %             '%s.zacc_Comb_160Uncorrected_Dim7.metric','%s.zacc_Comb_160Uncorrected_Dim8.metric'};
        
        %          infilename={'%s.con.metric',...
        %             '%s.acc_Comb_160Uncorrected_Dim1.metric','%s.acc_Comb_160Uncorrected_Dim2.metric',...
        %             '%s.acc_Comb_160Uncorrected_Dim3.metric','%s.acc_Comb_160Uncorrected_Dim4.metric',...
        %             '%s.acc_Comb_160Uncorrected_Dim5.metric','%s.acc_Comb_160Uncorrected_Dim6.metric',...
        %             '%s.acc_Comb_160Uncorrected_Dim7.metric','%s.acc_Comb_160Uncorrected_Dim8.metric'};
        
        %         varname={'ttest','dim1','dim2','dim3','dim4','dim5','dim6','dim7','dim8'};
        
        
        
        reg='all'; %which regions
        %         reg=3; %{'???'  'S1'  'M1'  'PMd'  'PMv'  'SMA'  'V12'  'SPLa'  'SPLp'}
        %                reg=29; %DESIKAN: which regions? 5: corpus callosum (baseline); 18: paracentral; 19: parsopercularis; 23: postcentral; 25: precentral; 26: precuneus; 29: superiorfrontal; 30: superiorparietal
        %%%31: superior temporal; 32: supramarginal; 35: transversetemporal (Heschl);
        %                  reg=3;
        funMask=99; %1; set to 99 if no functional mask
        
        D=[];
        A=[];
        vararginoptions(varargin,{'thres','infilename','atlas','ROI_name','reg'});
        
        for he=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas}  filesep hemName{he} ];
            cd(surfaceGroupDir);
            
            
            
            
            %%Anat ROI:
            %                              ROI_name=[hem{he} '.desikan.paint'];
            %                                  ROI_name=[hem{he} '.mask.paint'];
            %                         ROI_name=[hem{he} '.ROIs_BAs.paint'];
            %                                              ROI_name=[hem{he} '.M1handknob.paint'];
            %                                                      ROI_name=[hem{he} '.SMA.paint'];
            %                                                        ROI_name=[hem{he} '.PMd_new.paint'];
            
            ROI_name=['ROI.paint'];
            
            
            ROI=caret_load(ROI_name);
            
            if strcmp(reg,'all')
                RegName='all';
            else
                index=unique(ROI.data);
                RegName=ROI.paintnames(reg);
                %                 RegName=ROI.column_name(reg);
            end;
            
            if (strcmp(reg,'all'))
                Mask=ROI.data>0;
            else
                Mask=isincluded(index(reg),ROI.data);
            end;
            
            %%Functional ROI:
            
            %(M>R):
            %             ROIfun_name=[hem{he} '.con.metric'];
            %             ROIfun=caret_load(ROIfun_name);
            %             MaskFun=ROIfun.data>thres(1);
            
            %zacc Comb uncorrected:
            %             ROIfun_name=[hem{he} '.zacc_Comb_160Uncorrected.metric'];
            %             ROIfun=caret_load(ROIfun_name);
            %             MaskFun=ROIfun.data>thres(2);
            
            
            %             %cSPM contrast
            %             ROIfun_name=['cSPM_s' hem{he} '.con.metric'];
            %             ROIfun=caret_load(ROIfun_name);
            %             ROIfun.data=ROIfun.data(:,end);%if cSPM
            %             MaskFun=ROIfun.data>thres(1);
            %             MaskFun=repmat(MaskFun,1,16); %if cSPM
            
            
            %           % cSPM zacc
            ROIfun_name=['cSPM_s' hem{he} '._zacc_Comb_160_Movright.metric'];
            
            %                            ROIfun_name=['cSPM_s' hem{he} '.zacc_OrdOneout_160.metric'];
            ROIfun=caret_load(ROIfun_name);
            ROIfun.data=ROIfun.data(:,end);%if cSPM
            MaskFun=ROIfun.data>thres(2);
            MaskFun=repmat(MaskFun,1,(length(subj_name))); %if cSPM
            %
            %
            %             %%write out subject and hemi:
            vec=ones((length(subj_name)),1); %9 subjects
            D.SN=subj_name';
            D.hem=vec*he;
            
            
            
            for i=1:length(infilename) % image for analysis
                M=caret_load(sprintf(fullfile(['cSPM_s' hem{he} '.' infilename{i} '%s' '.metric']),prefix{he}));
                
                
                %%Mask: Conjunction Anat and Func ROI (M>R) only if zacc,
                %%but not if contrast (M>R), since that would be double dipping:
                
                if i>funMask  %%1; set to 5 if just ANAT
                    M.data=M.data .* MaskFun;
                    M.data(M.data==0)=NaN;
                end;
                
                %%%%***Take only the highest 20%***
                %                 dataMax=max(M.data);
                %                 dataMax=repmat(dataMax,length(M.data),1);
                %                 percMax=M.data ./ dataMax;
                %                 percMask=percMax>0.7;
                %                 M.data=M.data .* percMask;
                %                 M.data(M.data==0)=NaN;
                
                %Mean ROI:
                C.([varname{i} '_mean'])=mean(M.data(Mask,:))';
                %                 C.([varname{i} '_mean'])=nanmean(M.data(Mask,:))';
                
                
                
                D= addstruct(D,C,'column');
                C=[];
            end;
            
            A=addstruct(A,D);
            D=[];
        end;
        
        
        
        if dim>0
            %%%Dimensionality analysis:
            Y=[A.dim1_mean;A.dim2_mean;A.dim3_mean;A.dim4_mean;A.dim5_mean;A.dim6_mean;A.dim7_mean;A.dim8_mean];
            Hemi=repmat(A.hem,8,1);
            Cond1=ones(length(A.SN),1);
            Cond2=ones(length(A.SN),1)*2;
            Cond3=ones(length(A.SN),1)*3;
            Cond4=ones(length(A.SN),1)*4;
            Cond5=ones(length(A.SN),1)*5;
            Cond6=ones(length(A.SN),1)*6;
            Cond7=ones(length(A.SN),1)*7;
            Cond8=ones(length(A.SN),1)*8;
            Cond=[Cond1;Cond2;Cond3;Cond4;Cond5;Cond6;Cond7;Cond8];
            Subj=(1:length(subj_name))';
            Subj=repmat(Subj,16,1);
            
            figure;
            lineplot([Hemi Cond],Y,'style_thickline');
            ylim([-0.1 1.5]);
            title(ROI_name);
        end;
        
        
        %%repeated measures ANOVA: COND(inter/ord/temp) X HEMI(left/right)
        %  Y=[A.zacc_comb_mean;A.zacc_ord_mean;A.zacc_temp_mean];
        %  Hemi=repmat(A.hem,3,1);
        %  Cond1=ones(length(A.SN),1);
        %  Cond2=ones(length(A.SN),1)*2;
        %  Cond3=ones(length(A.SN),1)*3;
        %  Cond=[Cond1;Cond2;Cond3];
        %  Subj=(1:16)';
        %  Subj=repmat(Subj,6,1);
        %  stats = rm_anova2(Y,Subj,Cond,Hemi,{'Cond','Hemi'})
        %
        %  %%Joern's function:
        %  Factors=[Cond,Hemi];
        %  results=anovaMixed(Y,Subj,'within',Factors,{'condition','hemi'})
        %
        %
        %  %%repeated measures ANOVA: COND(ord/temp) X HEMI(left/right)
        %   Y=[A.zacc_ord_mean;A.zacc_temp_mean];
        %  Hemi=repmat(A.hem,2,1);
        %  Cond1=ones(length(A.SN),1);
        %  Cond2=ones(length(A.SN),1)*2;
        %  Cond=[Cond1;Cond2];
        %  Subj=(1:16)';
        %  Subj=repmat(Subj,4,1);
        %  stats = rm_anova2(Y,Subj,Cond,Hemi,{'Cond','Hemi'})
        %
        %  %%Joern's function:
        %  Factors=[Cond,Hemi];
        %  results=anovaMixed(Y,Subj,'within',Factors,{'condition','hemi'})
        %
        %  fontSize=14;
        
        
        %         figure;
        %
        %         a(1)=subplot(1,5,1);
        %         barplot([A.hem],A.ttest_mean,'split',A.hem); %%signal change M>R
        %         axis([0 3.5 -0.5 1.5]); title('M>R','FontSize',fontSize); ylabel('% sc (mean)','FontSize',fontSize);
        %
        %         a(2)=subplot(1,5,2);
        %         barplot([A.hem],A.zacc_combU_mean,'split',A.hem); title('Overall','FontSize',fontSize); ylabel('Zacc (mean)','FontSize',fontSize);
        %         a(3)=subplot(1,5,3);
        %         barplot([A.hem],A.zacc_comb_mean,'split',A.hem); title('Integrated','FontSize',fontSize); ylabel(' ','FontSize',fontSize);
        %         a(4)=subplot(1,5,4);
        %         barplot([A.hem],A.zacc_ord_mean,'split',A.hem); title('Ordinal','FontSize',fontSize); ylabel(' ','FontSize',fontSize);
        %         a(5)=subplot(1,5,5);
        %         barplot([A.hem],A.zacc_temp_mean,'split',A.hem); title('Temporal','FontSize',fontSize); ylabel(' ','FontSize',fontSize);
        %
        %
        %
        %         linkaxes(a(2:5),'xy');
        %          axis([0 3.5 -0.5 1.5]);
        
        
        %%organized by hemi across conditions;
        C=[];
        for i=1:15
            
            if i==1
                D=tapply(A,{'SN','hem'},{A.zacc_comb_mov_mean,'nanmean','name','zacc'});
            elseif i==2
                D=tapply(A,{'SN','hem'},{A.zacc_comb_prep_mean,'nanmean','name','zacc'});
            elseif i==3
                D=tapply(A,{'SN','hem'},{A.zacc_int_mov_mean,'nanmean','name','zacc'});
            elseif i==4
                D=tapply(A,{'SN','hem'},{A.zacc_int_prep_mean,'nanmean','name','zacc'});
            elseif i==5
                D=tapply(A,{'SN','hem'},{A.zacc_spat_mov_mean,'nanmean','name','zacc'});
            elseif i==6
                D=tapply(A,{'SN','hem'},{A.zacc_spat_prep_mean,'nanmean','name','zacc'});
            elseif i==7
                D=tapply(A,{'SN','hem'},{A.zacc_int_prep_mean,'nanmean','name','zacc'});
            elseif i==8
                D=tapply(A,{'SN','hem'},{A.zacc_temp_mov_mean,'nanmean','name','zacc'});
            elseif i==9
                D=tapply(A,{'SN','hem'},{A.zacc_temp_prep_mean,'nanmean','name','zacc'});
            elseif i==10
                D=tapply(A,{'SN','hem'},{A.ttest_mean(:,1),'nanmean','name','zacc'});
            elseif i==11
                D=tapply(A,{'SN','hem'},{A.ttest_mean(:,2),'nanmean','name','zacc'});
            elseif i==12
                D=tapply(A,{'SN','hem'},{A.ttest_mean(:,3),'nanmean','name','zacc'});
            elseif i==13
                D=tapply(A,{'SN','hem'},{A.ttest_mean(:,4),'nanmean','name','zacc'});
            elseif i==14
                D=tapply(A,{'SN','hem'},{A.ttest_mean(:,5),'nanmean','name','zacc'});
            elseif i==15
                D=tapply(A,{'SN','hem'},{A.ttest_mean(:,6),'nanmean','name','zacc'});
            end;
            
            cond=i; D.cond=repmat(cond,length(D.zacc),1);
            C=addstruct(C,D);
            D=[];
        end;
        figure;
        barplot([C.hem],C.zacc,'split',C.cond,'leg',{'Overall','Timing','Order','Integrated'},'style_bold',...
            'errorcolor',{[0 0 0],[1 0 0],[0 0 1],[0 1 0]},'facecolor',{[0 0 0],[1 0 0],[0 0 1],[0 1 0]},...
            'edgecolor',{[0 0 0],[1 0 0],[0 0 1],[0 1 0]},'gapwidth',[1 0.3 0.3]);
        
        %         barplot([C.hem],C.zacc,'split',C.cond,'style_bold','plotfcn','mean',...
        %             'errorcolor',{[0 0 0],[1 0 0],[0 0 1],[0 1 0]},'facecolor',{[0 0 0],[1 0 0],[0 0 1],[0 1 0]},...
        %             'edgecolor',{[0 0 0],[1 0 0],[0 0 1],[0 1 0]},'gapwidth',[1 0.3 0.3]);
        
        fontSize=12;
        ylabel('Zacc (mean)','FontSize',fontSize);
        axis([0 12 0 1.5]);
        
        %         title(RegName);
        
        set(gcf,'PaperPosition',[1 1 5 2]); %just played around with the figure size
        
        wysiwyg;
        
        
        bonf=4; %Bonferroni correction
        tailsTest=1;
        
        T=tapply(A,{'SN','hem'},{A.zacc_comb_mov_mean,'nanmean','name','zacc'},'subset',A.hem==1);
        [t_combMov_left,p_combU_left]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        T=tapply(A,{'SN','hem'},{A.zacc_comb_mov_mean,'nanmean','name','zacc'},'subset',A.hem==2);
        [t_combMov_right,p_combU_right]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        
        T=tapply(A,{'SN','hem'},{A.zacc_comb_prep_mean,'nanmean','name','zacc'},'subset',A.hem==1);
        [t_combPrep_left,p_comb_left]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        T=tapply(A,{'SN','hem'},{A.zacc_comb_prep_mean,'nanmean','name','zacc'},'subset',A.hem==2);
        [t_combPrep_right,p_comb_right]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        
        T=tapply(A,{'SN','hem'},{A.zacc_int_mov_mean,'nanmean','name','zacc'},'subset',A.hem==1);
        [t_ord_left,p_ord_left]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        T=tapply(A,{'SN','hem'},{A.zacc_int_mov_mean,'nanmean','name','zacc'},'subset',A.hem==2);
        [t_ord_right,p_ord_right]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        
        T=tapply(A,{'SN','hem'},{A.zacc_int_prep_mean,'nanmean','name','zacc'},'subset',A.hem==1);
        [t_temp_left,p_temp_left]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        T=tapply(A,{'SN','hem'},{A.zacc_int_prep_mean,'nanmean','name','zacc'},'subset',A.hem==2);
        [t_temp_right,p_temp_right]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        
        T=tapply(A,{'SN','hem'},{A.zacc_temp_mov_mean,'nanmean','name','zacc'},'subset',A.hem==1);
        [t_temp_left,p_temp_left]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        T=tapply(A,{'SN','hem'},{A.zacc_temp_mov_mean,'nanmean','name','zacc'},'subset',A.hem==2);
        [t_temp_right,p_temp_right]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        
        T=tapply(A,{'SN','hem'},{A.zacc_temp_prep_mean,'nanmean','name','zacc'},'subset',A.hem==1);
        [t_temp_left,p_temp_left]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        T=tapply(A,{'SN','hem'},{A.zacc_temp_prep_mean,'nanmean','name','zacc'},'subset',A.hem==2);
        [t_temp_right,p_temp_right]=ttest(T.zacc,T.zacc,tailsTest,'onesample');
        
        stats=[t_combU_left,p_combU_left*bonf,t_combU_right,p_combU_right*bonf;...
            t_comb_left,p_comb_left*bonf,t_comb_right,p_comb_right*bonf;...
            t_ord_left,p_ord_left*bonf,t_ord_right,p_ord_right*bonf;...
            t_temp_left,p_temp_left*bonf,t_temp_right,p_temp_right*bonf];
        disp(ROI_name);
        disp(stats);
        
        % save(fullfile(regDir,'area_stat.mat'),'-struct','D');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Correlation with Behaviour!
        cd(fullfile(behDir));
        Transfer=dload('transfer.txt');
        Transfer.comb=Transfer.comb*1000;
        Transfer.ord=Transfer.ord*1000;
        Transfer.temp=Transfer.temp*1000;
        
        %%Mean across hemis:
        T.combU=(A.zacc_combU_mean(1:16)+A.zacc_combU_mean(17:32))/2;
        T.combU(1,:)=[];
        T.ord=(A.zacc_ord_mean(1:16)+A.zacc_ord_mean(17:32))/2;
        T.ord(1,:)=[];
        T.temp=(A.zacc_temp_mean(1:16)+A.zacc_temp_mean(17:32))/2;
        T.temp(1,:)=[];
        %
        %         %%Left hemi:
        %         T.combU=A.zacc_combU_mean(1:16);
        %         T.combU(1,:)=[];
        %         T.ord=A.zacc_ord_mean(1:16);
        %         T.ord(1,:)=[];
        %         T.temp=A.zacc_temp_mean(1:16);
        %         T.temp(1,:)=[];
        %
        %%Right hemi:
        T.combU=A.zacc_combU_mean(17:32);
        T.combU(1,:)=[];
        T.ord=A.zacc_ord_mean(17:32);
        T.ord(1,:)=[];
        T.temp=A.zacc_temp_mean(17:32);
        T.temp(1,:)=[];
        
        %
        bonf=2; %bonferroni correction
        figure;
        b(1)=subplot(1,3,1); [r_combU,p_combU]=corr(T.combU,Transfer.comb,'tail','right')
        scatterplot(T.combU,Transfer.comb,'regression','linear'); title('Overall','FontSize',fontSize);
        ylabel('Transfer (diff. to bsl in ms)','FontSize',fontSize);
        xlabel('Zacc (mean)','FontSize',fontSize);
        
        b(2)=subplot(1,3,2); [r_ord,p_ord]=corr(T.ord,Transfer.ord,'tail','right')
        scatterplot(T.ord,Transfer.ord,'regression','linear');
        title('Ord','FontSize',fontSize);
        ylabel('Transfer (diff. to bsl in ms)','FontSize',fontSize);
        xlabel('Zacc (mean)','FontSize',fontSize);
        
        b(3)=subplot(1,3,3); [r_temp,p_temp]=corr(T.temp,Transfer.temp,'tail','right')
        scatterplot(T.temp,Transfer.temp,'regression','linear');
        title('Temp','FontSize',fontSize);
        ylabel('Transfer (diff. to bsl in ms)','FontSize',fontSize);
        xlabel('Zacc (mean)','FontSize',fontSize);
        
        stats=[r_combU,p_combU*bonf;...
            r_ord,p_ord*bonf;...
            r_temp,p_temp*bonf];
        disp(stats);
        
        linkaxes(b,'xy');
        axis([-1 3 -20 80]);
        set(gcf,'PaperPosition',[1 1 10 2]); %just played around with the figure size
        
        wysiwyg;
        
    case 'surf_overlap'
        zthresh=1.97;
        O=caret_load('cSPM_srh.zacc_OrdOneout_160.metric');
        T=caret_load('cSPM_srh.zacc_TempOneout_160.metric');
        A=O;
        A.data=zeros(length(O.data),15); %define structure
        %A.data(O.data(:,15)>zthresh&T.data(:,15)>zthresh,15)=(O.data(:,15)+T.data(:,15))/2; %smaller Z-value
        
        %A.data=(O.data(O.data(:,15)>zthresh&T.data(:,15)>zthresh,15)+T.data(O.data(:,15)>zthresh&T.data(:,15)>zthresh,15))/2; %ave Z-value
        aboveThresh=find(O.data(:,15)>zthresh&T.data(:,15)>zthresh);
        A.data(aboveThresh,15)=(O.data(aboveThresh,15)+T.data(aboveThresh,15))/2; %ave Z-value
        
        caret_save('cSPM_srh.zacc_overlapTempOrd_160.metric',A)
        
    case 'ROI_make_BAs_paint' % Define ROIS from probabilistic atlas
        for h=1:2 %out of memory error, but doing them separately helps
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.propatlas.metric']));
            
            % M1
            i=(strcmp('BA4a',C.column_name)  | strcmp('BA4p',C.column_name));
            D(:,1)=sum(C.data(:,i),2);
            
            % S1
            i=(strcmp('BA6',C.column_name));
            D(:,2)=sum(C.data(:,i),2);
            
            % V1
            i=(strcmp('V1',C.column_name));
            D(:,3)=sum(C.data(:,i),2);
            
            [a,b]=max(D,[],2);
            b(a<0.35)=0;
            P=caret_struct('paint','data',b,'paintnames',{'???','M1','PM','V1'});
            caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.ROIs_BAs.paint']),P);
            
            
            %             [a,b]=max(D,[],2);
            %             b(a<0.35)=0;
            %             P=caret_struct('paint','data',b,'paintnames',{'???','M1'});
            %             caret_save(fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.M1.paint']),P);
            
            
            
            
        end;
        
    case 'ROI_BA_define'  % Project these ROIs into individual space
        s=varargin{1};
        numregions=length(regname);
        for h=1:2
            C=caret_load(fullfile(caretDir,'fsaverage_sym',hemName{h},[hem{h} '.ROIs_BAs.paint']));
            caretSubjDir=fullfile(caretDir,['x' subj_name{s}]);
            file=fullfile(glmDir,subj_name{s},'mask.nii');
            for i=1:numregions
                R{i+(h-1)*numregions}.type='surf_nodes';
                R{i+(h-1)*numregions}.location=find(C.data==i);
                R{i+(h-1)*numregions}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                R{i+(h-1)*numregions}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                R{i+(h-1)*numregions}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                R{i+(h-1)*numregions}.linedef=[5,0,1];
                R{i+(h-1)*numregions}.image=file;
                R{i+(h-1)*numregions}.name=[subj_name{s} '_' regname{i} '_' hem{h}];
            end;
            % Cerebellum
            % if (~isnan(coord{s,h}(1)));
            %     R{h*4}=region('sphere',coord{s,h},9,[subj_name{s} '_' name{4} '_' hem{h}]);
            % end;
        end;
        R=region_calcregions(R);
        cd(regDir);
        for r=1:length(R);
            if (~isempty(R{r}))
                region_saveasimg(R{r},spm_vol(file));
            end;
        end;
        save([subj_name{s} '_regions.mat'],'R');
        varargout={R};
        
    case 'ROI_data'     % Extracts data for the 3*2 ROIs from the volume
        sn=varargin{1};
        T=[];
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            E.tt=repmat([1:4],1,numRuns(s));
            load(fullfile(regDir,[subj_name{s} '_regions.mat']));
            
            % Pgrayprop=fullfile(baseDir,'anatomicals',subj_name{s},['c1' subj_name{s} '_anatomical.nii']);
            P={};
            N=length(SPM.Vbeta)-numRuns(s);
            for i=1:N
                P{i}=sprintf('beta_%4.4d.img',i);
            end;
            V=spm_vol_RY(char(P));
            data = region_getdata(V,R);
            
            % T/FContrasts
            P={};
            P{1}='spmT_0006.img';
            P{2}='spmF_0007.img';
            V=spm_vol_RY(char(P));
            stats = region_getdata(V,R);
            
            for i=1:length(R)
                if (~isempty(R{i}))
                    D=[];
                    vec=ones(size(data{i},2),1);
                    D.beta=ones(length(vec),32)*NaN;
                    D.beta(1:length(vec),1:N)=data{i}';
                    D.SN=vec*s;
                    D.regNum=vec*i;
                    D.xyz=R{i}.data;
                    for t=1:max(E.tt);
                        D.meanAct(:,t)=mean(D.beta(:,E.tt==t),2);
                    end;
                    D.statsTslope=stats{i}(1,:)';
                    D.statsF=stats{i}(2,:)';
                    T=addstruct(T,D);
                end
            end;
        end;
        T.tt=E.tt;
        save(fullfile(regDir,['reg_data.mat']),'-struct','T');
        varargout={T};
        
    case 'surf_makeRGB'
        for h=1:2
            groupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
            C=caret_load([groupDir '/' hem{h} '.summary.' dataext '_' numvox '.metric']);
            
            RGBdata(:,3)=C.data(:,5); % 1-3 are left only
            RGBdata(:,4)=C.data(:,6); % 4-6 are right only
            RGBdata(:,8)=C.data(:,7); % 7-9 are both only
            RGBdata(:,12)=C.data(:,8); % 10-12 are both only
            RGBdata(:,13)=C.data(:,9); % 13-15 are both only
            RGBdata(:,15)=0;
            RGBdata(:,16)=C.data(:,7+h); % 16-18 are both only
            RGBdata(:,18)=C.data(:,4+h); % 16-18 are both only
            
            thresh=1.1;
            RGBdata(RGBdata<thresh)=0;
            sc=[thresh 3;thresh 3;thresh 3];
            C=caret_struct('RGBpaint','data',RGBdata,'scales',{sc,sc,sc,sc,sc,sc},'column_name',{'left','right','both','BL','BR'});
            caret_save([groupDir '/' hem{h} '.accuracy.' dataext '_' numvox '.RGB_paint'],C);
        end
        
    case 'surf_crosssection'
        nodeS=[65446 59843;32184 36587];
        T=[];
        Side={'left','right'};
        kind={'zacc','con'};
        for h=1:2
            groupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
            cd(groupDir);
            for s=1:2
                
                D.hem=ones(5,1)*h;
                D.SN=[1:5]';
                D.Side=ones(5,1)*s;
                for k=1:2
                    [M]=caret_crosssection(nodeS(h,:),7,[hem{h} '.FLAT.coord'],[hem{h} '.' kind{k} '_' Side{s} '.metric']);
                    %                     [M,x]=caret_crosssection(nodeS(h,:),7,[hem{h} '.FLAT.coord'],[hem{h} '.' kind{k} '_' Side{s} '.metric']);
                    D.(kind{k})=M;
                end
                T=addstruct(T,D);
            end;
            [M,x]=caret_crosssection(nodeS(h,:),7,[hem{h} '.FLAT.coord'],[hem{h} '.surface_shape']);
            T.x=repmat(x,10,1);
            T.sulcaldepth=repmat(M(2,:),10,1);
            
            figure(h);
            subplot(3,1,1);
            traceplot([1:50],T.sulcaldepth);
            subplot(3,1,2);
            traceplot([1:50],T.zacc,'split',T.Side,'errorfcn','stderr','subset',T.hem==h);
            set(gca','YLim',[0 12]);
            subplot(3,1,3);
            traceplot([1:50],T.con,'split',T.Side,'errorfcn','stderr','subset',T.hem==h);
            set(gca','YLim',[-3 7]);
            drawline(0,'dir','horz');
        end;
        
        varargout={T};
        
    case 'surf_a'
        
        %       for h=1:2
        %           areasExclude=[0,1,3]; %include only area 2 (premotor cortex)
        %
        %           groupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
        %           cd(groupDir);
        %
        %           P=caret_load([hem{h} '.ROIs_BAs.paint']);
        %
        %           for i=1:length(areasExclude)
        %               P.data(P.data==areasExclude(i))=NaN; %region number
        %           end;
        %
        %             caret_savepaint([hem{h} '.PM.paint'],P);
        %       end;
        
        
        for h=1:2
            
            dataCol=[2,2;5,4;4,3]; %SMA;PMd;PMv
            areaName={'SMA','PMd','PMv'};
            
            %%from PM_lateral
            dataCol=[2,2;3,3]; %SMA;PMd;PMv
            areaName={'PMd_new','PMv_new'};
            
            %     for area=1:3
            for area=1:2
                %         areasExclude=(0:6); %include only areas spec in index
                %         areasExclude=areasExclude(areasExclude~=index(area,h));
                
                groupDir=[caretDir filesep 'fsaverage'  filesep hemName{h} ];
                cd(groupDir);
                
                %         P=caret_load([hem{h} '.PMall.paint']);
                P=caret_load([hem{h} '.PM_lateral.paint']);
                P.data=P.data(:,dataCol(area,h));
                
                %         for i=1:length(areasExclude)
                %             P.data(P.data==areasExclude(i))=0; %region number
                %         end;
                
                P.column_name={'Column_01'};
                P.paintnames=areaName(area);
                P.num_paintnames=1;
                P.num_cols=1;
                
                caret_savepaint([hem{h} '.' areaName{area} '.paint'],P);
            end;
        end;
        
    case 'surf_ROI_extract'
        atlas=2;
        anaSubj=[1:length(subj_name)]; %subjects used in the group analysis
        
        %         INname={'accuracy_Comb_160Uncorrected','accuracy_Comb_160','accuracy_OrdOneout_160','accuracy_TempOneout_160'};
        %         OUTname={'acc_Comb_160Uncorrected','acc_Comb_160','acc_OrdOneout_160','acc_TempOneout_160'};
        
        INname={'accuracy_Comb_160Uncorrected'};
        OUTname={'acc_Comb_160Uncorrected'};
        
        for h=1:2
            surfaceGroupDir=[caretDir filesep atlasname{atlas} filesep hemName{h} ];
            cd(surfaceGroupDir);
            for j=1:length(INname); %----loop over each input metric file and make a group metric file
                k=0;
                for i=anaSubj; %----define names of subj metric files
                    k=k+1;
                    infilenames{k}=[caretDir filesep 'x' subj_name{i} filesep hemName{h} filesep subj_name{i} '_' INname{j} '.metric'];
                end;
                
            end;
            
            R.type='paint'
            R.location=fullfile([hem{h} '.PM.paint']);
            
            T=caret_regionstat('images',infilenames,'regions',R)
            T
        end;
        
    case 'surf_stat_mask'
        atlas=2;
        for h=1:2
            areasExclude=[1,5];
            
            groupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
            cd(groupDir);
            
            P=caret_load([hem{h} '.desikan.paint']);
            
            
            for i=1:length(areasExclude)
                P.data(P.data==areasExclude(i))=NaN; %region number
                %P.index=P.index(P.data~=areasExclude(i)); %location index
            end;
            
            %P.num_rows=length(P.data);
            caret_savepaint([hem{h} '.mask.paint'],P);
        end
        
    case 'surf_stat_list'
        
        %previously: caret_list(S,cSPM,0.005,0.05,'label','lh.desikan.paint')
        
        % Produces all statistical maps - corrected
        atlas=2;
        %         contraIpsi={'Contra','Ipsi'};
        
        %         SPMname={'con','zacc_Comb_160Uncorrected','zacc_Comb_160','zacc_TempOneout_160','zacc_OrdOneout_160','zacc_Comb_160UncorrectedMeanSearchlight'};
        SPMname={'_zacc_Comb_160_Mov'... %1
            '_zacc_Comb_160_Prep'...     %2
            '_zacc_Int_160_Mov'...       %3
            '_zacc_Int_160_Prep'...      %4
            '_zacc_Spat_160_Mov'...      %5
            '_zacc_Spat_160_Prep'...     %6
            '_zacc_Temp_160_Mov'...      %7
            '_zacc_Temp_160_Prep',...    %8
            'con_0001',...               %9
            'con_0002',...               %10
            'con_0003',...               %11
            'con_0004',...               %12
            'con_0005',...               %13
            'con_0006'};                 %14
        
        for i=1:length(SPMname) %to loop through all metric files
            %
            %            i = varargin{1};
            
            %             u=0.01;        % JD Uncorrected threshold
            %             k=0.2;          % JD Cluster threshold
            
            u=0.001;        % Uncorrected threshold
            k=0.05;
            
            
            for h=1:2
                groupDir=[caretDir filesep atlasname{atlas}  filesep hemName{h} ];
                cd(groupDir);
                
                %             P=caret_load([hem{h} '.lobes.paint']); %JD
                %             searchArea=(P.data==1 | P.data==2); % Frontal or parietal %JD
                
                S = load([hem{h} '.avrgsurface.mat']);
                P=caret_load([hem{h} '.mask.paint']);
                searchArea=(P.data>0); % Frontal or parietal
                
                
                load(['cSPM_s', hem{h}, '.',  SPMname{i}, '.mat']);
                
                %             load([hem{h} '.avrgsurface.mat']);
                %             load(['cSPM_s', contraIpsi{h}, '.',  SPMname{c}, '.mat']);
                
                %             M=caret_load([hem{h} '.desikan.paint']);
                %                 M=caret_load(['ROI.paint']);
                
                %                         caret_list(S,cSPM,u,k,'contrast',1,'sign',1,'mask',searchArea,...
                %                             'sort_by','area');
                %             caret_list(S,cSPM,u,k,'contrast',1,'sign',1,'mask',searchArea,...
                %                 'sort_by','area');
                disp([SPMname{i} ' ' hem{h}])
                caret_list_RY(S,cSPM,u,k,'contrast',1,'sign',1,'mask',searchArea,...
                    'sort_by','area');
                
                %                             caret_list_RY(S,cSPM,u,k,'contrast',1,'sign',1,'mask',searchArea,...
                %                                 'label',M,'sort_by','area');
            end
        end
        
    case 'surf_results'
        % Instructions on how to use the caret_toolbox for multiple comp.
        % correction
        h=1;
        conname='zacc_diff';
        maskname='zacc_avrg';
        u=3;
        k=0.2;
        u_mask=1.64; % p<0.05 uncorrected
        
        % FIRST MAKE SURFACE (DO THIS BEFORE)
        surfaceGroupDir=[caretDir filesep 'fsaverage_sym'  filesep hemName{h} ];
        cd(surfaceGroupDir)
        load([hem{h} '.surface.mat']);
        % Fist load the masking image and get the mask
        load(['cSPM_' maskname '.mat']);
        mask=cSPM.con(1).Z>u_mask;
        % Now get the contrast and make the list
        load(['cSPM_' conname '.mat']);
        T=caret_list(S,cSPM,u,k,'sort_by','p_cluster','mask',mask,...
            'save_thresholded',[hem{h} '.summary.metric,17']);
        
    case 'make_mask'  % Makes restricted analysis mask for MVA
        s=varargin{1};
        mask=fullfile(glmDir, subj_name{s},'mask.img');
        omask=fullfile(glmDir, subj_name{s},'maskbrain.nii');
        P1=fullfile(baseDir,'anatomicals',subj_name{s},sprintf('c1%s_anatomical.nii',subj_name{s}));
        P2=fullfile(baseDir,'anatomicals',subj_name{s},sprintf('c2%s_anatomical.nii',subj_name{s}));
        sP1=fullfile(baseDir,'anatomicals',subj_name{s},sprintf('sc1%s_anatomical.nii',subj_name{s}));
        sP2=fullfile(baseDir,'anatomicals',subj_name{s},sprintf('sc2%s_anatomical.nii',subj_name{s}));
        spm_smooth(P1,sP1,[4 4 4]);
        spm_smooth(P2,sP2,[4 4 4]);
        spm_imcalc_ui({mask,sP1,sP2},omask,'i1 & (i2 +i3)>0.01',{});
        
    case 'make_mask_suit'
        s=varargin{1};
        mask=fullfile(glmDir, subj_name{s},'mask.img');
        suit=fullfile(anatDir, subj_name{s},['c_', subj_name{s}, '_anatomical_pcereb.nii']);
        omask=fullfile(suitDir, subj_name{s},'maskbrainSUIT.nii');
        
        spm_imcalc_ui({mask,suit},omask,'i1>0 & i2>0.7',{});
        
    case 'make_mask_BG' %Basal ganglia mask (based on FSL/FIRST)
        s=varargin{1};
        mask=fullfile(glmDir, subj_name{s},'mask.img');
        BG=fullfile(anatDir, subj_name{s},['fslfirst_', subj_name{s}, '_subcort_all_fast_firstseg.nii']);
        
        % %         %isolate individual BG structures
        L_Caud=fullfile(anatDir, subj_name{s},'L_Caud.nii');
        %          spm_imcalc_ui(BG,L_Caud,'i1>10 & i1<12',{});
        
        R_Caud=fullfile(anatDir, subj_name{s},'R_Caud.nii');
        %          spm_imcalc_ui(BG,R_Caud,'i1>49 & i1<51',{});
        
        L_Put=fullfile(anatDir, subj_name{s},'L_Put.nii');
        %          spm_imcalc_ui(BG,L_Put,'i1>11 & i1<13',{});
        
        R_Put=fullfile(anatDir, subj_name{s},'R_Put.nii');
        %          spm_imcalc_ui(BG,R_Put,'i1>50 & i1<52',{});
        
        L_Pallid=fullfile(anatDir, subj_name{s},'L_Pallid.nii');
        %          spm_imcalc_ui(BG,L_Pallid,'i1>12 & i1<14',{});
        
        R_Pallid=fullfile(anatDir, subj_name{s},'R_Pallid.nii');
        %          spm_imcalc_ui(BG,R_Pallid,'i1>51 & i1<53',{});
        
        
        %%%%%%%%%% Put together into one BG image:
        %         BGall=fullfile(anatDir, subj_name{s},'BGall.nii');
        %          spm_imcalc_ui({L_Caud,R_Caud,L_Put,R_Put,L_Pallid,R_Pallid},BGall,'i1>0 | i2>0 | i3>0 | i4>0 | i5>0 | i6>0',{});
        
        
        %          % Put together into one BG image:
        BGleft=fullfile(anatDir, subj_name{s},'BGleft.nii');
        spm_imcalc_ui({L_Caud,L_Put,L_Pallid},BGleft,'i1>0 | i2>0 | i3>0',{});
        
        
        % Put together into one BG image:
        BGright=fullfile(anatDir, subj_name{s},'BGright.nii');
        spm_imcalc_ui({R_Caud,R_Put,R_Pallid},BGright,'i1>0 | i2>0 | i3>0',{});
        
        
        %%%Leave out Pallidum&Caudatus:
        % Put together into one BG image:
        
        %          BGall=fullfile(anatDir, subj_name{s},'BGall.nii');
        %           spm_imcalc_ui({L_Put,R_Put},BGall,'i1>0 | i2>0',{});
        
        
        %          BGleft=fullfile(anatDir, subj_name{s},'BGleft.nii');
        %          spm_imcalc_ui({L_Put},BGleft,'i1>0',{});
        %
        %
        %          % Put together into one BG image:
        %          BGright=fullfile(anatDir, subj_name{s},'BGright.nii');
        %          spm_imcalc_ui({R_Put},BGright,'i1>0',{});
        %%%
        
        
        %          % Mask with EPI
        %          omask=fullfile(glmDir, subj_name{s},'maskbrainBG.nii');
        %          spm_imcalc_ui({mask,BGall},omask,'i1>0 & i2>0',{});
        
        % Mask with EPI
        omask=fullfile(glmDir, subj_name{s},'maskbrainBGleft.nii');
        spm_imcalc_ui({mask,BGleft},omask,'i1>0 & i2>0',{});
        
        
        % Mask with EPI
        omask=fullfile(glmDir, subj_name{s},'maskbrainBGright.nii');
        spm_imcalc_ui({mask,BGright},omask,'i1>0 & i2>0',{});
        
    case 'MNI_make_rfx' %specify scans for the one-samples t-test, e.g. z-transformed MVA images
        %sl4_imana('MNI_make_rfx',[1:3])
        sn=varargin{1};
        %        folders={'data_Temp_oneout','data_Ord_oneout','data_Comb_corrected','data_Comb_Uncorrected'};
        foldersin={'data_Temp_oneout','data_Ord_oneout','data_Comb_Uncorrected'};
        
        %folders={'data_Temp_oneout/Right','data_Ord_oneout/Right','data_Comb_Uncorrected/Right'};
        folders={'data_Temp_oneout/Left','data_Ord_oneout/Left','data_Comb_Uncorrected/Left'};
        %         images= {'smva160_zacc_Temp_oneout','smva160_zacc_Ord_oneout','smva160_zacc_comb'};
        %        images= {'smva160_zacc_Temp_oneoutBG','smva160_zacc_Ord_oneoutBG','smva160_zacc_combCorrectedBG','smva160_zacc_combUncorrectedBG'};
        %         images= {'smva100_zacc_Temp_oneoutBG','smva100_zacc_Ord_oneoutBG','smva100_zacc_combUncorrectedBG'};
        
        
        %images= {'smva160_zacc_Temp_oneoutBGright','smva160_zacc_Ord_oneoutBGright','smva160_zacc_combUncorrectedBGright'};
        images= {'smva160_zacc_Temp_oneoutBGleft','smva160_zacc_Ord_oneoutBGleft','smva160_zacc_combUncorrectedBGleft'};
        
        %         mask=fullfile(groupDir,'data/mask_thresh.nii');
        %mask={fullfile(groupDir,'data_BG/mask_thresright.nii,1')};
        mask={fullfile(groupDir,'data_BG/mask_thresleft.nii,1')};
        
        nam=[];
        for j=1:numel(images)
            %mkdir(fullfile(groupDir,[images{j}])); % folder for each contrast
            matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(groupDir,'data_BG',[folders{j}])} %,'_4reg'
            scans=[];
            for s=sn:18
                scans{end+1} = fullfile(groupDir, [foldersin{j}], [images{j},'_',subj_name{s}, '.nii,1']) %_4reg
            end
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans= scans;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = mask;
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            spm_jobman('run',matlabbatch);
            
        end;
        
    case 'MNI_rfx_spm' %do one-samples t-test
        %sl4_imana('MNI_rfx_spm', 1:2)
        %con={'data_Temp_oneout','data_Ord_oneout','data_Comb_corrected','data_Comb_Uncorrected'};
        %con={'data_Temp_oneout','data_Ord_oneout','data_Comb_Uncorrected'};
        %con={'data_Temp_oneout/Right','data_Ord_oneout/Right','data_Comb_Uncorrected/Right'};
        con={'data_Temp_oneout/Left','data_Ord_oneout/Left','data_Comb_Uncorrected/Left'};
        
        contr=varargin{1};
        for k=contr
            dirname=fullfile(groupDir,'data_BG',[con{k}]);
            cd(dirname);
            load SPM;
            SPM=spm_spm(SPM);
            spmj_rfx_contrast(SPM);
            cd ..
        end;
        
    case 'SUIT_make_rfx' %specify scans for the one-samples t-test, e.g. z-transformed MVA images
        %sl4_imana('MNI_make_rfx',[1:3])
        sn=varargin{1};
        folders={'data_Temp_oneout','data_Ord_oneout','data_Comb_corrected','data_Comb_Uncorrected'};
        %         images= {'smva160_zacc_Temp_oneout','smva160_zacc_Ord_oneout','smva160_zacc_comb'};
        images= {'swcmmva160_zacc_Temp_oneout','swcmmva160_zacc_Ord_oneout','swcmmva160_zacc_comb_oneout','swcmmva160_zacc_comb_Uncorrected'};
        %         mask=fullfile(groupDir,'data/mask_thresh.nii');
        mask={fullfile(suitDir,'mask_threshGrey.nii,1')}; %group mask
        
        nam=[];
        for j=1:numel(images)
            %mkdir(fullfile(groupDir,[images{j}])); % folder for each contrast
            matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(suitDir,[folders{j}])}; %,'_4reg'
            scans=[];
            for s=sn:18
                scans{end+1} = fullfile(suitDir, [folders{j}], [images{j},'_',subj_name{s}, 'SUIT.nii,1']) %_4reg
            end
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans= scans;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = mask;
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            spm_jobman('run',matlabbatch);
            
        end;
        
    case 'SUIT_rfx_spm' %do one-samples t-test
        %sl4_imana('MNI_rfx_spm', 1:2)
        con={'data_Temp_oneout','data_Ord_oneout','data_Comb_corrected','data_Comb_Uncorrected'};
        contr=varargin{1};
        for k=contr
            dirname=fullfile(suitDir,[con{k}]);
            cd(dirname);
            load SPM;
            SPM=spm_spm(SPM);
            spmj_rfx_contrast(SPM);
            cd ..
        end;
        
    case 'suit_smooth' %%%Smoothing after suit reslice
        
        s=varargin{1};
        
        INfile={'_zacc_Comb_160Corrected_surfvol',...
            '_zacc_Comb_160OrdOneOut_surfvol',...
            '_zacc_Comb_160TempOneOut_surfvol',...
            '_zacc_Comb_160Uncorrected_surfvol'};
        
        outFolder={'CombCorrected',...
            'Ord',...
            'Temp',...
            'Comb'};
        
        %         INfile={'_zacc_Comb_160Corrected_surfvol','_zacc_Comb_160Corrected_surfvol_flip',...
        %             '_zacc_Comb_160OrdOneOut_surfvol','_zacc_Comb_160OrdOneOut_surfvol_flip',...
        %             '_zacc_Comb_160TempOneOut_surfvol','_zacc_Comb_160TempOneOut_surfvol_flip',...
        %             '_zacc_Comb_160Uncorrected_surfvol','_zacc_Comb_160Uncorrected_surfvol_flip'};
        %         outFolder={'CombCorrected','CombCorrected',...
        %             'Ord','Ord',...
        %             'Temp','Temp',...
        %             'Comb','Comb'};
        
        for i=1:length(INfile)
            comb=fullfile(glmDir, subj_name{s}, ['wcm' subj_name{s} INfile{i} '.nii']);
            scomb=fullfile(suitDir, outFolder{i}, ['swcm' subj_name{s} INfile{i} '.nii']);
            spm_smooth(comb,scomb,[4 4 4]);
        end;
        
        
        %         temp=fullfile(glmDir, subj_name{s},['Temp_oneout/wcmmva160_zacc_Temp_oneoutSUIT.nii']);
        %         stemp=fullfile(suitDir, ['data_Temp_oneout/swcmmva160_zacc_Temp_oneout_' subj_name{s} 'SUIT.nii']);
        %         spm_smooth(temp,stemp,[4 4 4]);
        % %
        %         ord=fullfile(glmDir, subj_name{s},['Ord_oneout/wcmmva160_zacc_Ord_oneoutSUIT.nii']);
        %         sord=fullfile(suitDir, ['data_Ord_oneout/swcmmva160_zacc_Ord_oneout_' subj_name{s} 'SUIT.nii']);
        %         spm_smooth(ord,sord,[4 4 4]);
        % %
        %
        %         comb=fullfile(glmDir, subj_name{s},['wcmmva160_zacc_combCorrectedSUIT.nii']);
        %         scomb=fullfile(suitDir, ['data_Comb_corrected/swcmmva160_zacc_Comb_oneout_' subj_name{s} 'SUIT.nii']);
        %         spm_smooth(comb,scomb,[4 4 4]);
        
        
        %         comb=fullfile(glmDir, subj_name{s},['wcmmva160_zacc_combUncorrectedSUIT.nii']);
        %         scomb=fullfile(suitDir, ['data_Comb_Uncorrected/swcmmva160_zacc_Comb_Uncorrected_' subj_name{s} 'SUIT.nii']);
        %         spm_smooth(comb,scomb,[4 4 4]);
        
        
        %     case 'MVA_mask' %%%%OLD
        %
        %         % Make a mask between greay matter and the epi in folder
        %         anatomical_brain4mask %(After registration!), produces better
        %         results than % caret adjustments per hand
        %
        %         %%% Smooth the grey matter image with 2mm to include voxels that
        %         are slightly %%% outside the brain, but do still have grey matter
        %         signal.
        %
        %
        %
        %         s=varargin{1}; cd(fullfile(glmDir, subj_name{s}));
        %
        %         matlabbatch{1}.spm.util.imcalc.input = {
        %             ['E:\Data\MotorControlServer\project\tempord\tempord3\GLM_firstlevel_run\',
        %             subj_name{s}, '\mask.img,1']
        %             ['E:\Data\MotorControlServer\project\tempord\tempord3\anatomicals\',
        %             subj_name{s}, '\', subj_name{s} '_mask4glm.img,1'] };
        %         matlabbatch{1}.spm.util.imcalc.output = 'gmmask.img'; %grey
        %         matter mask matlabbatch{1}.spm.util.imcalc.outdir = {''};
        %         matlabbatch{1}.spm.util.imcalc.expression = 'i2.*(i1>0)';
        %         matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        %         matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        %         matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        %         matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        %
        %
        %         spm_jobman('run',matlabbatch);
        
    case 'MVA_search' % Define the search lights for the MVA analysis
        
        
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            V=spm_vol('maskbrain.nii'); %if preceded by case MVA_mask
            X=spm_read_vols(V);
            [i,j,k]=ind2sub(size(X),find(X~=0));
            vox=[i j k];
            
            
            [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[16 160],V.mat,V.dim,[],'mva160_numvox.nii');
            save volsearch160.mat vox LI voxmin voxmax n
            
        end;
        
    case 'MVA_searchSUIT' % Define the search lights for the MVA analysis
        
        
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            %V=spm_vol('maskbrain.nii'); %if preceded by case MVA_mask
            V=spm_vol('maskbrainSUIT.nii'); %if preceded by case MVA_mask
            X=spm_read_vols(V);
            [i,j,k]=ind2sub(size(X),find(X~=0));
            vox=[i j k];
            
            
            [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[16 160],V.mat,V.dim,[],'mva160_numvoxSUIT.nii');
            save volsearch160SUIT.mat vox LI voxmin voxmax n
            %             [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[16 80],V.mat,V.dim,[],'mva80_numvoxSUIT.nii');
            %             save volsearch80SUIT.mat vox LI voxmin voxmax n
            
        end;
        
    case 'ROI_mean_MVA'
        sn=varargin{1};
        nrregr=10;
        
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            reg={'L_M1sphere.img'};
            
            %         R=region('roi_image',reg{k},1);
            %         R=region_calcregions(R);
            
            R.type='sphere'
            R.location(s)=[-34,-23,59]; %L_M1 M>R %%%% !!!!!!!TO DO: transform to subject space
            R.radius=160;
            
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<10)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<10); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:120).fname}';
            P(2:2:120)=[]; %Leave out betas for 1st derivative
            P(10:10:60)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            
            
            V=spm_vol(char(P));
            Y=region_getdata(V,R);
            
            [D,M]=lmva_randomsubspace(Y',@tempord3_eventclass,'params',{cTemp,run,train,test},...
                'num_iter',1,'size_subspace',160);
            acc(:,1)=M.acc;
            
        end;
        varargout={acc};
        acc=[];
        
    case 'MVA_searchBG' % Define the search lights for the MVA analysis
        
        sn=varargin{1};
        nr=varargin{2};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            %V=spm_vol('maskbrain.nii'); %if preceded by case MVA_mask
            %V=spm_vol('maskbrainBG.nii'); %if preceded by case MVA_mask
            image={'maskbrainBGright.nii','maskbrainBGleft.nii'};
            out={'mva160_numvoxBGright.nii','mva160_numvoxBGleft.nii'};
            V=spm_vol(image{nr}); %if preceded by case MVA_mask
            X=spm_read_vols(V);
            [i,j,k]=ind2sub(size(X),find(X~=0));
            vox=[i j k];
            
            if nr==1
                [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[16 160],V.mat,V.dim,[],out{1});
                save volsearch160BGright.mat vox LI voxmin voxmax n
            else
                [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[16 160],V.mat,V.dim,[],out{2});
                save volsearch160BGleft.mat vox LI voxmin voxmax n
            end;
            
            %             [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[16 160],V.mat,V.dim,[],'mva160_numvoxBG.nii');
            %             save volsearch160BG.mat vox LI voxmin voxmax n
            %             [LI,voxmin,voxmax,n]=lmva_voxelselection(vox(:,:)',vox',[16 80],V.mat,V.dim,[],'mva80_numvoxSUIT.nii');
            %             save volsearch80SUIT.mat vox LI voxmin voxmax n
            
        end;
        
    case 'MVA_BG_random_subspaceTiming' %random subspace - left - right BG
        
        sn=varargin{1};
        nrregr=10;
        
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            reg={'maskbrainBGleft.nii','maskbrainBGright.nii'};
            for k=1:2
                disp(reg{k});
                R=region('roi_image',reg{k},1);
                R=region_calcregions(R);
                
                
                load SPM;
                nrruns=length(SPM.nscan);
                
                % Load SPM info file (to label beta values)
                D=dload('SPM_info.ana'); %% no derivatives
                c=D.condition(D.condition<nrregr)'; % extract conditions (leave out errors!)
                cTemp=floor((c-1)/3)+1;
                run=D.RN(D.condition<nrregr); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
                
                
                % Generate column indices for cross-validation, where
                % cell i contains column indices of the respective test and
                % train set
                
                % Reminder:
                %%% T1_O1 T1_O2 T1_O3 %%% T2_O1 T2_O2 T2_O3 %%% T3_O1 T3_O2 T3_O3
                oneout=[1 0 0; 0 1 0; 0 0 1];
                oneout=repmat(oneout,1,18)';
                j=0;
                for i=1:3:nrruns*3
                    j=j+1;
                    test{i}   =find(run==j & oneout(:,1)==1); % Classify Temp with O1
                    test{i+1} =find(run==j & oneout(:,2)==1); % Classify Temp with O2
                    test{i+2} =find(run==j & oneout(:,3)==1); % Classify Temp with O3
                    
                    train{i}  =find(run~=j & oneout(:,1)~=1); % Train on Temp with O2 and O3
                    train{i+1}=find(run~=j & oneout(:,2)~=1); % Train on Temp with O1 and O3
                    train{i+2}=find(run~=j & oneout(:,3)~=1); % Train on Temp with O1 and O2
                end;
                
                %Read in P
                P={SPM.Vbeta(1:120).fname}';
                P(2:2:120)=[]; %Leave out betas for 1st derivative
                P(10:10:60)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
                
                V=spm_vol(char(P));
                Y=region_getdata(V,R);
                
                
                [D,M]=lmva_randomsubspace(Y',@tempord3_eventclass,'params',{cTemp,run,train,test},...
                    'num_iter',2000,'size_subspace',100);
                acc(:,k)=M.acc;
                
            end;
            varargout={acc};
            acc=[];
        end;
        
    case 'MVA_BG_random_subspaceOrder' %random subspace - left - right BG
        
        sn=varargin{1};
        nrregr=10;
        
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            reg={'maskbrainBGleft.nii','maskbrainBGright.nii'};
            for k=1:2
                disp(reg{k});
                R=region('roi_image',reg{k},1);
                R=region_calcregions(R);
                
                
                load SPM;
                nrruns=length(SPM.nscan);
                
                % Load SPM info file (to label beta values)
                D=dload('SPM_info.ana'); %% no derivatives
                c=D.condition(D.condition<nrregr)'; % extract conditions (leave out errors!)
                cOrd=mod((c-1),3)+1;
                run=D.RN(D.condition<nrregr); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
                
                mkdir(fullfile(glmDir,subj_name{s},'Ord_oneout')); % folder for each contrast
                
                % Generate column indices for cross-validation, where
                % cell i contains column indices of the respective test and
                % train set
                
                % Reminder:
                %%% T1_O1 T1_O2 T1_O3 %%% T2_O1 T2_O2 T2_O3 %%% T3_O1 T3_O2 T3_O3
                oneout=[0 0 1; 0 0 1; 0 0 1; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 1 0 0];
                oneout=repmat(oneout,6,1);
                j=0;
                for i=1:3:nrruns*3
                    j=j+1;
                    test{i}   =find(run==j & oneout(:,1)==1); % Classify Ord with T3
                    test{i+1} =find(run==j & oneout(:,2)==1); % Classify Ord with T2
                    test{i+2} =find(run==j & oneout(:,3)==1); % Classify Ord with T1
                    
                    train{i}  =find(run~=j & oneout(:,1)~=1); % Train on Ord with T1 and T2
                    train{i+1}=find(run~=j & oneout(:,2)~=1); % Train on Ord with T1 and T3
                    train{i+2}=find(run~=j & oneout(:,3)~=1); % Train on Ord with T2 and T3
                end;
                
                P={SPM.Vbeta(1:120).fname}';
                P(2:2:120)=[]; %Leave out betas for 1st derivative
                P(10:10:60)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
                
                
                V=spm_vol(char(P));
                Y=region_getdata(V,R);
                
                
                [D,M]=lmva_randomsubspace(Y',@tempord3_eventclass,'params',{cTemp,run,train,test},...
                    'num_iter',2000,'size_subspace',100);
                acc(:,k)=M.acc;
                
            end;
            varargout={acc};
            acc=[];
        end;
        
    case 'MVA_do' % performing the multivariate analysis on beta values
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<10)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<10); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:60).fname}';
            P(10:10:60)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            out={'mva160_acc_comb.nii','mva160_acc_TempComb.nii','mva160_acc_OrdComb.nii'};
            lmva_spm('volsearch160.mat',P,out,@tempord3_combinedclass,'params',{c,run,train,test});
        end;
        
    case 'MVA_do_d' % performing the multivariate analysis on beta values in the presence of 1st DERIVATIVE!!!
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<10)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<10); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:120).fname}';
            P(2:2:120)=[]; %Leave out betas for 1st derivative
            P(10:10:60)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            
            %%%New analysis: Mean patterns of T1/T2,T3 and O1,O2,O3
            %% Corrected (main patterns subtracted out):
            out={'mva160_acc_combCorrected.nii'};
            lmva_spm('volsearch160.mat',P,out,@tempord3_combinedclass_corrected4Main,'params',{c,run,train,test});
            %
            %            out={'mva160_acc_combCorrectedSUIT.nii'};
            %            lmva_spm('volsearch160SUIT.mat',P,out,@tempord3_combinedclass_corrected4Main,'params',{c,run,train,test});
            
            %             out={'mva160_acc_combCorrectedBG.nii'};
            %            lmva_spm('volsearch160BG.mat',P,out,@tempord3_combinedclass_corrected4Main,'params',{c,run,train,test});
            
            
        end;
        
    case 'MVA_do_d' % performing the multivariate analysis on beta values in the presence of 1st DERIVATIVE!!!
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<10)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<10); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:120).fname}';
            P(2:2:120)=[]; %Leave out betas for 1st derivative
            P(10:10:60)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            
            %%%New analysis: Mean patterns of T1/T2,T3 and O1,O2,O3
            %% Corrected (main patterns subtracted out):
            %             out={'mva160_acc_combCorrected.nii'};
            %             lmva_spm('volsearch160.mat',P,out,@tempord3_combinedclass_corrected4Main,'params',{c,run,train,test});
            %
            %            out={'mva160_acc_combCorrectedSUIT.nii'};
            %            lmva_spm('volsearch160SUIT.mat',P,out,@tempord3_combinedclass_corrected4Main,'params',{c,run,train,test});
            
            %             out={'mva160_acc_combCorrectedBG.nii'};
            %            lmva_spm('volsearch160BG.mat',P,out,@tempord3_combinedclass_corrected4Main,'params',{c,run,train,test});
            
            
            %% Uncorrected (main patterns not subtracted out):
            %             out={'mva160_acc_combUncorrected.nii'};
            %             lmva_spm('volsearch160.mat',P,out,@tempord3_combinedclass,'params',{c,run,train,test});
            %
            %            out={'mva160_acc_combUncorrectedSUIT.nii'};
            %            lmva_spm('volsearch160SUIT.mat',P,out,@tempord3_combinedclass,'params',{c,run,train,test});
            %
            %             out={'mva160_acc_combUncorrectedBG.nii'};
            %            lmva_spm('volsearch160BG.mat',P,out,@tempord3_combinedclass,'params',{c,run,train,test});
            
            %            out={'mva100_acc_combUncorrectedBG.nii'};
            %            lmva_spm('volsearch100BG.mat',P,out,@tempord3_combinedclass,'params',{c,run,train,test});
            %
            out={'mva160_acc_combUncorrectedBGright.nii'};
            lmva_spm('volsearch160BGright.mat',P,out,@tempord3_combinedclass,'params',{c,run,train,test});
            
            out={'mva160_acc_combUncorrectedBGleft.nii'};
            lmva_spm('volsearch160BGleft.mat',P,out,@tempord3_combinedclass,'params',{c,run,train,test});
            
            
            
            
        end;
        
    case 'MVA_do_dOrder' % performing the multivariate analysis on beta values in the presence of 1st DERIVATIVE!!!
        sn=varargin{1};
        nrregr=4; %nr of regressors including error regressor
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<nrregr)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<nrregr); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:48).fname}'; %nr of runs corresponds to
            P(2:2:48)=[]; %Leave out betas for 1st derivative
            P(4:4:24)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            out={'mva160_acc_Ord.nii'};
            lmva_spm('volsearch160.mat',P,out,@tempord3_eventclass,'params',{c,run,train,test});
            
        end;
        
    case 'MVA_do_dTemp' % performing the multivariate analysis on beta values in the presence of 1st DERIVATIVE!!!
        sn=varargin{1};
        nrregr=4; %nr of regressors including error regressor
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<nrregr)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<nrregr); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:48).fname}'; %nr of runs corresponds to
            P(2:2:48)=[]; %Leave out betas for 1st derivative
            P(4:4:24)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            out={'mva160_acc_Temp.nii'};
            lmva_spm('volsearch160.mat',P,out,@tempord3_eventclass,'params',{c,run,train,test});
            
        end;
        
    case 'MVA_do_dTemp_mean' % Mean across Betas, compare to dedicated regressor; performing the multivariate analysis on beta values in the presence of 1st DERIVATIVE!!!
        sn=varargin{1};
        nrregr=10; %nr of regressors including error regressor
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<nrregr)'; % extract conditions (leave out errors!)
            cTemp=floor((c-1)/3)+1;
            cTemp=cTemp(1:3:end); %only one regressor, because mean across O1,O2,O3
            run=D.RN(D.condition<nrregr); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            run=run(1:3:end);
            
            P={SPM.Vbeta(1:120).fname}'; Pall={SPM.Vbeta(1:120).fname}';
            P(2:2:120)=[];  %Leave out betas for 1st derivative
            P(10:10:60)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            j=0;
            
            mkdir(fullfile(glmDir,subj_name{s},'Temp')); % folder for each contrast
            for i=1:3:54
                j=j+1;
                TxO1=fullfile(glmDir, subj_name{s},P{i});
                TxO2=fullfile(glmDir, subj_name{s},P{i+1});
                TxO3=fullfile(glmDir, subj_name{s},P{i+2});
                TxMean=fullfile(glmDir, subj_name{s},['Temp/' Pall{j}]);
                
                spm_imcalc_ui({TxO1,TxO2,TxO3},TxMean,'(i1+i2+i3)/3',{});
                
            end;
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            copyfile('volsearch160.mat',[fullfile(glmDir, subj_name{s},'Temp')])
            cd(fullfile(glmDir, subj_name{s},'Temp'));
            
            P={SPM.Vbeta(1:18).fname}'
            
            out={'mva160_acc_TempMean.nii'};
            lmva_spm('volsearch160.mat',P,out,@tempord3_eventclass,'params',{cTemp,run,train,test});
            
        end;
        
    case 'MVA_do_TempOneout' % Neil Burgess suggestion!!!
        sn=varargin{1};
        nrregr=10; %nr of regressors including error regressor
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<nrregr)'; % extract conditions (leave out errors!)
            cTemp=floor((c-1)/3)+1;
            run=D.RN(D.condition<nrregr); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            mkdir(fullfile(glmDir,subj_name{s},'Temp_oneout')); % folder for each contrast
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            
            % Reminder:
            %%% T1_O1 T1_O2 T1_O3 %%% T2_O1 T2_O2 T2_O3 %%% T3_O1 T3_O2 T3_O3
            oneout=[1 0 0; 0 1 0; 0 0 1];
            oneout=repmat(oneout,1,18)';
            j=0;
            for i=1:3:nrruns*3
                j=j+1;
                test{i}   =find(run==j & oneout(:,1)==1); % Classify Temp with O1
                test{i+1} =find(run==j & oneout(:,2)==1); % Classify Temp with O2
                test{i+2} =find(run==j & oneout(:,3)==1); % Classify Temp with O3
                
                train{i}  =find(run~=j & oneout(:,1)~=1); % Train on Temp with O2 and O3
                train{i+1}=find(run~=j & oneout(:,2)~=1); % Train on Temp with O1 and O3
                train{i+2}=find(run~=j & oneout(:,3)~=1); % Train on Temp with O1 and O2
            end;
            %copyfile('volsearch160.mat',[fullfile(glmDir, subj_name{s},'Temp_oneout')]);
            %copyfile('beta_*',[fullfile(glmDir, subj_name{s},'Temp_oneout')]);
            %cd(fullfile(glmDir, subj_name{s},'Temp_oneout'));
            
            P={SPM.Vbeta(1:120).fname}';
            P(2:2:120)=[]; %Leave out betas for 1st derivative
            P(10:10:60)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            
            %             out={'Temp_oneout/mva160_acc_Temp_oneout.nii'};
            %             lmva_spm('volsearch160.mat',P,out,@tempord3_eventclass,'params',{cTemp,run,train,test});
            
            %               out={'Temp_oneout/mva160_acc_Temp_oneoutSUIT.nii'};
            %             lmva_spm('volsearch160SUIT.mat',P,out,@tempord3_eventclass,'params',{cTemp,run,train,test});
            
            %             out={'Temp_oneout/mva160_acc_Temp_oneoutBG.nii'};
            %             lmva_spm('volsearch160BG.mat',P,out,@tempord3_eventclass,'params',{cTemp,run,train,test});
            
            %             out={'Temp_oneout/mva100_acc_Temp_oneoutBG.nii'};
            %             lmva_spm('volsearch100BG.mat',P,out,@tempord3_eventclass,'params',{cTemp,run,train,test});
            %
            out={'Temp_oneout/mva160_acc_Temp_oneoutBGright.nii'};
            lmva_spm('volsearch160BGright.mat',P,out,@tempord3_eventclass,'params',{cTemp,run,train,test});
            out={'Temp_oneout/mva160_acc_Temp_oneoutBGleft.nii'};
            lmva_spm('volsearch160BGleft.mat',P,out,@tempord3_eventclass,'params',{cTemp,run,train,test});
            
        end;
        
    case 'MVA_do_OrdOneout' % Neil Burgess suggestion!!!
        sn=varargin{1};
        nrregr=10; %nr of regressors including error regressor
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<nrregr)'; % extract conditions (leave out errors!)
            cOrd=mod((c-1),3)+1;
            run=D.RN(D.condition<nrregr); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            mkdir(fullfile(glmDir,subj_name{s},'Ord_oneout')); % folder for each contrast
            
            % Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            
            % Reminder:
            %%% T1_O1 T1_O2 T1_O3 %%% T2_O1 T2_O2 T2_O3 %%% T3_O1 T3_O2 T3_O3
            oneout=[0 0 1; 0 0 1; 0 0 1; 0 1 0; 0 1 0; 0 1 0; 1 0 0; 1 0 0; 1 0 0];
            oneout=repmat(oneout,6,1);
            j=0;
            for i=1:3:nrruns*3
                j=j+1;
                test{i}   =find(run==j & oneout(:,1)==1); % Classify Ord with T3
                test{i+1} =find(run==j & oneout(:,2)==1); % Classify Ord with T2
                test{i+2} =find(run==j & oneout(:,3)==1); % Classify Ord with T1
                
                train{i}  =find(run~=j & oneout(:,1)~=1); % Train on Ord with T1 and T2
                train{i+1}=find(run~=j & oneout(:,2)~=1); % Train on Ord with T1 and T3
                train{i+2}=find(run~=j & oneout(:,3)~=1); % Train on Ord with T2 and T3
            end;
            
            P={SPM.Vbeta(1:120).fname}';
            P(2:2:120)=[]; %Leave out betas for 1st derivative
            P(10:10:60)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            
            %             out={'Ord_oneout/mva160_acc_Ord_oneout.nii'};
            %             lmva_spm('volsearch160.mat',P,out,@tempord3_eventclass,'params',{cOrd,run,train,test});
            
            %             out={'Ord_oneout/mva160_acc_Ord_oneoutSUIT.nii'};
            %             lmva_spm('volsearch160SUIT.mat',P,out,@tempord3_eventclass,'params',{cOrd,run,train,test});
            
            %              out={'Ord_oneout/mva160_acc_Ord_oneoutBG.nii'};
            %             lmva_spm('volsearch160BG.mat',P,out,@tempord3_eventclass,'params',{cOrd,run,train,test});
            %              out={'Ord_oneout/mva100_acc_Ord_oneoutBG.nii'};
            %             lmva_spm('volsearch100BG.mat',P,out,@tempord3_eventclass,'params',{cOrd,run,train,test});
            
            out={'Ord_oneout/mva160_acc_Ord_oneoutBGright.nii'};
            lmva_spm('volsearch160BGright.mat',P,out,@tempord3_eventclass,'params',{cOrd,run,train,test});
            out={'Ord_oneout/mva160_acc_Ord_oneoutBGleft.nii'};
            lmva_spm('volsearch160BGleft.mat',P,out,@tempord3_eventclass,'params',{cOrd,run,train,test});
            
            
        end;
        
    case 'MVA_do_dFinger' % performing the multivariate analysis on beta values in the presence of 1st DERIVATIVE!!!
        
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<6)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<6); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:72).fname}';
            P(2:2:72)=[]; %Leave out betas for 1st derivative
            P(6:6:36)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            out={'mva160_acc_finger.nii'};
            lmva_spm('volsearch160.mat',P,out,@tempord3_eventclass,'params',{c,run,train,test});
            
        end;
        
    case 'MVA_do_dTiming' % performing the multivariate analysis on beta values in the presence of 1st DERIVATIVE!!!
        
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<6)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<6); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:72).fname}';
            P(2:2:72)=[]; %Leave out betas for 1st derivative
            P(6:6:36)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            out={'mva160_acc_timing.nii'};
            lmva_spm('volsearch160.mat',P,out,@tempord3_eventclass,'params',{c,run,train,test});
            
        end;
        
    case 'MVA_do_dTiming_1stResponse' % performing the multivariate analysis on beta values in the presence of 1st DERIVATIVE!!!
        
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<4)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<4); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:48).fname}';
            P(2:2:48)=[]; %Leave out betas for 1st derivative
            P(4:4:24)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            out={'mva160_acc_timing.nii'};
            lmva_spm('volsearch160.mat',P,out,@tempord3_eventclass,'params',{c,run,train,test});
            
        end;
        
    case 'MVA_do_d_only' % performing the multivariate analysis on beta values for 1st DERIVATIVE!!!
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<10)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<10); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:120).fname}';
            P(1:2:120)=[]; %Leave out betas for basis function, keep 1st derivative only
            P(10:10:60)=[]; % Leave out 1st derivative betas estimating Error regressor! -> include only those that correspond to c!
            out={'mva100_acc_comb_d.nii','mva100_acc_TempComb_d.nii','mva100_acc_OrdComb_d.nii'};
            lmva_spm('volsearch100.mat',P,out,@tempord3_combinedclass,'params',{c,run,train,test});
            
        end;
        
    case 'MVA_do_d_both' % performing the multivariate analysis on beta values WITH 1st DERIVATIVE!!!
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            nrruns=length(SPM.nscan);
            
            % Load SPM info file (to label beta values)
            D=dload('SPM_info.ana'); %% no derivatives
            c=D.condition(D.condition<10)'; % extract conditions (leave out errors!)
            run=D.RN(D.condition<10); % extract run nr; or generate: run=kron([1:nrruns],ones(1,9)); %or just run=D.RN
            
            
            % OPTIONAL? Generate column indices for cross-validation, where
            % cell i contains column indices of the respective test and
            % train set
            for i=1:nrruns
                test{i}=find(run==i); % fprintf('test:'); display(test{i});
                train{i}=find(run~=i); % fprintf('train:'); display(train{i});
            end;
            
            
            P={SPM.Vbeta(1:120).fname}';
            P(19:20:120)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            P(19:19:114)=[]; % Leave out betas estimating Error regressor! -> include only those that correspond to c!
            out={'mva100_acc_comb_both.nii','mva100_acc_TempComb_both.nii','mva100_acc_OrdComb_both.nii'};
            lmva_spm_d('volsearch100.mat',P,out,@tempord3_combinedclass,'params',{c,run,train,test});
            
        end;
        
    case 'component_decomp'
        sn=varargin{1};
        for s=sn
            cd(fullfile(glmDir, subj_name{s}));
            load SPM;
            %%TODO!
            
        end;
        
    case 'glm_group_suit' %%% all subjects
        INfile={'zacc_Comb_160Corrected_surfvol','zacc_Comb_160Uncorrected_surfvol',...
            }
        for i=3:32
            INfile=   fullfile(suitDir);
        end;
        matlabbatch{1}.spm.stats.factorial_design.dir = {'/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/'};
        %%
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s03SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s04SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s05SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s06SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s07SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s08SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s10SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s12SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s13SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s14SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s15SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s16SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s17SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s18SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s19SUIT.nii,1'
            '/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/data_Temp_oneout/swcmmva160_zacc_Temp_oneout_s20SUIT.nii,1'
            };
        %%
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/Users/katja/Documents/MotorControlServer/data/tempord/tempord3/group_GLM_SUIT_trial_derivative/suit_mask4group.img,1'};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        
        spm_jobman('run',matlabbatch);
        
    case 'estimate_group_suit'
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {'E:\Data\MotorControlServer\project\tempord\tempord3\group_GLM_SUIT\SPM.mat'};
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        spm_jobman('run',matlabbatch);
        
    case 'contrast_group_suit'
        
    case 'MVA_overlap'  % Makes average across accuracies
        s=varargin{1};
        mask=fullfile(glmDir, subj_name{s},'mask.img');
        omask=fullfile(glmDir, subj_name{s},'combi_sum');
        for s=4:8
            P1=fullfile(glmDir,'anatomicals',subj_name{s},sprintf('c1%s_anatomical.nii',subj_name{s}));
            P2=fullfile(glmDir,'anatomicals',subj_name{s},sprintf('c2%s_anatomical.nii',subj_name{s}));
            sP1=fullfile(glmDir,'anatomicals',subj_name{s},sprintf('sc1%s_anatomical.nii',subj_name{s}));
            sP2=fullfile(glmDir,'anatomicals',subj_name{s},sprintf('sc2%s_anatomical.nii',subj_name{s}));
        end;
        spm_imcalc_ui({mask,sP1,sP2},omask,'i1 & (i2 +i3)>0.01',{});
        
    case 'MVA_zValue'
        
        s=varargin{1};
        cd(fullfile(glmDir,subj_name{s}));
        
        %%% Combi: 9 categories
        mu=1/9; %mu=0.11; not redone yet!
        N=6*9;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        images= {'_accuracy_Comb_160Corrected_surfvol',...
            '_accuracy_Comb_160Uncorrected_surfvol',...
            '_accuracy_Comb_160Corrected_surfvol_flip',...
            '_accuracy_Comb_160Uncorrected_surfvol_flip'};
        
        outimages={'_zacc_Comb_160Corrected_surfvol',...
            '_zacc_Comb_160Uncorrected_surfvol',...
            '_zacc_Comb_160Corrected_surfvol_flip',...
            '_zacc_Comb_160Uncorrected_surfvol_flip'};
        
        
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},[subj_name{s} images{j} '.nii']);
            output_image= fullfile(glmDir,subj_name{s},[subj_name{s} outimages{j} '.nii']);
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask!
        end;
        
        
        
        
        %%% TempOneOut and OrdOneOut
        mu=1/3; %previously: mu=0.33; not redone yet!
        N=6*3*3; %%% 3 times per Previously erroneously:N=6*3;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        images= {'_accuracy_Comb_160OrdOneOut_surfvol',...
            '_accuracy_Comb_160TempOneOut_surfvol',...
            '_accuracy_Comb_160OrdOneOut_surfvol_flip',...
            '_accuracy_Comb_160TempOneOut_surfvol_flip'};
        
        outimages={'_zacc_Comb_160OrdOneOut_surfvol',...
            '_zacc_Comb_160TempOneOut_surfvol',...
            '_zacc_Comb_160OrdOneOut_surfvol_flip',...
            '_zacc_Comb_160TempOneOut_surfvol_flip'};
        
        
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},[subj_name{s} images{j} '.nii']);
            output_image= fullfile(glmDir,subj_name{s},[subj_name{s} outimages{j} '.nii']);
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask!
        end;
        
    case 'MVA_zValue_combined'
        s=varargin(1);
        %%% Combi: 9 categories
        mu=1/9; %mu=0.11; not redone yet!
        N=6*9;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        cd(fullfile(glmDir,subj_name{s}))
        
        %%%
        %% Corrected:
        %         images= {'mva160_acc_combCorrected.nii'};
        %         outimages={'mva160_zacc_combCorrected.nii'};
        
        %         images= {'mva160_acc_combCorrectedSUIT.nii'};
        %         outimages={'mva160_zacc_combCorrectedSUIT.nii'};
        
        %         images= {'mva160_acc_combCorrectedBG.nii'};
        %         outimages={'mva160_zacc_combCorrectedBG.nii'};
        
        
        %% Uncorrected:
        %         images= {'mva160_acc_combUncorrected.nii'};
        %         outimages={'mva160_zacc_combUncorrected.nii'};
        
        %         images= {'mva160_acc_combUncorrectedSUIT.nii'};
        %         outimages={'mva160_zacc_combUncorrectedSUIT.nii'};
        
        %         images= {'mva160_acc_combUncorrectedBG.nii'};
        %         outimages={'mva160_zacc_combUncorrectedBG.nii'};
        
        %         images= {'mva100_acc_combUncorrectedBG.nii'};
        %         outimages={'mva100_zacc_combUncorrectedBG.nii'};
        
        %         images= {'mva160_acc_combUncorrectedBGright.nii'};
        %         outimages={'mva160_zacc_combUncorrectedBGright.nii'};
        
        %          images= {'mva160_acc_combUncorrectedBGleft.nii'};
        %         outimages={'mva160_zacc_combUncorrectedBGleft.nii'};
        
        
        %         s=varargin{1};
        %         for j=1:numel(images)
        %             input_image= fullfile(glmDir,subj_name{s},images{j});
        %             output_image= fullfile(glmDir,subj_name{s},outimages{j});
        %             spmj_imcalc_mtx(input_image, output_image,...
        %             sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask!
        %         end;
        
    case 'MVA_zValue_single'
        %%% Tempo and Order: 3 categories
        s=varargin{1};
        
        mu=1/3;%previously: mu=0.33; not redone yet!
        N=6*9; %%% Previously erroneously: N=6*3;!!!
        sigma=sqrt(mu*(1-mu)*1/N);
        %         images= {'mva160_acc_Ord.nii'}
        %         outimages= {'mva160_zacc_Ord.nii'}
        images= {'mva160_acc_TempComb.nii','mva160_acc_OrdComb.nii'}
        outimages= {'mva160_zacc_TempComb.nii','mva160_zacc_OrdComb.nii'}
        
        s=varargin{1};
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},images{j});
            output_image= fullfile(glmDir,subj_name{s},outimages{j});
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask!
        end;
        
    case 'MVA_zValue_TempOneout'
        %%% Tempo and Order: 3 categories
        s=varargin{1};
        
        mu=1/3; %previously: mu=0.33; not redone yet!
        N=6*3*3; %%% 3 times per Previously erroneously:N=6*3;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        %         images= {'mva160_acc_Temp_oneout.nii'}
        %         outimages= {'mva160_zacc_Temp_oneout.nii'}
        
        %         images= {'mva160_acc_Temp_oneoutSUIT.nii'}
        %         outimages= {'mva160_zacc_Temp_oneoutSUIT.nii'}
        
        %          images= {'mva160_acc_Temp_oneoutBG.nii'}
        %         outimages= {'mva160_zacc_Temp_oneoutBG.nii'}
        
        %           images= {'mva100_acc_Temp_oneoutBG.nii'}
        %         outimages= {'mva100_zacc_Temp_oneoutBG.nii'}
        
        %           images= {'mva160_acc_Temp_oneoutBGright.nii'}
        %         outimages= {'mva160_zacc_Temp_oneoutBGright.nii'}
        
        
        images= {'mva160_acc_Temp_oneoutBGleft.nii'}
        outimages= {'mva160_zacc_Temp_oneoutBGleft.nii'}
        
        s=varargin{1};
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},'Temp_oneout',images{j});
            output_image= fullfile(glmDir,subj_name{s},'Temp_oneout',outimages{j});
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask!
        end;
        
    case 'MVA_zValue_OrdOneout'
        %%% Tempo and Order: 3 categories
        s=varargin{1};
        
        mu=1/3; %previously: mu=0.33; not redone yet!
        N=6*3*3; %%% 3 times per Previously erroneously:N=6*3;
        sigma=sqrt(mu*(1-mu)*1/N);
        
        %         images= {'mva160_acc_Ord_oneout.nii'}
        %         outimages= {'mva160_zacc_Ord_oneout.nii'}
        
        %         images= {'mva160_acc_Ord_oneoutSUIT.nii'}
        %         outimages= {'mva160_zacc_Ord_oneoutSUIT.nii'}
        
        %          images= {'mva160_acc_Ord_oneoutBG.nii'}
        %         outimages= {'mva160_zacc_Ord_oneoutBG.nii'}
        
        %           images= {'mva100_acc_Ord_oneoutBG.nii'}
        %         outimages= {'mva100_zacc_Ord_oneoutBG.nii'}
        
        %           images= {'mva160_acc_Ord_oneoutBGright.nii'}
        %         outimages= {'mva160_zacc_Ord_oneoutBGright.nii'}
        
        images= {'mva160_acc_Ord_oneoutBGleft.nii'}
        outimages= {'mva160_zacc_Ord_oneoutBGleft.nii'}
        
        s=varargin{1};
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},'Ord_oneout',images{j});
            output_image= fullfile(glmDir,subj_name{s},'Ord_oneout',outimages{j});
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask!
        end;
        
    case 'MVA_zValue_Timing'
        %%% Tempo and Order: 5 temporal intervals
        mu=0.2;
        N=6*5;
        sigma=sqrt(mu*(1-mu)*1/N);
        images= {'mva160_acc_timing.nii'}
        outimages= {'mva160_zacc_timing.nii'}
        s=varargin{1};
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},images{j});
            output_image= fullfile(glmDir,subj_name{s},outimages{j});
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask!
        end;
        
    case 'MVA_zValue_Finger'
        %%% Tempo and Order: 5 fingers
        mu=0.2;
        N=6*5;
        sigma=sqrt(mu*(1-mu)*1/N);
        images= {'mva160_acc_finger.nii'}
        outimages= {'mva160_zacc_finger.nii'}
        s=varargin{1};
        for j=1:numel(images)
            input_image= fullfile(glmDir,subj_name{s},images{j});
            output_image= fullfile(glmDir,subj_name{s},outimages{j});
            spmj_imcalc_mtx(input_image, output_image,...
                sprintf('(X./X).*((X-%d)/%d)',mu, sigma)); %(X./X) acts like a mask!
        end;
        
    case 'MVA_smooth' %%%Smoothing before suit reslice %%Subject space
        
        s=varargin{1};
        
        %Comb:
        %         comb=fullfile(glmDir, subj_name{s},['mva160_zacc_combCorrected.nii']);
        %         scomb=fullfile(glmDir, subj_name{s},['smva160_zacc_combCorrected.nii']);
        %
        %         %Temp:
        %         temp=fullfile(glmDir, subj_name{s},'Temp_oneout',['mva160_zacc_Temp_oneout.nii']);
        %         stemp=fullfile(glmDir, subj_name{s},'Temp_oneout',['smva160_zacc_Temp_oneout.nii']);
        %
        %         %Ord:
        %         ord=fullfile(glmDir, subj_name{s},'Ord_oneout',['mva160_zacc_Ord_oneout.nii']);
        %         sord=fullfile(glmDir, subj_name{s},'Ord_oneout',['smva160_zacc_Ord_oneout.nii']);
        
        
        
        %% Comb Corrected
        %         comb=fullfile(glmDir, subj_name{s},['mva160_zacc_combCorrectedBG.nii']);
        %         scomb=fullfile(glmDir, subj_name{s},['smva160_zacc_combCorrectedBG.nii']);
        
        %% Comb uncorrected
        %          comb=fullfile(glmDir, subj_name{s},['mva160_zacc_combUncorrectedBG.nii']);
        %         scomb=fullfile(glmDir, subj_name{s},['smva160_zacc_combUncorrectedBG.nii']);
        %
        %         %Temp:
        %         temp=fullfile(glmDir, subj_name{s},'Temp_oneout',['mva160_zacc_Temp_oneoutBG.nii']);
        %         stemp=fullfile(glmDir, subj_name{s},'Temp_oneout',['smva160_zacc_Temp_oneoutBG.nii']);
        %
        %         %Ord:
        %         ord=fullfile(glmDir, subj_name{s},'Ord_oneout',['mva160_zacc_Ord_oneoutBG.nii']);
        %         sord=fullfile(glmDir, subj_name{s},'Ord_oneout',['smva160_zacc_Ord_oneoutBG.nii']);
        
        
        comb=fullfile(glmDir, subj_name{s},['mva160_zacc_combUncorrectedBGleft.nii']);
        scomb=fullfile(glmDir, subj_name{s},['smva160_zacc_combUncorrectedBGleft.nii']);
        
        %Temp:
        temp=fullfile(glmDir, subj_name{s},'Temp_oneout',['mva160_zacc_Temp_oneoutBGleft.nii']);
        stemp=fullfile(glmDir, subj_name{s},'Temp_oneout',['smva160_zacc_Temp_oneoutBGleft.nii']);
        
        %Ord:
        ord=fullfile(glmDir, subj_name{s},'Ord_oneout',['mva160_zacc_Ord_oneoutBGleft.nii']);
        sord=fullfile(glmDir, subj_name{s},'Ord_oneout',['smva160_zacc_Ord_oneoutBGleft.nii']);
        
        spm_smooth(comb,scomb,[4 4 4]);
        spm_smooth(temp,stemp,[4 4 4]);
        spm_smooth(ord,sord,[4 4 4]);
        
    case 'MVA_smooth_Timing' %%%Smoothing before suit reslice %%Subject space
        
        s=varargin{1};
        
        temp=fullfile(glmDir, subj_name{s},['mva160_zacc_timing.nii']);
        
        
        stemp=fullfile(glmDir, subj_name{s},['smva160_zacc_timing.nii']);
        
        spm_smooth(temp,stemp,[4 4 4]);
        
    case 'MVA_smooth_Finger' %%%Smoothing before suit reslice %%Subject space
        
        s=varargin{1};
        
        ord=fullfile(glmDir, subj_name{s},['mva160_zacc_finger.nii']);
        
        
        sord=fullfile(glmDir, subj_name{s},['smva160_zacc_finger.nii']);
        
        spm_smooth(ord,sord,[4 4 4]);
        
    case 'MNI_normalization_write'
        %        images= {'smva160_zacc_combCorrected.nii','Temp_oneout/smva160_zacc_Temp_oneout.nii','Ord_oneout/smva160_zacc_Ord_oneout.nii'};
        %         images= {'smva160_zacc_combUncorrectedBG.nii','smva160_zacc_combCorrectedBG.nii','Temp_oneout/smva160_zacc_Temp_oneoutBG.nii','Ord_oneout/smva160_zacc_Ord_oneoutBG.nii','maskbrainBG.nii'};
        %         outDir= {'data_Comb_Uncorrected','data_Comb_corrected','data_Temp_oneout','data_Ord_oneout','data_BG'}
        images= {'smva100_zacc_combUncorrectedBG.nii','Temp_oneout/smva100_zacc_Temp_oneoutBG.nii','Ord_oneout/smva100_zacc_Ord_oneoutBG.nii','maskbrainBG.nii'};
        outDir= {'data_Comb_Uncorrected','data_Temp_oneout','data_Ord_oneout','data_BG'}
        
        %%LEFT RIGHT BG SEPARATELY:
        
        %         images= {'smva160_zacc_combUncorrectedBGright.nii','Temp_oneout/smva160_zacc_Temp_oneoutBGright.nii','Ord_oneout/smva160_zacc_Ord_oneoutBGright.nii','maskbrainBGright.nii'};
        %         outDir= {'data_Comb_Uncorrected','data_Temp_oneout','data_Ord_oneout','data_BG'}
        
        images= {'smva160_zacc_combUncorrectedBGleft.nii','Temp_oneout/smva160_zacc_Temp_oneoutBGleft.nii','Ord_oneout/smva160_zacc_Ord_oneoutBGleft.nii','maskbrainBGleft.nii'};
        outDir= {'data_Comb_Uncorrected','data_Temp_oneout','data_Ord_oneout','data_BG'}
        
        %
        
        s=varargin{1};
        defor= fullfile(baseDir, 'anatomicals', subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
        for j=1:numel(images)
            [dir,name,ext]=spm_fileparts(images{j});
            sn_images{j}= fullfile(glmDir,subj_name{s}, images{j});
            out_images{j}= fullfile(groupDir,outDir{j},[name '_' subj_name{s} '.nii']);
        end
        spmj_normalization_write(defor, sn_images,'outimages',out_images);
        
    case 'MNI_mask_anatomical' %%%%__________________________________________________
        %sl4_imana('MNI_mask_anatomical', 3:8)
        sn=varargin{1};
        cd (fullfile(groupDir,'data'));
        i=0;
        for s=sn
            i=i+1;
            images{i}= fullfile(glmDir,subj_name{s},['maskbrain.nii']);
        end
        spmj_imcalc_mtx(images,'mask_avrg.nii','nanmean(X)'); %Mean across subject
        spmj_imcalc_mtx('mask_avrg.nii','mask_thres.nii','X>0.4'); %Threshold mask
        
    case 'MNI_BGmask'   %normalize BG mask for group mask
        sn=varargin{1};
        cd (fullfile(groupDir,'data_BG'));
        i=0;
        for s=sn
            i=i+1;
            images{i}= fullfile(groupDir,'data_BG',['maskbrainBGleft_', subj_name{s}, '.nii']);
        end
        spmj_imcalc_mtx(images,'mask_avrgleft.nii','nanmean(X)'); %Mean across subject
        spmj_imcalc_mtx('mask_avrgleft.nii','mask_thresleft.nii','X>0.7'); %Threshold mask
        
    case 'SUIT_mask'   %normalize BG mask for group mask
        %         sn=varargin{1};
        %         cd (fullfile(suitDir));
        %         i=0;
        %         for s=sn
        %             i=i+1;
        %             images{i}= fullfile(glmDir,subj_name{s},['wcmmva160_numvoxSUIT.nii']);
        %         end
        %         spmj_imcalc_mtx(images,'mask_avrg.nii','nanmean(X)'); %Mean across subject
        %         spmj_imcalc_mtx('mask_avrg.nii','mask_thres.nii','X>0.7'); %Threshold mask
        
        %Mask only grey matter:
        suitGrey=fullfile(anatDir,['suit_templates/Cerebellum-SUIT-maxprob.nii']);
        %         in=['mask_thres.nii'];
        
        out=['mask_threshGrey.nii'];
        %         spm_imcalc_ui({suitGrey,in},out,'(i1>0.1).*i2',{}); %mask suit grey matter prob above 0.1 with group mask
        spm_imcalc_ui(suitGrey,out,'(i1>0.1)',{}); %mask suit grey matter prob above 0.1 with group mask
            
    case 'ttest_suit'
        
        dirIn={'Comb','Temp','Ord','CombCorrected'};
        accIn={'_zacc_Comb_160Uncorrected_surfvol','_zacc_Comb_160TempOneOut_surfvol',...
            '_zacc_Comb_160OrdOneOut_surfvol','_zacc_Comb_160Corrected_surfvol'};
        for i=1:length(dirIn)
            scansIn=[];
            for s=3:34
                if s<19
                    scans={fullfile(suitDir,dirIn{i},['swcm' subj_name{s} accIn{i} '.nii'])};
                    scansIn=[scansIn;scans];
                else
                    scans={fullfile(suitDir,dirIn{i},['swcm' subj_name{s} accIn{i} 'flip.nii'])};
                    scansIn=[scansIn;scans];
                end;
            end;
            
            matlabbatch{1}.spm.stats.factorial_design.dir = fullfile(suitDir,dirIn{i});
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scansIn;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = fullfile(suitDir, ['mask_threshGrey.nii,1']);
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            spm_jobman('run',matlabbatch);
            
            %             matlabbatch{1}.spm.stats.fmri_est.spmmat = fullfile(suitDir,dirIn{i},['SPM.mat']);
            %             matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            %             spm_jobman('run',matlabbatch);
        end;
        
    case 'MNI_normalization_write_Timing'
        mkdir(fullfile(groupDir,'data')); % folder for each contrast
        images= {'smva160_zacc_timing.nii'};
        s=varargin{1};
        defor= fullfile(baseDir, 'anatomicals', subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
        for j=1:numel(images)
            [dir,name,ext]=spm_fileparts(images{j});
            sn_images{j}= fullfile(glmDir,subj_name{s},images{j});
            out_images{j}= fullfile(groupDir,'data',[name '_' subj_name{s} '.nii']);
        end
        spmj_normalization_write(defor, sn_images,'outimages',out_images);
        
    case 'MNI_normalization_write_Finger'
        mkdir(fullfile(groupDir,'data')); % folder for each contrast
        images= {'smva160_zacc_finger.nii'};
        s=varargin{1};
        defor= fullfile(baseDir, 'anatomicals', subj_name{s}, [subj_name{s}, '_anatomical_seg_sn.mat']);
        for j=1:numel(images)
            [dir,name,ext]=spm_fileparts(images{j});
            sn_images{j}= fullfile(glmDir,subj_name{s},images{j});
            out_images{j}= fullfile(groupDir,'data',[name '_' subj_name{s} '.nii']);
        end
        spmj_normalization_write(defor, sn_images,'outimages',out_images);
        
    case 'MNI_mask'
        spmj_imcalc_mtx([],'mask_avrg.nii','nansum(X)/10');
        spmj_imcalc_mtx('mask_avrg.nii','mask_thres.nii','X>0.4');
        
        
        
        
        %%%%%%% ROI analysis
        %%% Volume based:
    
    case 'make_ROI_SUIT'
        %%Group:
        %Make conjunction between suit ROIs and Group acc data:
        %Make suit ROIs with imagecalc first, e.g. Left_CrusI -> i1==8 on
        %Cerebellum_SUIT.nii
        
        
        cd(suitDir);
        accT={'Con/spmT_0001.img','data_Comb_Uncorrected/spmT_0001.img','data_Ord_oneout/spmT_0001.img','data_Temp_oneout/spmT_0001.img','data_Comb_corrected/spmT_0001.img'}; %group
        acc={'con','comb_uncorrected','ord','temp','comb'};
        %          roiName={'Right_VI'};
        %            roiName={'Left_VI'}; %%%sign.
        %            roiName={'RH'};
        roiName={'LH'};
        %          roiName={'Right_VI_CrusI'};
        %          roiName={'Right_CrusI'};
        %         roiName={'Vermis_VI'};
        %          roiName={'Left_CrusI'};
        
        for j=1:length(roiName)
            for i=1:length(accT)
                suitROI=fullfile(anatDir,'suit_templates',[roiName{j},'.img']);
                out=fullfile(regDir,['ROI_', roiName{j}, acc{i},'.nii']);
                %                  spm_imcalc_ui({accT{i},suitROI},out,'(i1>2.6).*i2',{}); %mask group acc with suit ROI
                spm_imcalc_ui({accT{i},suitROI},out,'(i1>2.6)',{}); %mask group acc with suit ROI
            end;
        end;
        
        %     case 'acc_ROI'
        %         %extract location for regions:
        %         %roiName={'Right_VI'};
        %         suitROI=fullfile(anatDir,'suit_templates',['Right_VI_CrusI.img']);
        %
        %         accDir={'data_Comb_Uncorrected'};
        %         accName={'Comb_Uncorrected'};
        %
        %
        %         R=region('roi_image',suitROI,1);
        %         R=region_calcregions(R);
        %
        %
        %         %Func z-trasformed data (in SUIT space):
        %         for i=3:18
        %         zAccSubj{i-2}=fullfile(suitDir,accDir,['swcmmva160_zacc_Comb_Uncorrected_',subj_name{i}, 'SUIT.nii']);
        %         end;
        %         V=spm_vol(zAccSubj);
        %         D=region_getdata(V,R);
        
    case 'acc_ROI_SUIT'
        %%Individual:
        
        %subj_name{4:10};
        cd(regDir);
        %%%%%%%%%%%%%%%%%%%%
        %%Func(acc) ROI
        %         %loop over k:
        %         roiName={'Right_VI','Vermis_VI','Left_CrusI'}; %group peak acc on cerebellum-suit
        %         accT={'data_Ord_oneout/spmT_0001.img','data_Temp_oneout/spmT_0001.img','data_Comb_corrected/spmT_0001.img'}; %group
        %
        %         %loop over j:
        %         accDir={'data_Comb_corrected','data_Ord_oneout','data_Temp_oneout'};
        %         accName={'Comb_oneout','Ord_oneout','Temp_oneout'};
        %        acc={'wcmmva160_zacc_combCorrectedSUIT.nii','Ord_oneout/wcmmva160_zacc_Ord_oneoutSUIT.nii','Temp_oneout/wcmmva160_zacc_Temp_oneoutSUIT.nii'};
        %
        %
        %         R=[];
        %         for j=1:length(accDir) % all acc types
        %             for i=3:length(subj_name) % all subjects (3:11)
        %                 %zAccSubj=fullfile(suitDir,accDir{j},['swcmmva160_zacc_',accName{j},'_',subj_name{i}, '.nii']);
        %                 zAccSubj=fullfile(glmDir,subj_name{i},acc{j}); %non_smoothed
        %                 for k=1:length(roiName) % all ROIs
        %                     out=fullfile(['swcmmva160_zacc_', accName{j},'_', subj_name{i},'_', roiName{k}, '.nii']);
        %                     accROI=fullfile(suitDir,accT{k});
        %                     spm_imcalc_ui({zAccSubj,accROI},out,'(i2>2.23).*i1',{});
        %
        %                     %write out ROI accuracies in a file %highest 10% zacc
        %                     V=spm_vol(out);
        %                     X=spm_read_vols(V);
        %                     X=X(X>0);
        %                     X=sort(X,'descend');
        %                     lengthX=length(X);
        %                     highestX= X(1:round(lengthX*0.2),1); %highest 10% zacc
        %                     highestXmean=median(highestX); %mean
        %                     D.subj=i;
        %                     D.roi=k;
        %                     D.accName=j;
        %                     D.mean=highestXmean;
        %                     R=addstruct(R,D,'row','force');
        %                     clear D;
        %                 end;
        %             end;
        %         end;
        %         dsave('SUIT_roi.txt',R);
        
        
        
        %%%%%%%%%%%%%%%%%%%%
        %%anat ROI
        %         %loop over k:
        %         roiName={'Right_VI','Vermis_VI','Left_CrusI'}; %as defined on cerebellum-suit
        %         %loop over j:
        %         accDir={'data_Comb_corrected','data_Ord_oneout','data_Temp_oneout'};
        %         accName={'Comb_oneout','Ord_oneout','Temp_oneout'};
        %         acc={'wcmmva160_zacc_combCorrectedSUIT.nii','Ord_oneout/wcmmva160_zacc_Ord_oneoutSUIT.nii','Temp_oneout/wcmmva160_zacc_Temp_oneoutSUIT.nii'};
        %
        %
        %         R=[];
        %         for j=1:length(accDir) % all acc types
        %             for i=3:length(subj_name) % all subjects (3:11)
        %                 %zAccSubj=fullfile(suitDir,accDir{j},['swcmmva160_zacc_',accName{j},'_',subj_name{i}, '.nii']);
        %                  zAccSubj=fullfile(glmDir,subj_name{i},acc{j}); %non_smoothed
        %                 for k=1:length(roiName) % all ROIs
        %                     out=fullfile(['swcmmva160_zacc_', accName{j},'_', subj_name{i},'_', roiName{k}, '.nii']);
        %                     suitROI=fullfile(anatDir,'suit_templates',[roiName{k},'.img']);
        %                     spm_imcalc_ui({zAccSubj,suitROI},out,'(i2>0).*i1',{});
        %
        %                     %write out ROI accuracies in a file %highest 10% zacc
        %                     V=spm_vol(out);
        %                     X=spm_read_vols(V);
        %                     X=X(X>0);
        %                     X=sort(X,'descend');
        %                     lengthX=length(X);
        %                     highestX= X(1:round(lengthX*0.2),1); %highest 10% zacc
        %                     highestXmean=mean(highestX); %mean
        %                     D.subj=i;
        %                     D.roi=k;
        %                     D.accName=j;
        %                     D.mean=highestXmean;
        %                     R=addstruct(R,D,'row','force');
        %                     clear D;
        %                 end;
        %             end;
        %         end;
        %         dsave('SUIT_roi.txt',R);
        
        %%%%%%%%%%%%%%%%%%%%
        %%Anat + Func(acc) ROI
        
        %           %loop over k:
        %         roiMask={'ROI_Right_VIord','ROI_Vermis_VItemp','ROI_Left_CrusIcomb'}; %as defined on group
        %
        %         %loop over j:
        %         accDir={'data_Comb_corrected','data_Ord_oneout','data_Temp_oneout'};
        %         accName={'Comb_oneout','Ord_oneout','Temp_oneout'};
        %         acc={'wcmmva160_zacc_combCorrected.nii','Ord_oneout/wcmmva160_zacc_Ord_oneout.nii','Temp_oneout/wcmmva160_zacc_Temp_oneout.nii'};
        %         R=[];
        %         for j=1:length(accDir) % all acc types
        %             for i=3:length(subj_name) % all subjects (3:11)
        %                 %zAccSubj=fullfile(suitDir,accDir{j},['swcmmva160_zacc_',accName{j},'_',subj_name{i}, '.nii']);
        %                 zAccSubj=fullfile(glmDir,subj_name{i},acc{j}); %non_smoothed
        %
        %                 for k=1:length(roiMask) % all ROIs
        %                     out=fullfile(['swcmmva160_zacc_', accName{j},'_', subj_name{i},'_', roiMask{k}, '.nii']);
        %                     mask=fullfile([roiMask{k}, '.nii']);
        %                     spm_imcalc_ui({zAccSubj,mask},out,'(i2>0).*i1',{});
        %
        %                     %write out ROI accuracies in a file
        %                     %write out ROI accuracies in a file %highest 10% zacc
        %                     V=spm_vol(out);
        %                     X=spm_read_vols(V);
        %                     X=X(X>0);
        %                     X=sort(X,'descend');
        %                     lengthX=length(X);
        %                     highestX= X(1:round(lengthX*0.05),1); %highest 10% zacc
        %                     highestXmean=median(highestX); %mean
        %                     D.subj=i;
        %                     D.roi=k;
        %                     D.accName=j;
        %                     D.mean=highestXmean;
        %                     R=addstruct(R,D,'row','force');
        %                     clear D;
        %                 end;
        %             end;
        %         end;
        %         dsave('SUIT_roi.txt',R);
        
        %loop over k:
        %            roiMask={'ROI_Right_VIcon','ROI_Right_VIcomb_uncorrected','ROI_Right_VIord','ROI_Right_VItemp'}; %as defined on group
        %            fileout={'ROI_Right_VIcon.txt','ROI_Right_VIcomb_uncorrected.txt','ROI_Right_VIord.txt','ROI_Right_VItemp.txt'};
        
        %            roiMask={'ROI_Left_VIcon','ROI_Left_VIcomb_uncorrected','ROI_Left_VIord','ROI_Left_VItemp'}; %as defined on group
        %            fileout={'ROI_Left_VIcon.txt','ROI_Left_VIcomb_uncorrected.txt','ROI_Left_VIord.txt','ROI_Left_VItemp.txt'};
        
        %            roiMask={'ROI_RHcon','ROI_RHcomb_uncorrected','ROI_RHord','ROI_RHtemp'}; %as defined on group
        %            fileout={'ROI_RHcon.txt','ROI_RHcomb_uncorrected.txt','ROI_RHord.txt','ROI_RHtemp.txt'};
        
        roiMask={'ROI_LHcon','ROI_LHcomb_uncorrected','ROI_LHord','ROI_LHtemp'}; %as defined on group
        fileout={'ROI_LHcon.txt','ROI_LHcomb_uncorrected.txt','ROI_LHord.txt','ROI_LHtemp.txt'};
        
        
        %            roiMask={'ROI_Vermis_VIcon','ROI_Vermis_VIcomb_uncorrected','ROI_Vermis_VIord','ROI_Vermis_VItemp'}; %as defined on group
        %            fileout={'ROI_Vermis_VIcon.txt','ROI_Vermis_VIcomb_uncorrected.txt','ROI_Vermis_VIord.txt','ROI_Vermis_VItemp.txt'};
        %
        %             roiMask={'ROI_Left_CrusIcon','ROI_Left_CrusIcomb_uncorrected','ROI_Left_CrusIord','ROI_Left_CrusItemp'}; %as defined on group
        %            fileout={'ROI_Left_CrusIcon.txt','ROI_Left_CrusIcomb_uncorrected.txt','ROI_Left_CrusIord.txt','ROI_Left_CrusItemp.txt'};
        
        
        %             roiMask={'ROI_Right_VI_CrusIcon','ROI_Right_VI_CrusIcomb_uncorrected','ROI_Right_VI_CrusIord','ROI_Right_VI_CrusItemp'}; %as defined on group
        %            fileout={'ROI_Right_VI_CrusIcon.txt','ROI_Right_VI_CrusIcomb_uncorrected.txt','ROI_Right_VI_CrusIord.txt','ROI_Right_VI_CrusItemp.txt'};
        
        %             roiMask={'ROI_Right_CrusIcon','ROI_Right_CrusIcomb_uncorrected','ROI_Right_CrusIord','ROI_Right_CrusItemp'}; %as defined on group
        %            fileout={'ROI_Right_CrusIcon.txt','ROI_Right_CrusIcomb_uncorrected.txt','ROI_Right_CrusIord.txt','ROI_Right_CrusItemp.txt'};
        %
        
        
        %           roiMask={'ROI_Right_VIord'}; %as defined on group
        %          roiMask={'ROI_Right_VIcon'}; %as defined on group
        %         roiMask={'ROI_Left_CrusIcomb_uncorrected'};
        %         roiMask={'ROI_Vermis_VIcomb_uncorrected'};
        
        %loop over j:
        accDir={'Con','data_Comb_Uncorrected','data_Comb_corrected','data_Ord_oneout','data_Temp_oneout'};
        accName={'Con','Comb_Uncorrected','Comb_oneout','Ord_oneout','Temp_oneout'};
        acc={'wcmcon_0001.img','wcmmva160_zacc_combUncorrectedSUIT.nii','wcmmva160_zacc_combCorrectedSUIT.nii','Ord_oneout/wcmmva160_zacc_Ord_oneoutSUIT.nii','Temp_oneout/wcmmva160_zacc_Temp_oneoutSUIT.nii'};
        R=[];
        for k=1:length(roiMask)
            %         for k=1
            %         for j=1:length(accDir) % all acc types
            for j=1:length(accName) % all acc types
                for i=3:length(subj_name) % all subjects (3:11)
                    %zAccSubj=fullfile(suitDir,accDir{j},['swcmmva160_zacc_',accName{j},'_',subj_name{i}, '.nii']);
                    zAccSubj=fullfile(glmDir,subj_name{i},acc{j}); %non_smoothed
                    
                    
                    %                 %%% with region Toolbox:
                    %                 ROI=region('roi_image',roiMask,1); %Read in ROI defined as 1
                    %                 ROI=region_calcregions(ROI);
                    %
                    %
                    %                 for k=1:length(roiMask) % all ROIs
                    %                     out=fullfile(['swcmmva160_zacc_', accName{j},'_', subj_name{i},'_', roiMask{k}, '.nii']);
                    %                     %write out ROI accuracies in a file
                    %                     %write out ROI accuracies in a file %highest 10% zacc
                    %                     V=spm_vol(out);
                    %                     X=region_getdata(V,ROI);
                    %                     X=X(X>0);
                    %                     X=sort(X,'descend');
                    %                     D.subj=i;
                    %                     D.roi=k;
                    %                     D.accName=j;
                    %                     D.mean=mean(X);
                    %                     R=addstruct(R,D,'row','force');
                    %                     clear D;
                    %                 end;
                    %             end;
                    %         end;
                    
                    
                    out=fullfile(['wcmmva160_zacc_', accName{j},'_', subj_name{i},'_', roiMask{k}, '.nii']);
                    mask=fullfile(regDir, [roiMask{k}, '.nii']);
                    spm_imcalc_ui({zAccSubj,mask},out,'(i2>0).*i1',{});
                    
                    %write out ROI accuracies in a file
                    %write out ROI accuracies in a file %highest 10% zacc
                    V=spm_vol(out);
                    X=spm_read_vols(V);
                    X=X(X~=0);
                    
                    X=sort(X,'descend');
                    D.subj=i;
                    D.roi=k;
                    D.accName=j;
                    D.mean=mean(X);
                    R=addstruct(R,D,'row','force');
                    clear D;
                    
                end;
            end;
            dsave(fileout{k},R);
        end;
        
    case 'ROI_SUIT_stats'
        
        cd(regDir);
        %                   fileIn={'ROI_Left_VIcon.txt','ROI_Left_VIcomb_uncorrected.txt','ROI_Left_VIord.txt','ROI_Left_VItemp.txt'};
        %               fileIn={'ROI_Right_VIcon.txt','ROI_Right_VIcomb_uncorrected.txt','ROI_Right_VIord.txt','ROI_Right_VItemp.txt'};
        fileIn={'ROI_RHcon.txt','ROI_RHcomb_uncorrected.txt','ROI_RHord.txt','ROI_RHtemp.txt'};
        %    fileIn={'ROI_LHcon.txt','ROI_LHcomb_uncorrected.txt','ROI_LHord.txt','ROI_LHtemp.txt'};
        %         fileIn={'ROI_Vermis_VIcon.txt','ROI_Vermis_VIcomb_uncorrected.txt','ROI_Vermis_VIord.txt','ROI_Vermis_VItemp.txt'};
        %           fileIn={'ROI_Left_CrusIcon.txt','ROI_Left_CrusIcomb_uncorrected.txt','ROI_Left_CrusIord.txt','ROI_Left_CrusItemp.txt'};
        %  fileIn={'ROI_Right_VI_CrusIcon.txt','ROI_Right_VI_CrusIcomb_uncorrected.txt','ROI_Right_VI_CrusIord.txt','ROI_Right_VI_CrusItemp.txt'};
        %   fileIn={'ROI_Right_CrusIcon.txt','ROI_Right_CrusIcomb_uncorrected.txt','ROI_Right_CrusIord.txt','ROI_Right_CrusItemp.txt'};
        
        fontSize=  14;
        A=dload(fileIn{1});
        
        
        % T=tapply(A,{'subj','accName'},{A.mean,'mean','name','acc'},'subset',A.accName==2)
        T=tapply(A,{'subj','accName'},{A.mean,'mean','name','acc'},'subset',A.accName>1)
        Factors=[T.accName];
        results=anovaMixed(T.acc,T.subj,'within',Factors,{'condition'});
        
        figure;
        
        %Reorder conditions
        T.accName=T.accName*10;
        T.accName(T.accName==20)=2; %overall
        T.accName(T.accName==30)=5; %comb
        T.accName(T.accName==40)=4; %ord
        T.accName(T.accName==50)=3; %temp
        
        %Reorder conditions
        %       T.accName=T.accName*10;
        %       T.accName(T.accName==20)=2+4; %overall
        %       T.accName(T.accName==30)=5+4; %comb
        %       T.accName(T.accName==40)=4+4; %ord
        %       T.accName(T.accName==50)=3+4; %temp
        
        % acc={'wcmcon_0001.img','wcmmva160_zacc_combUncorrectedSUIT.nii','wcmmva160_zacc_combCorrectedSUIT.nii','Ord_oneout/wcmmva160_zacc_Ord_oneoutSUIT.nii','Temp_oneout/wcmmva160_zacc_Temp_oneoutSUIT.nii'};
        
        barplot(T.accName,T.acc,'split',T.accName,'leg',{'overall','integrated','order','timing'},'leglocation','northeast');
        
        axis([0 10 -0.5 1.5]);
        
        
        barplot([T.accName],T.acc,'split',T.accName,'style_bold','plotfcn','mean',...
            'errorcolor',{[0 0 0],[1 0 0],[0 0 1],[0 1 0]},'facecolor',{[0 0 0],[1 0 0],[0 0 1],[0 1 0]},...
            'edgecolor',{[0 0 0],[1 0 0],[0 0 1],[0 1 0]},'gapwidth',[0.3 0.1 0.1]);
        
        ylabel('Zacc (mean)','FontSize',fontSize);
        axis([0 12 0 1.5]);
        
        %         title(RegName);
        
        set(gcf,'PaperPosition',[1 1 5 2]); %just played around with the figure size
        
        wysiwyg;
        
        %T-test:
        
        bonf=4; %Bonferroni correction
        tailsTest=1;
        
        D=tapply(A,{'subj'},{T.acc,'nanmean','name','zacc'},'subset',T.accName==2);
        [t_combU,p_combU]=ttest(D.zacc,D.zacc,tailsTest,'onesample');
        
        D=tapply(A,{'subj'},{T.acc,'nanmean','name','zacc'},'subset',T.accName==3);
        [t_comb,p_comb]=ttest(D.zacc,D.zacc,tailsTest,'onesample');
        
        D=tapply(A,{'subj'},{T.acc,'nanmean','name','zacc'},'subset',T.accName==4);
        [t_ord,p_ord]=ttest(D.zacc,D.zacc,tailsTest,'onesample');
        
        
        D=tapply(A,{'subj'},{T.acc,'nanmean','name','zacc'},'subset',T.accName==5);
        [t_temp,p_temp]=ttest(D.zacc,D.zacc,tailsTest,'onesample');
        
        disp('One-sample t-test')
        stats=[t_combU,p_combU*bonf;...
            t_comb,p_comb*bonf;...
            t_ord,p_ord*bonf;...
            t_temp,p_temp*bonf];
        
        disp(stats);
        
        
        
        
        
        %%Corr with beh.
        
        cd(fullfile(behDir));
        Transfer=dload('transfer.txt');
        Transfer.comb=Transfer.comb*1000;
        Transfer.ord=Transfer.ord*1000;
        Transfer.temp=Transfer.temp*1000;
        
        
        
        
        bonf=1; %Bonferroni correction
        figure;
        fontSize=14;
        %exclude 1 subj
        acc=T.acc(T.accName==2);
        acc=acc(2:16);
        b(1)=subplot(1,3,1); [r_combU,p_combU]=corr(acc,Transfer.comb,'tail','right')
        scatterplot(acc,Transfer.comb,'regression','linear','printcorr'); title('Overall','FontSize',fontSize);
        
        acc=T.acc(T.accName==4);
        acc=acc(2:16);
        b(2)=subplot(1,3,2); [r_ord,p_ord]=corr(acc,Transfer.ord,'tail','right')
        scatterplot(acc,Transfer.ord,'regression','linear','printcorr'); title('Overall','FontSize',fontSize);
        
        acc=T.acc(T.accName==5);
        acc=acc(2:16);
        b(3)=subplot(1,3,3); [r_temp,p_temp]=corr(acc,Transfer.temp,'tail','right')
        scatterplot(acc,Transfer.temp,'regression','linear','printcorr'); title('Overall','FontSize',fontSize);
        
        
        stats=[r_combU,p_combU*bonf;...
            r_ord,p_ord*bonf;...
            r_temp,p_temp*bonf];
        disp(stats);
        
        linkaxes(b,'xy');
        %         axis([-0.5 3 -20 80]);
        set(gcf,'PaperPosition',[1 1 15 5]); %just played around with the figure size
        
        wysiwyg;
        
        ylabel('Transfer (diff. to bsl in ms)','FontSize',fontSize);
        xlabel('Zacc (mean)','FontSize',fontSize);
        
        
end


