function prepProd2_forceAna_RY(what)

% addpath(genpath('G:\projectsBackup\rhys\prepProd2\matlab')); %adjust
% addpath(genpath('D:\projects\toolboxes\tools')); %joern's extensions for spm
% addpath(genpath('D:\projects\toolboxes\userfun')); %joern's util tools (open source)
% addpath(genpath('D:\projects\toolboxes\RainCloudPlots-master')) %code for raincloud plots

subj_name={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
    's11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23',...
    's24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34','s35','s36',...
    's37','s38','s39','s40','s41','s42','s43','s44','s45','s46','s47','s48','s49',...
    's50','s51','s52','s53','s54','s55','s56','s57','s58','s59','s60'}; %% chronological without missing subject numbers, for later vector references

%%% Blocks:
% BN01-BN04 probePre blocks
% BN05-BN22 training Day1
% BN23-BN40 training Day2
% BN29-BN40 Production from memory only (performance check!)
% BN41-BN44 probePost blocks
% BN45-BN46 training refresher
% BN47-BN52 fMRI

%%% trialType coding:
% trialType==1 sequence (instructed & mem)
% trialType==2 catch
%%% mode coding:
% mode==1 sequence instructed
% mode==2 sequence mem
% mode==0 catch

%%% All usable
subj=[3,5,6,7,9,10,13,16,17,18,20,21,22,25,26,31,32,34,36,38,39,40,41,42]; %meet both criteria (interaction & error rate)

% subj=[2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,20,21,22,25,26,31,32,33,34,36,37,38,39,40,41,42]; %includes those that don't meet criteria but were scanned
% subj=[2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42]; %includes all, even those not scanned

% subj=40; %for single participants, edit accordingly

subjExp=ones(1,length(subj));

savePath='G:\projectsBackup\rhys\prepProd2\data\behavioural\forces\processed'; %output path
baseDir='G:\projectsBackup\rhys\prepProd2\data\behavioural\forces'; %location of data path

%% Blocks to Analyse
minBN=1;
maxBN=52;

minBNTest=47;
maxBNTest=52;

% analysis of the spatial and temporal accuracy of subject performance
B=[];
switch(what)
    
    case 'createGroupForceFile' %produces group file containing all forces during test phase - takes a while
        for i=subj
            for j=minBNTest:maxBNTest
                fname=sprintf('exp_BN%02d.mat',(j));
                fileIn=fullfile(baseDir,subj_name{i},fname);
                load(fileIn);
                
                B = addstruct(B,E);
            end
            
            disp(['subj ' num2str(i) ' loaded.'])
        end
        
        save([baseDir '\groupForceTest.mat'],'B','-v7.3')
        disp(['Group test file saved under ', baseDir '\groupForceTest.mat'])
        
    case 'noGoBaselineVForce' %T test of force in no-go trials, baseline versus fractal onset to last possible go cue appearance
        
        load([baseDir '\groupForceTest.mat'])
        
        loopCounter = 1;
        
        for i=subj
            probeForces = B.forces(B.trialType == 2 & B.BN>=47 & B.subj == i);
            
            for j=1:length(probeForces)
                baseline = probeForces{j}(:,1:5);
                baseline = baseline(1:500,:);
                baseMean(j,:) = mean(baseline);
                baseMeanOvr(j,:) = mean(baseMean(j,:));
                
                active = probeForces{j}(:,1:5);
                active = active(501:2483,:);
                activeMean(j,:) = mean(active);
                activeMeanOvr(j,:) = mean(activeMean(j,:));
            end
            
            subjBaseMean(loopCounter,:) = mean(baseMeanOvr);
            subjActiveMean(loopCounter,:) = mean(activeMeanOvr);
            
            loopCounter = loopCounter+1;
            
        end
        
        ttest(subjBaseMean, subjActiveMean,2,'paired')
%         computeCohen_d(subjBaseMean, subjActiveMean,'paired')
        subjMeans = [subjBaseMean subjActiveMean];
        xlswrite([savePath '\noGoBaselineVForce'], subjMeans)
        
    case 'goBaselineVForce' %As above, but for go trials
        
        load([baseDir '\groupForceTest.mat'])
        
        loopCounter = 1;
        
        for i=subj
            probeForces = B.forces(B.trialType == 1 & B.BN>=47 & B.subj == i);
            
            for j=1:length(probeForces)
                baseline = probeForces{j}(:,1:5);
                baseline = baseline(1:500,:);
                baseMean(j,:) = mean(baseline);
                baseMeanOvr(j,:) = mean(baseMean(j,:));
                
                active = probeForces{j}(:,1:5);
                active = active(501:1000,:);
                activeMean(j,:) = mean(active);
                activeMeanOvr(j,:) = mean(activeMean(j,:));
            end
            
            subjBaseMean(loopCounter,:) = mean(baseMeanOvr);
            subjActiveMean(loopCounter,:) = mean(activeMeanOvr);
            
            loopCounter = loopCounter+1;
            
        end
        
        ttest(subjBaseMean, subjActiveMean,2,'paired')
        computeCohen_d(subjBaseMean, subjActiveMean,'paired')
        disp(['baseline stdev: ' num2str(std(subjBaseMean))])
        disp(['baseline stdev: ' num2str(std(subjActiveMean))])
        subjMeans = [subjBaseMean subjActiveMean];
        xlswrite([savePath '\GoBaselineVForce'], subjMeans)
        
end
