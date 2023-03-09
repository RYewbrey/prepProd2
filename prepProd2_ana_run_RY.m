function A=prepProd2_ana_run_RY(what)

%%% Add before starting any scripts (comment in when pasting into command line):
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

savePath='G:\projectsBackup\rhys\prepProd2\data\behavioural\'; %output path

A=[];
for i=1:length(subj)
    disp(['Subj' num2str(subj(i))] );
    D=prepProd2_ana_RY(subj(i),savePath);
    D.expertise=repmat(subjExp(i),size(D.subj,1),1);
    A=addstruct(A,D);
end;

condName={'trained','temporal transfer','spatial transfer','new'};

switch(what)
    
    case 'plot_exampleForce'
        
        cd 'E:\projects\rhys\prepProd2\data\behavioural\forces\s40'
        load('exp_BN47.mat')
        
        figure
        colourScheme = [0 0 .8; 1 .2 .6; 0.4 .4 .4; 1 .6 .4; 0 .8 .2; ...
            0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];
        
        for i=1:length(colourScheme)
            plot(E.thresholdedSmoothedForces{20}(1:5000,i), 'LineWidth', 2, 'color', colourScheme(i,:))
            hold on
        end
        
        ylabel('Force (N)', 'fontSize', 12, 'fontName', 'Calibri')
        xlabel('Time (ms)', 'fontSize', 12, 'fontName', 'Calibri')
        
        set(gca, 'fontSize', 12, 'fontName', 'Calibri')
        
        drawline(1,'dir','horz','linestyle','- -')
        drawline([500 (E.cueTime(20,:)+500+E.cueDur(20))])
        set(gca, 'XTick', 0:500:5000, 'XTickLabel', -1.5:0.5:3.5)
        
    case 'plot_exampleForceCatch'
        
        cd 'E:\projects\rhys\prepProd2\data\behavioural\forces\s40'
        load('exp_BN47.mat')
        
        figure
        colourScheme = [0 0 .8; 1 .2 .6; 0.4 .4 .4; 1 .6 .4; 0 .8 .2; ...
            0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];
        
        for i=1:length(colourScheme(:,1))
            plot(E.thresholdedSmoothedForces{40}(1:4000,i), 'LineWidth', 2, 'color', colourScheme(i,:))
            hold on
        end
        
        ylabel('Force (N)', 'fontSize', 12, 'fontName', 'Calibri')
        xlabel('Time (ms)', 'fontSize', 12, 'fontName', 'Calibri')
        
        set(gca, 'fontSize', 12, 'fontName', 'Calibri')
        
        ylim([-0.5, 2.5])
        
        drawline(1,'dir','horz','linestyle','- -')
        drawline(500)
        set(gca, 'XTick', 0:500:5000, 'XTickLabel', -1.5:0.5:2.5)
        
    case 'outlier' %error rates before and after fMRI - comment accordingly
        
        %% Error performance check BEFORE fMRI!
        trialNr=length(A.errorFinger(A.BN>=29&A.BN<=40&A.subj==subj(1)&A.trialType==1&A.mode==2&A.trialType==1&A.mode==2));
        T=tapply(A,{'subj'},{A.errorFinger,'sum','name','errorFinger'},'subset',A.BN>=29&A.BN<=40&A.trialType==1&A.mode==2);
        T.errorFingerPerc=T.errorFinger./trialNr*100;
        
        outlierData(:,1)=T.subj;
        outlierData(:,2)=T.errorFingerPerc;
        disp(outlierData(:,2))
        cd(savePath);
        xlswrite('outlier_errorRate.xlsx',outlierData);
        
        %% Error performance check during fMRI
        % Check for error rates in test phase (following training >=BN47)
        trialNr=length(A.errorFinger(A.BN>=47&A.BN<=52&A.subj==subj(1)&A.trialType==1&A.mode==2));
        T=tapply(A,{'subj'},{A.errorFinger,'sum','name','errorFinger'},'subset',A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2);
        T.errorFingerPerc=T.errorFinger./trialNr*100;
        meanErrorFinger=nanmean(T.errorFingerPerc);
        sdErrorFinger=nanstd(T.errorFingerPerc);
        cutOff=meanErrorFinger+sdErrorFinger*3;
        outlier=T.subj(T.errorFingerPerc>cutOff); % cuttOff=
        disp(outlier); % subject ... max error rate ...%
        
        outlierData(:,1)=T.subj;
        outlierData(:,2)=T.errorFingerPerc;
        
        cd(savePath);
        xlswrite('outlier_errorRate_fMRI.xlsx',outlierData);
        
        
        
        T; 
        
    case 'performanceCheck' %to look for a significant timing x interval interaction at end of training
        %Overall performance check before fMRI
        
        %%%%%%%%%%%%%%%%%%%%%%% ONLY SELECT 1 SUBJECT AT TOP OF SCRIPT %%%%%%%%%%%%%%%%%%%%%%%
        
        
        T=tapply(A,{'subj'},{A.points,'sum','name','pointsAll'},'subset',A.BN>0&A.BN<=44); %sum of points after 2 training days
        trainPoints=[T.subj T.pointsAll];
        disp(trainPoints)
        
        I1S1T1 = A.timingInterval1(A.spatID==1 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I1S1T2 = A.timingInterval1(A.spatID==1 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I1S2T1 = A.timingInterval1(A.spatID==2 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I1S2T2 = A.timingInterval1(A.spatID==2 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        
        I2S1T1 = A.timingInterval2(A.spatID==1 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I2S1T2 = A.timingInterval2(A.spatID==1 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I2S2T1 = A.timingInterval2(A.spatID==2 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I2S2T2 = A.timingInterval2(A.spatID==2 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        
        I3S1T1 = A.timingInterval3(A.spatID==1 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I3S1T2 = A.timingInterval3(A.spatID==1 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I3S2T1 = A.timingInterval3(A.spatID==2 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I3S2T2 = A.timingInterval3(A.spatID==2 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        
        I4S1T1 = A.timingInterval4(A.spatID==1 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I4S1T2 = A.timingInterval4(A.spatID==1 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I4S2T1 = A.timingInterval4(A.spatID==2 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        I4S2T2 = A.timingInterval4(A.spatID==2 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=29 & A.BN<=40 & A.points >=1 & A.jitterMask==0);
        
        anovaOutput = padcat(I1S1T1, I1S1T2, I1S2T1, I1S2T2, I2S1T1, I2S1T2, I2S2T1, I2S2T2, I3S1T1, I3S1T2, I3S2T1, I3S2T2, I4S1T1, I4S1T2, I4S2T1, I4S2T2);
        
        cd('E:\projects\rhys\prepProd2\docs')
        
        xlswrite([subj_name{subj} '_performanceAnova'], anovaOutput)
        
    case 'performanceCheckfMRI' %as above, at end of fMRI
        
        %Overall performance check after fMRI, decide inclusion in report.
        
        %%%%%%%%%%%%%%%%%%%%%%% ONLY SELECT 1 SUBJECT AT TOP OF SCRIPT %%%%%%%%%%%%%%%%%%%%%%%
        T=tapply(A,{'subj'},{A.points,'sum','name','pointsAll'},'subset',A.BN>=47&A.BN<=52); %sum of points after 2 training days
        trainPoints=[T.subj T.pointsAll];
        disp(trainPoints)
        
        I1S1T1 = A.timingInterval1(A.spatID==1 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I1S1T2 = A.timingInterval1(A.spatID==1 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I1S2T1 = A.timingInterval1(A.spatID==2 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I1S2T2 = A.timingInterval1(A.spatID==2 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        
        I2S1T1 = A.timingInterval2(A.spatID==1 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I2S1T2 = A.timingInterval2(A.spatID==1 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I2S2T1 = A.timingInterval2(A.spatID==2 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I2S2T2 = A.timingInterval2(A.spatID==2 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        
        I3S1T1 = A.timingInterval3(A.spatID==1 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I3S1T2 = A.timingInterval3(A.spatID==1 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I3S2T1 = A.timingInterval3(A.spatID==2 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I3S2T2 = A.timingInterval3(A.spatID==2 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        
        I4S1T1 = A.timingInterval4(A.spatID==1 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I4S1T2 = A.timingInterval4(A.spatID==1 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I4S2T1 = A.timingInterval4(A.spatID==2 & A.tempID==1 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        I4S2T2 = A.timingInterval4(A.spatID==2 & A.tempID==2 & A.trialType==1 & A.mode==2 & A.BN>=47 & A.BN<=52 & A.points >=1 & A.jitterMask==0);
        
        anovaOutput = padcat(I1S1T1, I1S1T2, I1S2T1, I1S2T2, I2S1T1, I2S1T2, I2S2T1, I2S2T2, I3S1T1, I3S1T2, I3S2T1, I3S2T2, I4S1T1, I4S1T2, I4S2T1, I4S2T2);
        
        cd('E:\projects\rhys\prepProd2\docs')
        
        xlswrite([subj_name{subj} '_performanceAnovafMRI'], anovaOutput) 
        
    case 'intervalProd' %interval production in from memory trials during test block, do they modulate timings as a group?
        T=[];
        colorScheme={[.9 0 0], [0 .9 0],[.4 0 0],[0 .4 0]};
        monoColourScheme={[.9 .6 0], [0 0 0], [0 .6 .5], [.8 .5 .7]};
        %Take mean duration of each IPI across sequences from memory in testing phase
                
        T1=tapply(A,{'subj','seqID'},{A.timingInterval1,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        T2=tapply(A,{'subj','seqID'},{A.timingInterval2,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        T3=tapply(A,{'subj','seqID'},{A.timingInterval3,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        T4=tapply(A,{'subj','seqID'},{A.timingInterval4,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
%         T1=tapply(A,{'subj','seqID'},{A.timingInterval1,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
%         T2=tapply(A,{'subj','seqID'},{A.timingInterval2,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
%         T3=tapply(A,{'subj','seqID'},{A.timingInterval3,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
%         T4=tapply(A,{'subj','seqID'},{A.timingInterval4,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
        T1.intervalPos=ones(size(T1.subj,1),1);
        T2.intervalPos=ones(size(T2.subj,1),1)*2;
        T3.intervalPos=ones(size(T3.subj,1),1)*3;
        T4.intervalPos=ones(size(T4.subj,1),1)*4;
        T=addstruct(T,T1);
        T=addstruct(T,T2);
        T=addstruct(T,T3);
        T=addstruct(T,T4);
        
        K1=tapply(A,{'subj','seqID'},{A.intervalPrc1,'nanmean','name','intervalPrc1'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        K2=tapply(A,{'subj','seqID'},{A.intervalPrc2,'nanmean','name','intervalPrc2'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        K3=tapply(A,{'subj','seqID'},{A.intervalPrc3,'nanmean','name','intervalPrc3'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        K4=tapply(A,{'subj','seqID'},{A.intervalPrc4,'nanmean','name','intervalPrc4'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
%         K1=tapply(A,{'subj','seqID'},{A.intervalPrc1,'nanmean','name','intervalPrc1'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
%         K2=tapply(A,{'subj','seqID'},{A.intervalPrc2,'nanmean','name','intervalPrc2'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
%         K3=tapply(A,{'subj','seqID'},{A.intervalPrc3,'nanmean','name','intervalPrc3'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
%         K4=tapply(A,{'subj','seqID'},{A.intervalPrc4,'nanmean','name','intervalPrc4'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
        J1.intervalPrc=K1.intervalPrc1;
        J2.intervalPrc=K2.intervalPrc2;
        J3.intervalPrc=K3.intervalPrc3;
        J4.intervalPrc=K4.intervalPrc4;
        T=addstruct(T,J1);
        T=addstruct(T,J2);
        T=addstruct(T,J3);
        T=addstruct(T,J4);
        
        K5=tapply(A,{'subj','seqID'},{A.targetIntervalPrc1,'nanmean','name','targetIntervalPrc1'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        K6=tapply(A,{'subj','seqID'},{A.targetIntervalPrc2,'nanmean','name','targetIntervalPrc2'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        K7=tapply(A,{'subj','seqID'},{A.targetIntervalPrc3,'nanmean','name','targetIntervalPrc3'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        K8=tapply(A,{'subj','seqID'},{A.targetIntervalPrc4,'nanmean','name','targetIntervalPrc4'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
%         K5=tapply(A,{'subj','seqID'},{A.targetIntervalPrc1,'nanmean','name','targetIntervalPrc1'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
%         K6=tapply(A,{'subj','seqID'},{A.targetIntervalPrc2,'nanmean','name','targetIntervalPrc2'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
%         K7=tapply(A,{'subj','seqID'},{A.targetIntervalPrc3,'nanmean','name','targetIntervalPrc3'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
%         K8=tapply(A,{'subj','seqID'},{A.targetIntervalPrc4,'nanmean','name','targetIntervalPrc4'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
        J5.targetIntervalPrc=K5.targetIntervalPrc1;
        J6.targetIntervalPrc=K6.targetIntervalPrc2;
        J7.targetIntervalPrc=K7.targetIntervalPrc3;
        J8.targetIntervalPrc=K8.targetIntervalPrc4;
        T=addstruct(T,J5);
        T=addstruct(T,J6);
        T=addstruct(T,J7);
        T=addstruct(T,J8);
        
        
        figure;
        
%         monoColourScheme = {[0 0 0], [0.9 0.9 0.9], [0.87 0.72 0.52], [0.5 0 0]};
        
        %Uncorrected
        lineplot(T.intervalPos,T.timing,'style_thickline',...
            'split',T.seqID,'linecolor',monoColourScheme,'errorcolor',monoColourScheme,'markerfill',monoColourScheme,'markercolor',monoColourScheme,...
            'leg', 'auto','leglocation','northeast','linewidth',3,'errorwidth',3);
        title('Uncorrected', 'FontSize',12, 'FontName','Calibri');
        ylabel('Interval (ms)', 'FontSize',12, 'FontName','Calibri');
        xlabel('Interval position', 'FontSize',12, 'FontName','Calibri');
        ylim([200 1000])
        set(gca,'FontSize',12, 'FontName','Calibri')
        
        
        %Mean corrected plots
        M=tapply(T,{'subj'},{T.timing,'nanmean','name','meanTiming'});
        Mmean=nanmean(M.meanTiming);
        C=T;
        for s=1:length(subj)
            C.timing(C.subj==subj(s))=C.timing(C.subj==subj(s))-M.meanTiming(M.subj==subj(s));
            C.timing(C.subj==subj(s),:)=C.timing(C.subj==subj(s),:)-M.meanTiming(M.subj==subj(s),:)+Mmean;
        end;
        
        figure;
        lineplot(C.intervalPos,C.timing,'style_thickline',...
            'split',C.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        ylabel('Interval (ms) Mean-corrected');
        xlabel('Interval position');
        ylim([-400 1200]);
        set(gca,'FontSize',16)
        
        %% Output for SPSS
        condTiming=4;
        seqID=4;
        condCol=condTiming*seqID;
        subjRows=length(subj);
        S=nan(subjRows,condCol);
        
        S(:,1)=T1.timing(T1.seqID==1,1);
        S(:,2)=T1.timing(T1.seqID==2,1);
        S(:,3)=T1.timing(T1.seqID==3,1);
        S(:,4)=T1.timing(T1.seqID==4,1);
        S(:,5)=T2.timing(T2.seqID==1,1);
        S(:,6)=T2.timing(T2.seqID==2,1);
        S(:,7)=T2.timing(T2.seqID==3,1);
        S(:,8)=T2.timing(T2.seqID==4,1);
        S(:,9)=T3.timing(T3.seqID==1,1);
        S(:,10)=T3.timing(T3.seqID==2,1);
        S(:,11)=T3.timing(T3.seqID==3,1);
        S(:,12)=T3.timing(T3.seqID==4,1);
        S(:,13)=T4.timing(T4.seqID==1,1);
        S(:,14)=T4.timing(T4.seqID==2,1);
        S(:,15)=T4.timing(T4.seqID==3,1);
        S(:,16)=T4.timing(T4.seqID==4,1);
        
        S=[ subj' S ];
        cd(savePath);
        
        %Save interval percent per IPI per delay:
        %         S(:,17)=T.intervalPrc(T.intervalPos==1&T.seqID==1,1);
        %         S(:,18)=T.intervalPrc(T.intervalPos==1&T.seqID==2,1);
        %         S(:,19)=T.intervalPrc(T.intervalPos==1&T.seqID==3,1);
        %         S(:,20)=T.intervalPrc(T.intervalPos==1&T.seqID==4,1);
        %         S(:,21)=T.intervalPrc(T.intervalPos==2&T.seqID==1,1);
        %         S(:,22)=T.intervalPrc(T.intervalPos==2&T.seqID==2,1);
        %         S(:,23)=T.intervalPrc(T.intervalPos==2&T.seqID==3,1);
        %         S(:,24)=T.intervalPrc(T.intervalPos==2&T.seqID==4,1);
        %         S(:,25)=T.intervalPrc(T.intervalPos==3&T.seqID==1,1);
        %         S(:,26)=T.intervalPrc(T.intervalPos==3&T.seqID==2,1);
        %         S(:,27)=T.intervalPrc(T.intervalPos==3&T.seqID==3,1);
        %         S(:,28)=T.intervalPrc(T.intervalPos==3&T.seqID==4,1);
        %         S(:,29)=T.intervalPrc(T.intervalPos==4&T.seqID==1,1);
        %         S(:,30)=T.intervalPrc(T.intervalPos==4&T.seqID==2,1);
        %         S(:,31)=T.intervalPrc(T.intervalPos==4&T.seqID==3,1);
        %         S(:,32)=T.intervalPrc(T.intervalPos==4&T.seqID==4,1);
        
        xlswrite('IntervalProd.xlsx',S);
        
        
        %Plot interval prc as a function of each IPI, split by delay:
        figure;
        colorScheme={[.7 .7 .7], [.45 .45 .45], [0 0 0]};
        lineplot(T.intervalPos,T.intervalPrc,'style_thickline',...
            'split',T.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        xlabel('Interval position');
        ylabel('Percent of mean interval in long delay');
        ylim([20 250]);
        set(gca,'FontSize',16)
        %Plot target interval prc as a function of each IPI, split by timing condition:
        figure;
        lineplot(T.intervalPos,T.targetIntervalPrc,'style_thickline',...
            'split',T.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        xlabel('Interval position');
        ylabel('Target percent of interval in long delay');
        ylim([20 250]);
        set(gca,'FontSize',16)
        
        T;
        
    case 'intervalProdBySubj' %interval production in from memory trials during test block, timing modulation by subject
        T=[];
        colorScheme={[.9 0 0], [0 .9 0],[.4 0 0],[0 .4 0]};
        %Take mean duration of each IPI across sequences from memory in testing phase
        figure;
        for i=1:length(unique(A.subj))
            T=[];
            figure %for separate plots
            
            %raw intervals
            T1.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T1.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T1.timing = A.timingInterval1(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            T2.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T2.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T2.timing = A.timingInterval2(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            T3.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T3.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T3.timing = A.timingInterval3(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            T4.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T4.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T4.timing = A.timingInterval4(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            T1.intervalPos=ones(size(T1.subj,1),1);
            T2.intervalPos=ones(size(T2.subj,1),1)*2;
            T3.intervalPos=ones(size(T3.subj,1),1)*3;
            T4.intervalPos=ones(size(T4.subj,1),1)*4;
            T=addstruct(T,T1);
            T=addstruct(T,T2);
            T=addstruct(T,T3);
            T=addstruct(T,T4);
            
            % intervals as a percentage
            K1.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K1.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K1.intervalPrc1 = A.intervalPrc1(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K2.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K2.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K2.intervalPrc2 = A.intervalPrc2(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K3.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K3.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K3.intervalPrc3 = A.intervalPrc3(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K4.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K4.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K4.intervalPrc4 = A.intervalPrc4(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            J1.intervalPrc=K1.intervalPrc1;
            J2.intervalPrc=K2.intervalPrc2;
            J3.intervalPrc=K3.intervalPrc3;
            J4.intervalPrc=K4.intervalPrc4;
            T=addstruct(T,J1);
            T=addstruct(T,J2);
            T=addstruct(T,J3);
            T=addstruct(T,J4);
            
            %target intervals as a percentage
            K5.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K5.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K5.targetIntervalPrc1 = A.targetIntervalPrc1(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K6.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K6.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K6.targetIntervalPrc2 = A.targetIntervalPrc2(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K7.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K7.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K7.targetIntervalPrc3 = A.targetIntervalPrc3(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K8.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K8.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K8.targetIntervalPrc4 = A.targetIntervalPrc4(A.trialType==1&A.mode==2&A.BN>=29 & A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            J5.targetIntervalPrc=K5.targetIntervalPrc1;
            J6.targetIntervalPrc=K6.targetIntervalPrc2;
            J7.targetIntervalPrc=K7.targetIntervalPrc3;
            J8.targetIntervalPrc=K8.targetIntervalPrc4;
            T=addstruct(T,J5);
            T=addstruct(T,J6);
            T=addstruct(T,J7);
            T=addstruct(T,J8);
            
            
            
            
            %Uncorrected
%             subplot(6,6,i) %comment out for multiple individual figures
            lineplot(T.intervalPos,T.timing,'style_thickline',...
                'split',T.seqID, 'errorfcn', 'std','linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
            title(['Subj ', num2str(subj(i))]);
            ylabel('Interval (ms)');
            xlabel('Interval position');
            ylim([0 1500]);
            set(gca,'FontSize',10)
            
        end;
        
        %Mean corrected plots
        M=tapply(T,{'subj'},{T.timing,'nanmean','name','meanTiming'});
        Mmean=nanmean(M.meanTiming);
        C=T;
        for s=1:length(subj)
            C.timing(C.subj==subj(s),:)=C.timing(C.subj==subj(s),:)-M.meanTiming(M.subj==subj(s),:);
            C.timing(C.subj==subj(s),:)=C.timing(C.subj==subj(s),:)-M.meanTiming(M.subj==subj(s),:)+Mmean;
        end;
        
        figure;
        lineplot(C.intervalPos,C.timing,'style_thickline',...
            'split',C.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        ylabel('Interval (ms) Mean-corrected');
        xlabel('Interval position');
        ylim([-400 1200]);
        set(gca,'FontSize',16)
        
        
        %Plot interval prc as a function of each IPI, split by delay:
        figure;
        colorScheme={[.7 .7 .7], [.45 .45 .45], [0 0 0]};
        lineplot(T.intervalPos,T.intervalPrc,'style_thickline',...
            'split',T.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        xlabel('Interval position');
        ylabel('Percent of mean interval in long delay');
        ylim([20 250]);
        set(gca,'FontSize',16)
        %Plot target interval prc as a function of each IPI, split by timing condition:
        figure;
        lineplot(T.intervalPos,T.targetIntervalPrc,'style_thickline',...
            'split',T.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        xlabel('Interval position');
        ylabel('Target percent of interval in long delay');
        ylim([20 250]);
        set(gca,'FontSize',16)
        
        T;
        
    case 'intervalProdBySubjfMRI' %interval production in from memory trials during test block, timing modulation by subject
        T=[];
        colorScheme={[.9 0 0], [0 .9 0],[.4 0 0],[0 .4 0]};
        %Take mean duration of each IPI across sequences from memory in testing phase
        figure;
        for i=1:length(unique(A.subj))
            T=[];
            %             figure %for separate plots
            
            %raw intervals
            T1.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T1.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T1.timing = A.timingInterval1(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            T2.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T2.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T2.timing = A.timingInterval2(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            T3.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T3.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T3.timing = A.timingInterval3(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            T4.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T4.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T4.timing = A.timingInterval4(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            T1.intervalPos=ones(size(T1.subj,1),1);
            T2.intervalPos=ones(size(T2.subj,1),1)*2;
            T3.intervalPos=ones(size(T3.subj,1),1)*3;
            T4.intervalPos=ones(size(T4.subj,1),1)*4;
            T=addstruct(T,T1);
            T=addstruct(T,T2);
            T=addstruct(T,T3);
            T=addstruct(T,T4);
            
            % intervals as a percentage
            K1.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K1.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K1.intervalPrc1 = A.intervalPrc1(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K2.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K2.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K2.intervalPrc2 = A.intervalPrc2(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K3.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K3.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K3.intervalPrc3 = A.intervalPrc3(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K4.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K4.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K4.intervalPrc4 = A.intervalPrc4(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            J1.intervalPrc=K1.intervalPrc1;
            J2.intervalPrc=K2.intervalPrc2;
            J3.intervalPrc=K3.intervalPrc3;
            J4.intervalPrc=K4.intervalPrc4;
            T=addstruct(T,J1);
            T=addstruct(T,J2);
            T=addstruct(T,J3);
            T=addstruct(T,J4);
            
            %target intervals as a percentage
            K5.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K5.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K5.targetIntervalPrc1 = A.targetIntervalPrc1(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K6.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K6.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K6.targetIntervalPrc2 = A.targetIntervalPrc2(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K7.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K7.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K7.targetIntervalPrc3 = A.targetIntervalPrc3(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            K8.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K8.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K8.targetIntervalPrc4 = A.targetIntervalPrc4(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            
            J5.targetIntervalPrc=K5.targetIntervalPrc1;
            J6.targetIntervalPrc=K6.targetIntervalPrc2;
            J7.targetIntervalPrc=K7.targetIntervalPrc3;
            J8.targetIntervalPrc=K8.targetIntervalPrc4;
            T=addstruct(T,J5);
            T=addstruct(T,J6);
            T=addstruct(T,J7);
            T=addstruct(T,J8);
            
            
            
            
            %Uncorrected
            subplot(5,5,i) %comment out for multiple individual figures
            lineplot(T.intervalPos,T.timing,'style_thickline',...
                'split',T.seqID, 'errorfcn', 'std','linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
            title(['Subj ', num2str(subj(i))]);
            ylabel('Interval (ms)');
            xlabel('Interval position');
            ylim([0 1500]);
            set(gca,'FontSize',10)
            
        end;
        
        %Mean corrected plots
        M=tapply(T,{'subj'},{T.timing,'nanmean','name','meanTiming'});
        Mmean=nanmean(M.meanTiming);
        C=T;
        for s=1:length(subj)
            C.timing(C.subj==subj(s),:)=C.timing(C.subj==subj(s),:)-M.meanTiming(M.subj==subj(s),:);
            C.timing(C.subj==subj(s),:)=C.timing(C.subj==subj(s),:)-M.meanTiming(M.subj==subj(s),:)+Mmean;
        end;
        
        figure;
        lineplot(C.intervalPos,C.timing,'style_thickline',...
            'split',C.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        ylabel('Interval (ms) Mean-corrected');
        xlabel('Interval position');
        ylim([-400 1200]);
        set(gca,'FontSize',16)
        
        
        %Plot interval prc as a function of each IPI, split by delay:
        figure;
        colorScheme={[.7 .7 .7], [.45 .45 .45], [0 0 0]};
        lineplot(T.intervalPos,T.intervalPrc,'style_thickline',...
            'split',T.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        xlabel('Interval position');
        ylabel('Percent of mean interval in long delay');
        ylim([20 250]);
        set(gca,'FontSize',16)
        %Plot target interval prc as a function of each IPI, split by timing condition:
        figure;
        lineplot(T.intervalPos,T.targetIntervalPrc,'style_thickline',...
            'split',T.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        xlabel('Interval position');
        ylabel('Target percent of interval in long delay');
        ylim([20 250]);
        set(gca,'FontSize',16)
        
        T;
        
    case 'performanceAnova' %interval production in from memory trials during training, timing modulation by subject
        T=[];
        colorScheme={[.9 0 0], [0 .9 0],[.4 0 0],[0 .4 0]};
        %Take mean duration of each IPI across sequences from memory in testing phase
        figure;
        for i=1:length(unique(A.subj))
            T1=tapply(A,{'subj','seqID'},{A.timingInterval1,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T2=tapply(A,{'subj','seqID'},{A.timingInterval2,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T3=tapply(A,{'subj','seqID'},{A.timingInterval3,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T4=tapply(A,{'subj','seqID'},{A.timingInterval4,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            T1.intervalPos=ones(size(T1.subj,1),1);
            T2.intervalPos=ones(size(T2.subj,1),1)*2;
            T3.intervalPos=ones(size(T3.subj,1),1)*3;
            T4.intervalPos=ones(size(T4.subj,1),1)*4;
            T=addstruct(T,T1);
            T=addstruct(T,T2);
            T=addstruct(T,T3);
            T=addstruct(T,T4);
            
            K1=tapply(A,{'subj','seqID'},{A.intervalPrc1,'nanmean','name','intervalPrc1'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K2=tapply(A,{'subj','seqID'},{A.intervalPrc2,'nanmean','name','intervalPrc2'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K3=tapply(A,{'subj','seqID'},{A.intervalPrc3,'nanmean','name','intervalPrc3'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K4=tapply(A,{'subj','seqID'},{A.intervalPrc4,'nanmean','name','intervalPrc4'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            J1.intervalPrc=K1.intervalPrc1;
            J2.intervalPrc=K2.intervalPrc2;
            J3.intervalPrc=K3.intervalPrc3;
            J4.intervalPrc=K4.intervalPrc4;
            T=addstruct(T,J1);
            T=addstruct(T,J2);
            T=addstruct(T,J3);
            T=addstruct(T,J4);
            
            K5=tapply(A,{'subj','seqID'},{A.targetIntervalPrc1,'nanmean','name','targetIntervalPrc1'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K6=tapply(A,{'subj','seqID'},{A.targetIntervalPrc2,'nanmean','name','targetIntervalPrc2'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K7=tapply(A,{'subj','seqID'},{A.targetIntervalPrc3,'nanmean','name','targetIntervalPrc3'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            K8=tapply(A,{'subj','seqID'},{A.targetIntervalPrc4,'nanmean','name','targetIntervalPrc4'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            J5.targetIntervalPrc=K5.targetIntervalPrc1;
            J6.targetIntervalPrc=K6.targetIntervalPrc2;
            J7.targetIntervalPrc=K7.targetIntervalPrc3;
            J8.targetIntervalPrc=K8.targetIntervalPrc4;
            T=addstruct(T,J5);
            T=addstruct(T,J6);
            T=addstruct(T,J7);
            T=addstruct(T,J8);
            
            
            
            
            %Uncorrected
            subplot(3,3,i)
            lineplot(T.intervalPos,T.timing,'style_thickline',...
                'split',T.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
            title(['Subj ', num2str(subj(i))]);
            ylabel('Interval (ms)');
            xlabel('Interval position');
            ylim([100 1400]);
            set(gca,'FontSize',16)
            
            
        end;
        
        %Mean corrected plots
        M=tapply(T,{'subj'},{T.timing,'nanmean','name','meanTiming'});
        Mmean=nanmean(M.meanTiming);
        C=T;
        for s=1:length(subj)
            C.timing(C.subj==subj(s),:)=C.timing(C.subj==subj(s),:)-M.meanTiming(M.subj==subj(s),:);
            C.timing(C.subj==subj(s),:)=C.timing(C.subj==subj(s),:)-M.meanTiming(M.subj==subj(s),:)+Mmean;
        end;
        
        figure;
        lineplot(C.intervalPos,C.timing,'style_thickline',...
            'split',C.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        ylabel('Interval (ms) Mean-corrected');
        xlabel('Interval position');
        ylim([-400 1200]);
        set(gca,'FontSize',16)
        
        
        %Plot interval prc as a function of each IPI, split by delay:
        figure;
        colorScheme={[.7 .7 .7], [.45 .45 .45], [0 0 0]};
        lineplot(T.intervalPos,T.intervalPrc,'style_thickline',...
            'split',T.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        xlabel('Interval position');
        ylabel('Percent of mean interval in long delay');
        ylim([20 250]);
        set(gca,'FontSize',16)
        %Plot target interval prc as a function of each IPI, split by timing condition:
        figure;
        lineplot(T.intervalPos,T.targetIntervalPrc,'style_thickline',...
            'split',T.seqID,'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
        xlabel('Interval position');
        ylabel('Target percent of interval in long delay');
        ylim([20 250]);
        set(gca,'FontSize',16)
        
        T;
        
    case 'jitter' %plot delay between actual and planned cue appearance - for cogent lag
        %         figure;
        %subplot(1,3,1);
        lineplot(A.TNall,A.timingJitter,'style_thickline',...
            'plotfcn','mean','leg','auto','leglocation','northeast');
        title('Deviation from target interval ');
        xlabel('TN');
        ylabel('Abs Jitter (ms)');
        
        figure;
        for i=1:length(unique(A.subj))
            subplot(length(unique(A.subj)),1,i);
            scatter(A.TNall(A.jitterMask==0&A.subj==subj(i)),A.timingJitter(A.jitterMask==0&A.subj==subj(i)),'blue');
            hold on;
            scatter(A.TNall(A.jitterMask==1&A.subj==subj(i)),A.timingJitter(A.jitterMask==1&A.subj==subj(i)),'red');
            ylim([0 50]);
            sumTN=length(find(A.subj==subj(i)&A.BN<=44)); %number of trials beh. sessions (trained seq only)
            affTN=length(find(A.subj==subj(i)&A.BN<=44&A.jitterMask==1));
            affTNperc=affTN/sumTN*100;
            
            title(['Deviation from target interval, Subj: ' num2str(subj(i)) ', perc:' num2str(affTNperc)]);
            xlabel('TN');
            drawline(1104,'dir','vert');
            xlim([0 768]);
            
            %     ylabel(['Abs Jitter (ms) Subj' num2str(subj(i)) ]);
        end;
        
        A;
        
    case 'prepost' %plot learning across conditions - trained, trained timing, trained order, new
        colorScheme={[0 1 0],[1 0 0],[0 0 1],[0 0 0]}; %split by seqCond - green=trained,red=trained timing,blue=trained finger order,black=new
        %         %RT interval percent deviation
        %         figure;
        %         for i=1:length(unique(A.subj))
        %             subplot(1,length(unique(A.subj)),i);
        %             lineplot(A.prepost,A.RTabsInterval,'split',A.seqCond,'errorfcn','std','subset',A.subj==subj(i)&A.trialType==1&A.prepost>0&A.points>=1&A.jitterMask==0,'style_thickline',...
        %                 'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        %             ylabel('RTabs Interval (%)');
        %         end;
        %         figure;
        %         T=tapply(A,{'subj','seqCond','prepost'},{A.RTabsInterval,'nanmean','name','RTabsInterval'},'subset',A.trialType==1&A.prepost>0&A.points>=1&A.jitterMask==0);
        %         lineplot(T.prepost,T.RTabsInterval,'split',T.seqCond,'style_thickline',...
        %             'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        %         ylabel('RTabs Interval (%)');
        
        %RTabs deviation ms after first 3 repetitions
        %         figure;
        %         for i=1:length(unique(A.subj))
        %             subplot(1,length(unique(A.subj)),i);
        %             lineplot(A.prepost,A.RTabs,'split',A.seqCond,'errorfcn','std','subset',A.subj==subj(i)&A.trialType==1&A.prepost>0&A.trialRep<4&A.points>=1&A.jitterMask==0,'style_thickline',...
        %                 'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        %             ylabel('RTabs ms');
        %         end;
        figure;
        T=tapply(A,{'subj','seqCond','prepost','trialRep'},{A.RTabs,'nanmean','name','RTabs'},'subset',A.trialType==1&A.prepost>0&A.trialRep>0&A.points>=1&A.jitterMask==0);
        lineplot([T.prepost T.trialRep],T.RTabs,'split',T.seqCond,'style_thickline',...
            'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        ylabel('RTabs deviation (ms)');
        
        
        figure;
        %         T=tapply(A,{'subj','seqCond','prepost','trialRep'},{A.RT,'nanmean','name','RT'},'subset',A.trialType==1&A.prepost>0&A.trialRep>0&A.points>=1&A.jitterMask==0);
        T=tapply(A,{'subj','seqCond','prepost','trialRep'},{A.RT,'nanmedian','name','RT'},'subset',A.trialType==1&A.prepost>0&A.trialRep>0&A.points>=1&A.jitterMask==0);
        lineplot([T.prepost T.trialRep],T.RT,'split',T.seqCond,'style_thickline',...
            'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        ylabel('RT (ms)');
        
        
        
        %%%Do a loop across subj
        figure;
        for i=1:length(subj)
            subplot(5,5,i);
            T=tapply(A,{'subj','seqCond','prepost','trialRep'},{A.RTabs,'nanmean','name','RTabs'},'subset',A.trialType==1&A.prepost>0&A.trialRep>0&A.points>=1&A.jitterMask==0&A.subj==subj(i));
            lineplot([T.prepost T.trialRep],T.RTabs,'split',T.seqCond,'style_thickline',...
                'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto','leglocation','northeast');
            ylabel('RTabs deviation (ms)');
        end
        
        
        
        %Split by repetition
        figure;
        T=tapply(A,{'subj','seqCond','prepost'},{A.RTabs,'nanmean','name','RTabs'},'subset',A.trialType==1&A.prepost>0&A.trialRep>3&A.points>=1&A.jitterMask==0);
        lineplot(T.prepost,T.RTabs,'split',T.seqCond,'style_thickline',...
            'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        title('After first 3 repetitions');
        ylabel('RTabs deviation (ms)');
        
        
        %errorRate after first 3 repetitions
        %         figure;
        %         for i=1:length(unique(A.subj))
        %             subplot(1,length(unique(A.subj)),i);
        %             lineplot(A.prepost,A.errorRate,'split',A.seqCond,'errorfcn','std','subset',A.subj==subj(i)&A.trialType==1&A.prepost>0&A.trialRep<4&A.points>=1&A.jitterMask==0,'style_thickline',...
        %                 'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        %             ylabel('Error (%)');
        %         end;
        %         figure;
        %         T=tapply(A,{'subj','seqCond','prepost'},{A.errorRate,'nanmean','name','errorRate'},'subset',A.trialType==1&A.prepost>0&A.trialRep<4&A.points>=1&A.jitterMask==0);
        %         lineplot(T.prepost,T.errorRate,'split',T.seqCond,'style_thickline',...
        %             'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        %             ylabel('Error (%)');
        
        %errorRate
        figure;
        for i=1:length(unique(A.subj))
            subplot(1,length(unique(A.subj)),i);
            lineplot(A.prepost,A.errorRate,'split',A.seqCond,'errorfcn','std','subset',A.subj==subj(i)&A.trialType==1&A.prepost>0&A.points>=1&A.jitterMask==0,'style_thickline',...
                'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
            ylabel('Error (%)');
        end;
        figure;
        T=tapply(A,{'subj','seqCond','prepost'},{A.errorRate,'nanmean','name','errorRate'},'subset',A.trialType==1&A.prepost>0&A.points>=1&A.jitterMask==0);
        lineplot(T.prepost,T.errorRate,'split',T.seqCond,'style_thickline',...
            'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        ylabel('Error (%)');
        T;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% This section switches timing and order around, to align with the rest of the figures in the paper%%%%%%%%%%%%%%%
        
        figure
        T=tapply(A,{'subj','seqCond'},{A.RTabs,'nanmean','name','RTabs'},'subset',A.trialType==1&A.prepost==2&A.trialRep>=4&A.points>=1&A.jitterMask==0);
        twos = T.seqCond==2; %SWAP ORDER AND TIMING AROUND IN STRUCTURE
        threes = T.seqCond==3;
        T.seqCond(twos) = 3;
        T.seqCond(threes) = 2;
        
        myboxplot(T.seqCond,T.RTabs, 'xtickoff')
        set(gca,'xticklabel',{'Trained','Trained Order', 'Trained Timing','New'})
        
        ylabel('Absolute timing deviation from target (ms)','FontName', 'Calibri', 'FontSize', 12)
        set(gca,'FontSize', 12, 'FontName', 'calibri')
        
        %to plot connecting lines between participant data
%         hold on
%         dataConcat = [T.RTabs(T.seqCond ==1), T.RTabs(T.seqCond ==2), T.RTabs(T.seqCond ==3), T.RTabs(T.seqCond ==4)];
%         
%         for i = 1:length(subj) %plot lines that connect each participant's data point
% %             plot([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'color', [.5, .5, .5], 'lineWidth', 1.2) %solid grey lines
%             patchline([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
% %             plot([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
%             patchline([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
% %             plot([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
%             patchline([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
%         end
        
    case 'prepostExtract'%as above but produce excel sheet for spss
        %%% output for SPSS analyses
        %synchonisation rep ANOVA
        maxSeqCond=4;   %trained, timing, order, new
        prepost=2;  %prePost
        subjRows=length(subj);
        condCol=maxSeqCond*prepost;
        S=nan(subjRows,condCol);
        k=0; %counter
        for j=1:prepost %seq condition
            for i=1:maxSeqCond
                k=k+1;
                K=tapply(A,{'subj'},{A.RTabs,'nanmedian','name','RTabs'},'subset',A.seqCond==i&A.prepost==j&A.trialType==1&A.points>=1&A.trialRep>=4&A.jitterMask==0);
                S(:,k)=K.RTabs;
            end;
        end;
        S=[ subj' S ];
        cd(savePath);
        xlswrite('PrePost_RTabs.xlsx',S);
        
    case 'prepostExtract_RT'%as above but RT instead of RTabs
        %%% output for SPSS analyses
        %synchonisation rep ANOVA
        maxSeqCond=4;   %trained, temp, spat, new
        prepost=2;  %prePost
        subjRows=length(subj);
        condCol=maxSeqCond*prepost;
        S=nan(subjRows,condCol);
        k=0; %counter
        for j=1:prepost %seq condition
            for i=1:maxSeqCond
                k=k+1;
                K=tapply(A,{'subj'},{A.RT,'nanmean','name','RT'},'subset',A.seqCond==i&A.prepost==j&A.trialType==1&A.points>=1&A.trialRep>=4&A.jitterMask==0);
                S(:,k)=K.RT;
            end;
        end;
        S=[ subj' S ];
        cd(savePath);
        xlswrite('PrePost_RT.xlsx',S);
        
    case 'training' % Day1 vs. Day2
        colorScheme={[0 0 0]}; %
        %RT interval percent deviation
        %         figure;
        %         for i=1:length(unique(A.subj))
        %             subplot(1,length(unique(A.subj)),i);
        %             lineplot(A.training,A.RTabsInterval,'errorfcn','std','subset',A.subj==subj(i)&A.trialType==1&A.mode==2&A.seqCond==1&A.training>0&A.points>=1&A.jitterMask==0,'style_thickline',...
        %                 'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        %             xlabel ('training');
        %             ylabel('RTabs Interval (%)');
        %         end;
        figure;
        T=tapply(A,{'subj','training'},{A.RTabsInterval,'nanmean','name','RTabsInterval'},'subset',A.trialType==1&A.mode==2&A.seqCond==1&A.training>0&A.points>=1&A.jitterMask==0);
        lineplot(T.training,T.RTabsInterval,'style_thickline',...
            'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        xlabel ('training');
        ylabel('RTabs Interval (%)');
        
        %RTabs deviation
        %         figure;
        %         for i=1:length(unique(A.subj))
        %             subplot(1,length(unique(A.subj)),i);
        %             lineplot(A.training,A.RTabs,'errorfcn','std','subset',A.subj==subj(i)&A.trialType==1&A.mode==2&A.seqCond==1&A.training>0&A.points>=1&A.jitterMask==0,'style_thickline',...
        %                 'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        %             xlabel ('training');
        %             ylabel('RTabs deviation (ms)');
        %         end;
        figure;
        T=tapply(A,{'subj','training'},{A.RTabs,'nanmean','name','RTabs'},'subset',A.trialType==1&A.mode==2&A.seqCond==1&A.training>0&A.points>=1&A.jitterMask==0);
        lineplot(T.training,T.RTabs,'style_thickline',...
            'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        xlabel ('training');
        ylabel('RTabs deviation (ms)');
        
        %errorRate
        %         figure;
        %         for i=1:length(unique(A.subj))
        %             subplot(1,length(unique(A.subj)),i);
        %             lineplot(A.training,A.errorRate,'errorfcn','std','subset',A.subj==subj(i)&A.trialType==1&A.mode==2&A.seqCond==1&A.training>0&A.points>=1&A.jitterMask==0,'style_thickline',...
        %                 'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        %             xlabel ('training');
        %             ylabel('Error (%)');
        %         end;
        figure;
        T=tapply(A,{'subj','training'},{A.errorRate,'nanmean','name','errorRate'},'subset',A.trialType==1&A.mode==2&A.seqCond==1&A.training>0&A.points>=1&A.jitterMask==0);
        lineplot(T.training,T.errorRate,'style_thickline',...
            'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        xlabel ('training');
        ylabel('Error (%)');
        T;
        
    case 'seqRep' % probePost: seq conditions across repetitions against RT deviation (ms)
        colorScheme={[0 1 0],[1 0 0],[0 0 1],[0 0 0]}; %split by seqCond - green=trained,red=trained timing,blue=trained finger order,black=new
        figure;
        for i=1:length(unique(A.subj))
            subplot(1,length(unique(A.subj)),i);
            lineplot(A.trialRep,A.RTabs,'split',A.seqCond,'errorfcn','std','subset',A.subj==subj(i)&A.trialType==1&A.prepost==2&A.points>=1&A.jitterMask==0,...
                'style_thickline','linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
            xlabel('Sequence repetition');
            ylabel('RT deviation (ms)');
        end;
        figure;
        T=tapply(A,{'subj','trialRep','seqCond'},{A.RTabs,'nanmean','name','RTabs'},'subset',A.trialType==1&A.prepost==2&A.points>=1&A.jitterMask==0);
        lineplot(T.trialRep,T.RTabs,'split',T.seqCond,'style_thickline',...
            'linecolor',colorScheme,'errorcolor',colorScheme,'markerfill',colorScheme,'markercolor',colorScheme, 'leg', 'auto');
        xlabel('Sequence repetition');
        ylabel('RT deviation (ms)');
        T;
        
    case 'seqCond' % probePost: seq conditions against RT deviation (ms) after 3th repetition
        colorScheme={[0 1 0],[1 0 0],[0 0 1],[0 0 0]}; %split by seqCond - green=trained,red=trained timing,blue=trained finger order,black=new
        %         figure;
        %         for i=1:length(unique(A.subj))
        %             subplot(1,length(unique(A.subj)),i);
        %             myboxplot(D.seqCond,D.RTabs,'subset',A.subj==subj(i)&A.trialType==1&A.prepost==2&A.points>=1&A.jitterMask==0,'style_tukey',...
        %             'leg','auto');
        %             xlabel('Sequence condition');
        %             ylabel('RT deviation (ms)');
        %         end;
        figure;
        T=tapply(A,{'subj','seqCond'},{A.RTabs,'nanmean','name','RTabs'},'subset',A.trialType==1&A.prepost==2&A.trialRep>=4&A.points>=1&A.jitterMask==0);
        myboxplot(T.seqCond,T.RTabs,'style_block','leg','auto');
        xlabel('Sequence condition');
        ylabel('RT deviation (ms)');
        T;
        
    case 'speedAcc' %speed-accuracy trade-off
        
        A.initTime=A.timing(:,1); %first column is initiation time;
        
        T1=tapply(A,{'subj','prepost'},{A.errorFinger,'nanmean','name','errorFinger'},'subset',A.trialType==1&A.mode==1&A.seqCond&A.jitterMask==0);
        T2=tapply(A,{'subj','prepost'},{A.initTime,'nanmedian','name','initTime'},'subset',A.trialType==1&A.mode==1&A.seqCond&A.jitterMask==0);
        
        for i=1:2
            subplot(1,2,i);
            scatter(T2.initTime(T2.prepost==i,:),T1.errorFinger(T1.prepost==i,:));
            ylim([0 .35]);
            xlim([200 500]);
            
        end;
        A;
        
        %
        
    case 'speedTiming' %initiation time x RTabs
        A.initTime=A.timing(:,1); %first column is initiation time;
        
        T1=tapply(A,{'subj','prepost'},{A.RTabs,'nanmean','name','RTabs'},'subset',A.trialType==1&A.mode==1&A.seqCond&A.jitterMask==0);
        T2=tapply(A,{'subj','prepost'},{A.initTime,'nanmedian','name','initTime'},'subset',A.trialType==1&A.mode==1&A.seqCond&A.jitterMask==0);
        
        for i=1:2
            subplot(1,2,i);
            scatter(T2.initTime(T2.prepost==i,:),T1.errorFinger(T1.prepost==i,:));
            ylim([200 400]);
            xlim([200 500]);
            
            
        end;
        
        %     case 'speedTimingGroup' %initiation time x RTabs
        %         A.initTime=A.timing(:,1); %first column is initiation time;
        %
        %         T1=tapply(A,{'subj','prepost'},{A.RTabs,'nanmean','name','errorFinger'},'subset',A.trialType==1&A.mode==1&A.seqCond&A.jitterMask==0);
        %         T2=tapply(A,{'subj','prepost'},{A.initTime,'nanmedian','name','initTime'},'subset',A.trialType==1&A.mode==1&A.seqCond&A.jitterMask==0);
        %
        %         for i=1:2
        %             subplot(1,2,i);
        %             scatter(T2.initTime(T2.prepost==i,:),T1.errorFinger(T1.prepost==i,:));
        %             ylim([200 400]);
        %             xlim([200 500]);
        %
        %
        %         end;
        
    case 'raster' %plot presses for participants in raster plot - horz lines separate sequences
        
        %%% Raster type graph
        
        %         BNstart=5; %training phase
        %         BNend  =40;
        %
        %         BNstart=29; %test phase
        %         BNend  =40;
        
        
        %         BNstart=01; %pre-training (control trials inc.)
        %         BNend  =02;
        %
        %         BNstart=39; %post-training (control trials inc.)
        %         BNend  =40;
        
        BNstart=47; %fMRI (from memory only)
        BNend=52;
        
        zeroPoint=0; %from first press; 0 is from Go cue; 1 from 1st press.
        
        figure;
        for sbj=1:length(subj)
            spikes=[];
            seq=[];
            spikesTarget=[];
            spikesPress=[];
%             subplot(5,5,sbj);
            
            for i=1:max(D.seqID) %
                %                 S.timing=A.timing(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.trialType==1&A.mode==2,:);
                %                 S.RTabsFingerInterval=A.RTabsFingerInterval(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.trialType==1&A.mode==2,:);
                %                 S.seqID=A.seqID(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0&A.trialType==1&A.mode==2,:);
                %                 S.cueTimePlanned=A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.trialType==1&A.mode==2,:);
                
                S.timing=A.timing(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.points>=1&A.trialType==1&A.mode==2,:);
                S.RTabsFingerInterval=A.RTabsFingerInterval(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.points>=1&A.trialType==1&A.mode==2,:);
                S.seqID=A.seqID(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0&A.points>=1&A.trialType==1&A.mode==2,:);
                S.cueTimePlanned=A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.points>=1&A.trialType==1&A.mode==2,:);
                S.press=A.press(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.points>=1&A.trialType==1&A.mode==2,:);
                
                
                S=dsort(S,S.RTabsFingerInterval);
                spikes=[spikes,S.timing'];
                seq=[seq,S.seqID'];
                spikesTarget=[spikesTarget,S.cueTimePlanned'];
                spikesPress=[spikesPress, S.press'];
            end;
            %
            %                 S.timing=A.timing(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0,:);
            %                 S.RTabsInterval=A.RTabsInterval(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0,:);
            %                 S.seqID=A.seqID(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0,:);
            %                 S.cueTimePlanned=A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0,:);
            %
            %                 S=dsort(S,S.RTabsInterval);
            %                 spikes=[S.timing'];
            %                 seq=[S.seqID'];
            %                 spikesTarget=[S.cueTimePlanned'];
            
            
            for i=1:max(D.seqID)
                if i==1
                    seqID(seq==i)=length(seq(seq==i));
                else
                    seqID(seq==i)=length(seq(seq==i))+length(seq(seq<i));
                end;
            end;
            
            spikes=spikes-repmat(spikes(1,:),size(spikes,1),1)*zeroPoint;
            
            %             nrPressesColor=[0 0 0;...
            %                 1 0 0;...
            %                 0 1 0;...
            %                 0 0 1;...
            %                 1 1 0];
            %              nrPressesColor=[1  0  1;...
            %                 1  1  1;...
            %                 2  1  0;...
            %                 3  1  0;...
            %                 4  0  1];
            
            %  nrPressesColor=[0 0 .8;.2 .2 .6;.4 .4 .4;.6 .6 .4;.8 .8 .2];
            nrPressesColor=[0 0 .8; 1 .2 .6; 0.4 .4 .4; 1 .6 .4; 0 .8 .2]; %new color scheme
            %
            
            %             t=(-200:5500); %old
            %             t=(-200:4000);
            %             t=(-300:4000);
            %             t=(-300:3500);
            %               t=(-300:3000);
            %             t=(-400:3000);
            %             t=(-500:3000);
            t=(-500:4000);
            for i=1:size(spikes,2)
                
                nrPressesColor=[0 0 .8; 1 .2 .6; 0.4 .4 .4; 1 .6 .4; 0 .8 .2]; %reset colour vector at start of every loop
                
                spikesTrial=spikes(:,i);
                spikesTrial(find(isnan(spikesTrial)))=[]; %remove NANs
                spikesTrial = round(spikesTrial);
                
                spikesPressesTrial = spikesPress(:,i);
                
                nrPressesColor = nrPressesColor(spikesPressesTrial,:);
                                
                [lia1,spikeLoc]=ismember(spikesTrial,t);
                %             tColor=(1:length(t));
                
                %Color encodes press number
                tColor=zeros(length(t),3);
                
                
                for k=1:length(spikeLoc)
                    if k==1
                        tColor(1:spikeLoc(k),:)=repmat(nrPressesColor(k,:),length((1:spikeLoc(k))),1);
                    else
                        tColor(spikeLoc(k),:)=repmat(nrPressesColor(k,:),length(spikeLoc(k)),1) ;
                    end;
                end;
                raster=nan(length(t),1); %NaN instead of zeroes for raster plot
                raster(spikeLoc)=1;
                rasterAll(:,i)=raster;
                ntrial=size(spikes,2)-i+1;
                colorCode=t;%tColor;
                %                 scatter(t,raster*ntrial,9,colorCode,'o','filled');
                scatter(t,raster*ntrial,9,tColor,'o','filled');
                
                xlim([min(t) max(t)]);
                hold on;
            end;
            
            drawline(seqID,'dir','horz','color',[0 0 0]);
            
            %% overlap Spike Target
            %             spikesTarget=[A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==1&A.errorFinger==0,:)',...
            %                 A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==2&A.errorFinger==0,:)',...
            %                 A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==3&A.errorFinger==0,:)',...
            %                 A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==4&A.errorFinger==0,:)'];
            %
            
            tColor=zeros(length(t),3);
            spikesTarget=spikesTarget-repmat(spikesTarget(1,:),size(spikesTarget,1),1)*zeroPoint;
            for i=1:size(spikesTarget,2)
                spikesTrial=spikesTarget(:,i);
                spikesTrial(find(isnan(spikesTrial)))=[]; %remove NANs
                [lia1,spikeLoc]=ismember(spikesTrial,t);
                raster=nan(length(t),1); %NaN instead of zeroes for raster plot
                raster(spikeLoc)=1;
                rasterAll(:,i)=raster;
                ntrial=size(spikes,2)-i+1;
                colorCode=[0 0 0];
                scatter(t,raster*ntrial,9,colorCode,'.');
                xlim([min(t) max(t)]);
                xlim([min(t) max(t)]);
                hold on;
            end;
            title(['Subj' num2str(subj(sbj)) ], 'fontSize', 12, 'fontName', 'Calibri');
            title(['Subj' num2str(sbj) ], 'fontSize', 12, 'fontName', 'Calibri');
            
            ylabel('Trial', 'fontSize', 12, 'fontName', 'Calibri')
            xlabel('Time (ms)', 'fontSize', 12, 'fontName', 'Calibri')
            
            set(gca, 'fontSize', 12, 'fontName', 'Calibri')
            
            %to make this figure's x axis match the force example x axis
%             xlim([-1500 3500])
            
        end;
        
    case 'rasterChrono' %as above but plot in chronological order, not ascending deviation
        %%% Raster type graph ordered chronologically to show compression
        
        %         BNstart=5; %training phase
        %         BNend  =40;
        
        %         BNstart=47; %test phase
        %         BNend  =52;
        %
        BNstart=5; %all from memory trials
        BNend  =52;
        
        %         BNstart=01; %pre-training (control trials inc.)
        %         BNend  =02;
        %
        %         BNstart=39; %post-training (control trials inc.)
        %         BNend  =40;
        zeroPoint=0; %from first press
        
        figure;
        for sbj=1:length(subj)
            spikes=[];
            seq=[];
            spikesTarget=[];
            subplot(4,5,sbj);
            
            for i=1:max(D.seqID) %
                S.timing=A.timing(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.trialType==1&A.mode==2,:);
                S.RTabsFingerInterval=A.RTabsFingerInterval(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.trialType==1&A.mode==2,:);
                S.seqID=A.seqID(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0&A.trialType==1&A.mode==2,:);
                S.cueTimePlanned=A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.trialType==1&A.mode==2,:);
                S.blockN=A.TNall(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==i&A.errorFinger==0&A.trialType==1&A.mode==2,:);
                
                S=dsort(S,S.blockN);
                spikes=[spikes,S.timing'];
                seq=[seq,S.seqID'];
                spikesTarget=[spikesTarget,S.cueTimePlanned'];
            end;
            %
            %                 S.timing=A.timing(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0,:);
            %                 S.RTabsInterval=A.RTabsInterval(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0,:);
            %                 S.seqID=A.seqID(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0,:);
            %                 S.cueTimePlanned=A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.errorFinger==0,:);
            %
            %                 S=dsort(S,S.RTabsInterval);
            %                 spikes=[S.timing'];
            %                 seq=[S.seqID'];
            %                 spikesTarget=[S.cueTimePlanned'];
            
            
            for i=1:max(D.seqID)
                if i==1
                    seqID(seq==i)=length(seq(seq==i));
                else
                    seqID(seq==i)=length(seq(seq==i))+length(seq(seq<i));
                end;
            end;
            spikes(spikes<-500)=NaN; %KK added in July 2019; if negative RTs (should not happen when points filter added, but does in a few subjects) set on NaN
            spikes=spikes-repmat(spikes(1,:),size(spikes,1),1)*zeroPoint;
            
            %             nrPressesColor=[0 0 0;...
            %                 1 0 0;...
            %                 0 1 0;...
            %                 0 0 1;...
            %                 1 1 0];
            %              nrPressesColor=[1  0  1;...
            %                 1  1  1;...
            %                 2  1  0;...
            %                 3  1  0;...
            %                 4  0  1];
            
            %  nrPressesColor=[0 0 .8;.2 .2 .6;.4 .4 .4;.6 .6 .4;.8 .8 .2];
            nrPressesColor=[0 0 .8; 1 .2 .6; 0.4 .4 .4; 1 .6 .4; 0 .8 .2]; %new color scheme
            %
            
            %             t=(-200:5500); %old
            %             t=(-200:4000);
            %             t=(-300:4000);
            %             t=(-300:3500);
            %               t=(-300:3000);
            %             t=(-400:3000);
            %             t=(-500:3000);
            t=(-500:4500);
            for i=1:size(spikes,2)
                spikesTrial=spikes(:,i);
                spikesTrial(isnan(spikesTrial))=[]; %remove NANs
                [lia1,spikeLoc]=ismember(spikesTrial,t);
                %             tColor=(1:length(t));
                
                %Color encodes press number
                tColor=zeros(length(t),3);
                
                
                for k=1:length(spikeLoc)
                    if k==1
                        tColor(1:spikeLoc(k),:)=repmat(nrPressesColor(k,:),length((1:spikeLoc(k))),1);
                    else
                        tColor(spikeLoc(k),:)=repmat(nrPressesColor(k,:),length(spikeLoc(k)),1) ;
                    end;
                end;
                raster=nan(length(t),1); %NaN instead of zeroes for raster plot
                
                raster(spikeLoc)=1;
                rasterAll(:,i)=raster;
                ntrial=size(spikes,2)-i+1;
                colorCode=t;%tColor;
                %                 scatter(t,raster*ntrial,9,colorCode,'o','filled');
                %                 scatter(t,raster*ntrial,9,tColor,'o','filled');
                scatter(t,raster*i,9,tColor,'o','filled');
                
                xlim([min(t) max(t)]);
                hold on;
            end;
            
            %             drawline(seqID,'dir','horz','color',[0 0 0]);
            
            %% overlap Spike Target
            %             spikesTarget=[A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==1&A.errorFinger==0,:)',...
            %                 A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==2&A.errorFinger==0,:)',...
            %                 A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==3&A.errorFinger==0,:)',...
            %                 A.cueTimePlanned(A.subj==subj(sbj)&A.BN>=BNstart&A.BN<=BNend&A.seqID==4&A.errorFinger==0,:)'];
            %
            
            tColor=zeros(length(t),3);
            spikesTarget=spikesTarget-repmat(spikesTarget(1,:),size(spikesTarget,1),1)*zeroPoint;
            for i=1:size(spikesTarget,2)
                spikesTrial=spikesTarget(:,i);
                spikesTrial(find(isnan(spikesTrial)))=[]; %remove NANs
                [lia1,spikeLoc]=ismember(spikesTrial,t);
                raster=nan(length(t),1); %NaN instead of zeroes for raster plot
                raster(spikeLoc)=1;
                rasterAll(:,i)=raster;
                ntrial=size(spikes,2)-i+1;
                colorCode=[0 0 0];
                %                 scatter(t,raster*ntrial,9,colorCode,'.');
                scatter(t,raster*i,9,colorCode,'.');
                xlim([min(t) max(t)]);
                xlim([min(t) max(t)]);
                hold on;
            end;
            title(['Subj' num2str(subj(sbj)) ]);
            title(['Subj' num2str(sbj) ]);
        end;
        
        
        
        
        %% Extract data for SPSS
        
        %%%%%% pre-post:
        %RTabsInterval
        %ErrorRate
        
        %%%%%% Reaction time to initiate the sequence:
        
        
        %T;
        
    case 'initiationTimes' %initiation times T test across two timing conditions
        
        %Across participants
        
        Y=tapply(A,{'subj','tempID'},{A.timing(:,1),'nanmedian','name','initTime'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0);
        
        t1=Y.initTime(Y.tempID==1);
        t2=Y.initTime(Y.tempID==2);
        [H P] = ttest(t1,t2,2,'paired')
        
        figure
        lineplot(Y.tempID,Y.initTime,'style_thickline','leg','auto');
        ylabel('Initiation Time (ms)');
        xlabel('Timing')
        ylim([250 600])
        
        
        if P < 0.05 && P >= 0.001
            title('*')
        elseif P < 0.001
            title('**')
        end
        
    case 'initiationTimesBySubj' %as above subj-by-subj
        
        loopCounter = 1;
        figure
        
        for i = subj
            
            disp(['Participant: ' num2str(i)])
            
            Y.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0&A.subj==i);
            Y.tempID = A.tempID(A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0&A.subj==i);
            Y.initTime = A.timing(:,1);
            Y.initTime = Y.initTime(A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0&A.subj==i);
            
            t1=Y.initTime(Y.tempID==1);
            t2=Y.initTime(Y.tempID==2);
            
            [P(loopCounter), H(loopCounter)] = ranksum(t1,t2); %Wilcoxon Ranked test performed as many participants don't have equal trials across conditions
            %             ttest(t1,t2,2,'paired')
            mean1(loopCounter) = mean(Y.initTime(Y.tempID==1));
            mean2(loopCounter) = mean(Y.initTime(Y.tempID==2));
            
            subplot(5,5,loopCounter)
            lineplot(Y.tempID,Y.initTime,'style_thickline','leg','auto');
            ylabel('Initiation Time (ms)');
            xlabel('Timing')
            ylim([250 600])
            
            if P(loopCounter) < 0.05 && P(loopCounter) >= 0.001
                title('*')
            elseif P(loopCounter) < 0.001
                title('**')
            end
            
            loopCounter = loopCounter + 1;
        end
        
        summary = [subj' mean1' mean2' P' H'];
        disp(summary)
        cd('E:\projects\rhys\prepProd2\docs')
        xlswrite('initiationTimesBySubj', summary)
        
    case 'initiationTimesBySubj_allseqs'
            
        Y = tapply(A,{'subj'},{A.timing(:,1),'median','name','initTime'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0);
        
        disp(['mean initiation time = ' num2str(mean(Y.initTime))])
        disp(['SD initiation time = ' num2str(std(Y.initTime))])
        disp(['SE initiation time = ' num2str(std(Y.initTime)/sqrt(length(Y.initTime)))])
        
    case 'movementTimeBySequence' %movement time across the four sequences
        
        %extract final and first press time for all subjs across all seqs
        TLast=tapply(A,{'subj','seqID'},{A.timing(:,5),'nanmean','name','moveTime'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0);
        TFirst=tapply(A,{'subj','seqID'},{A.timing(:,1),'nanmean','name','moveTime'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0);
        
        %subtract first from final press time
        TLast.moveTime = TLast.moveTime - TFirst.moveTime;
        T = TLast;
        
        %plot it
        figure
        myboxplot(T.seqID, T.moveTime, 'xtickoff')
        set(gca,'xticklabel',{'Sequence 1','Sequence 2', 'Sequence 3','Sequence 4'})
        
        ylabel('Total Movement Time (ms)','FontName', 'Calibri', 'FontSize', 12)
        set(gca,'FontSize', 12, 'FontName', 'calibri')
        
        %to plot connecting lines between participant data
        hold on
        dataConcat = [T.moveTime(T.seqID ==1), T.moveTime(T.seqID ==2), T.moveTime(T.seqID ==3), T.moveTime(T.seqID ==4)];
        
        for i = 1:length(subj) %plot lines that connect each participant's data point
%             plot([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'color', [.5, .5, .5], 'lineWidth', 1.2) %solid grey lines
            patchline([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
%             plot([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
            patchline([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
%             plot([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
            patchline([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
        end
        
        cd('G:\projectsBackup\rhys\prepProd2\docs')
        xlswrite('movementTimeBySequence', dataConcat)
        
    case 'movementTimeBySubj'
        
        fileHeader = {'O1T1', 'O1T2', 'O2T1', 'O2T2'};
        allMoveTime = A.timing(:,5) - A.timing(:,1);
        
        for i=subj
            
            moveTime1 = allMoveTime(A.seqID==1&A.subj==i&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
            moveTime2 = allMoveTime(A.seqID==2&A.subj==i&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
            moveTime3 = allMoveTime(A.seqID==3&A.subj==i&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
            moveTime4 = allMoveTime(A.seqID==4&A.subj==i&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
            
            summary = padcat(moveTime1, moveTime2, moveTime3, moveTime4);
            
            cd('G:\projectsBackup\rhys\prepProd2\docs\movementTime')
            xlswrite(['movementTime_' subj_name{i}], summary)
            
        end
        
    case 'absoluteDeviation'
        
        %Across participants
        
        Y=tapply(A,{'subj','tempID'},{A.RTabs(:,1),'nanmean','name','absDev'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0);
        
        t1=Y.absDev(Y.tempID==1);
        t2=Y.absDev(Y.tempID==2);
        [P H] = ttest(t1,t2,2,'paired')
        
        figure
        lineplot(Y.tempID,Y.absDev,'style_thickline','leg','auto');
        ylabel('Absolute Deviation (ms)');
        xlabel('Timing')
        ylim([100 450])
        
        if P < 0.05 && P >= 0.001
            title('*')
        elseif P < 0.001
            title('**')
        end
        
    case 'absoluteDeviationBySubj'
        
        loopCounter = 1;
        
        for i = subj
            
            disp(['Participant: ' num2str(i)])
            
            Y.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0&A.subj==i);
            Y.tempID = A.tempID(A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0&A.subj==i);
            Y.absDev = A.RTabs(:,1);
            Y.absDev = Y.absDev(A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.points>=1&A.jitterMask==0&A.subj==i);
            
            t1=Y.absDev(Y.tempID==1);
            t2=Y.absDev(Y.tempID==2);
            
%             [P(loopCounter), H(loopCounter)] = ranksum(t1,t2); %Wilcoxon Ranked test performed as many participants don't have equal trials across conditions
                        ttest(t1,t2,2,'paired')
            mean1(loopCounter) = mean(Y.absDev(Y.tempID==1));
            mean2(loopCounter) = mean(Y.absDev(Y.tempID==2));
            
            subplot(5,5,loopCounter)
            lineplot(Y.tempID,Y.absDev,'style_thickline','leg','auto','xcat',{'Timing'});
            ylabel('Absolute Deviation (ms)');
            ylim([100 450])
            
            if P(loopCounter) < 0.05 && P(loopCounter) >= 0.001
                title(['S' num2str(i), ' *'])
            elseif P(loopCounter) < 0.001
                title(['S' num2str(i), ' **'])
            end
            
%             title(num2str(i));
            
            loopCounter = loopCounter + 1;
        end
        
        summary = [subj' mean1' mean2' P' H'];
        disp(summary)
        cd('E:\projects\rhys\prepProd2\docs')
        xlswrite('initiationTimesBySubj', summary)
        
    case 'errorAcrossTimings'
        
        %Across participants
        
        trialNr=length(A.errorFinger(A.BN>=47&A.BN<=52&A.subj==i&A.trialType==1&A.mode==2&A.trialType==1));
        Y=tapply(A,{'subj','tempID'},{A.errorFinger,'sum','name','errorFinger'},'subset',A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2);
        Y.errorFingerPerc=Y.errorFinger./trialNr*100;
        
        %         Y=tapply(A,{'subj','tempID'},{A.RTabs(:,1),'nanmean','name','absDev'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.jitterMask==0);
        
        t1=Y.errorFingerPerc(Y.tempID==1);
        t2=Y.errorFingerPerc(Y.tempID==2);
        ttest(t1,t2,2,'paired')
        
        figure
        lineplot(Y.tempID,Y.errorFingerPerc,'style_thickline','leg','auto','xcat',{'Timing'});
        ylabel('Error (%)');
        
    case 'errorAcrossTimingsBySubj'
        
        loopCounter = 1;
        
        for i = subj
            
            disp(['Participant: ' num2str(i)])
            
            for j = 1:4
                trialNr(j)=length(A.errorFinger(A.BN>=47&A.BN<=52&A.subj==i&A.trialType==1&A.mode==2&A.trialType==1&A.seqID==j));
                Y.errorFinger = A.errorFinger(A.BN>=47&A.BN<=52&A.subj==i&A.trialType==1&A.mode==2&A.trialType==1&A.seqID==j);
                Y.errorFingerPerc(j) = Y.errorFinger./trialNr(j)*100;
                
            end
%             Y=tapply(A,{'subj','tempID'},{A.errorFinger,'sum','name','errorFinger'},'subset',A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2);
%             Y.errorFingerPerc=Y.errorFinger./trialNr*100;
            
            Y.subj = A.subj(A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.jitterMask==0&A.subj==i);
            Y.seqID = A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.BN<=52&A.jitterMask==0&A.subj==i);
            Y.errorFingerPerc = Y.errorFinger./trialNr*100;
            
            t1=Y.errorFingerPerc(Y.seqID==1);
            t2=Y.errorFingerPerc(Y.seqID==2);
            t3=Y.errorFingerPerc(Y.seqID==3);
            t4=Y.errorFingerPerc(Y.seqID==4);
            
%             [P(loopCounter), H(loopCounter)] = ranksum(t1,t2); %Wilcoxon Ranked test performed as many participants don't have equal trials across conditions
            %             ttest(t1,t2,2,'paired')
%             mean1(loopCounter) = mean(Y.errorFingerPerc(Y.tempID==1));
%             mean2(loopCounter) = mean(Y.errorFingerPerc(Y.tempID==2));
            
            subplot(5,5,loopCounter)
            lineplot(Y.seqID,Y.errorFingerPerc,'style_thickline','leg','auto','xcat',{'Timing'});
            ylabel('Error (%)');
            ylim([-0.1 0.2])
            
%             if P(loopCounter) < 0.05 && P(loopCounter) >= 0.001
%                 title('*')
%             elseif P(loopCounter) < 0.001
%                 title('**')
%             end
            
            loopCounter = loopCounter + 1;
        end
        
        summary = [subj' mean1' mean2' P' H'];
        disp(summary)
        cd('E:\projects\rhys\prepProd2\docs')
        xlswrite('initiationTimesBySubj', summary)
        
    case 'errorBySequence'
        
        figure
        monoColourScheme={[.9 .6 0], [0 0 0], [0 .6 .5], [.8 .5 .7]};
        
        trialNr=length(A.errorFinger(A.BN>=47&A.BN<=52&A.subj==subj(1)&A.seqID==1&A.trialType==1&A.mode==2&A.trialType==1));
        T=tapply(A,{'seqID','subj'},{A.errorFinger,'sum','name','errorFinger'},'subset',A.BN>=47&A.BN<=52&A.trialType==1);
        T.errorFingerPerc=T.errorFinger./trialNr*100;
        
        myboxplot(T.seqID,T.errorFingerPerc, 'xtickoff')
        set(gca,'xticklabel',{'Sequence 1','Sequence 2', 'Sequence 3','Sequence 4'})
        
        ylim([-5 40])
        ylabel('Error Rate (%)','FontName', 'Calibri', 'FontSize', 12)
        set(gca,'FontSize', 12, 'FontName', 'calibri')
        
        %to plot connecting lines between participant data
        hold on
        dataConcat = [T.errorFingerPerc(T.seqID ==1), T.errorFingerPerc(T.seqID ==2), T.errorFingerPerc(T.seqID ==3), T.errorFingerPerc(T.seqID ==4)];
        
        for i = 1:length(subj) %plot lines that connect each participant's data point
%             plot([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'color', [.5, .5, .5], 'lineWidth', 1.2) %solid grey lines
            patchline([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
%             plot([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
            patchline([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
%             plot([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
            patchline([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
        end
        
        cd('E:\projects\rhys\prepProd2\docs')
        xlswrite('errorBySequence', dataConcat)
        
    case 'errorBySequenceBySubj'
        
%         trialNr=length(A.errorFinger(A.BN>=47&A.BN<=52&A.subj==subj(1)&A.seqID==1&A.trialType==1&A.mode==2&A.trialType==1&A.mode==2));
%         T=tapply(A,{'seqID','subj'},{A.errorFinger,'sum','name','errorFinger'},'subset',A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2);
%         T.errorFingerPerc=T.errorFinger./trialNr*100;
        
        noSeq = 4;
        
        figure
        
        for i=1:length(subj)
            
            subplot(5,5,i)
            
            trialNr=length(A.errorFinger(A.seqID==1&A.trialType==1&A.mode==2&A.BN>=47 & A.BN<=52&A.subj==subj(i)));
            
            T.subj = A.subj(A.seqID==1&A.trialType==1&A.mode==2&A.BN>=47 & A.BN<=52&A.subj==subj(i));
            T.seqID = A.seqID(A.seqID==1&A.trialType==1&A.mode==2&A.BN>=47 & A.BN<=52&A.subj==subj(i));
            
            for j=1:noSeq
                
                T.errorFinger(T.seqID == j) = A.errorFinger(A.seqID==j&A.trialType==1&A.mode==2&A.BN>=47 & A.BN<=52&A.subj==subj(i));
                
                T.errorFingerPerc(T.seqID == j) = T.errorFinger(A.seqID==j)./trialNr*100;
            end
            
            myboxplot(T.seqID,T.errorFingerPerc, 'xtickoff')
            set(gca,'xticklabel',{'Sequence 1','Sequence 2', 'Sequence 3','Sequence 4'})
            
%             ylim([-5 40])
%             ylabel('Error Rate (%)','FontName', 'Calibri', 'FontSize', 25)
%             set(gca,'FontSize', 25, 'FontName', 'calibri')
            
        end
        
    case 'timingBySequence_prodlengthCorrected'
        
        %get the mean production time across sequences within participants
        seqLength = tapply(A,{'subj'}, {A.timing(:,5),'nanmean','name','prodLength'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        lc = 1; %loopCounter
        
        for i = subj
            for j = 1:max(A.seqID(A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0)) %for each of the sequences in the fMRI blocks
                
                %get the IPI times for each production within sequences (x4)
                timings = A.timingInterval(A.subj==i&A.seqID==j&A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0,:);
                
                %get the IPI target times for each production within sequences
                targets = A.targetIntervalDiff(A.subj==i&A.seqID==j&A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0,:);
                
                %use mean production time to calculate percentage lengths of produced intervals
                percTimings = bsxfun(@rdivide,timings,seqLength.prodLength(seqLength.subj==i));
                
                
                %use mean production time to calculate target interval percentages
                percTargets = bsxfun(@rdivide,targets,seqLength.prodLength(seqLength.subj==i));
                
                %calculate error as production time - target time, then show as absolute
                percIntervalDev = abs(percTimings - percTargets) * 100;
                totPercDev = sum(percIntervalDev,2);
                meanPercDev(lc,j) = mean(totPercDev);
                
            end
            lc = lc + 1; %up loop counter
        end
        
        %put our interval deviation numbers into a struct to be compatible with myboxplot function
        T1.meanPercDev = meanPercDev(:,1); T1.subj = seqLength.subj; T1.seqID = ones(length(T1.subj),1) *1; %mean perc, subj ID, and seq ID
        T2.meanPercDev = meanPercDev(:,2); T2.subj = seqLength.subj; T2.seqID = ones(length(T2.subj),1) *2;
        T3.meanPercDev = meanPercDev(:,3); T3.subj = seqLength.subj; T3.seqID = ones(length(T3.subj),1) *3;
        T4.meanPercDev = meanPercDev(:,4); T4.subj = seqLength.subj; T4.seqID = ones(length(T4.subj),1) *4;
        
        T = addstruct(T1, T2); T = addstruct(T, T3); T = addstruct(T, T4);
        
        %produce figure
        figure
        myboxplot(T.seqID, T.meanPercDev, 'xtickoff')
        set(gca,'xticklabel',{'Sequence 1','Sequence 2', 'Sequence 3','Sequence 4'})
        
        ylabel('Timing Error (% whole sequence)','FontName', 'Calibri', 'FontSize', 12)
        set(gca,'FontSize', 12, 'FontName', 'calibri')
        
        %to plot connecting lines between participant data
        hold on
        dataConcat = meanPercDev;
        
        for i = 1:length(subj) %plot lines that connect each participant's data point
%             plot([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'color', [.5, .5, .5], 'lineWidth', 1.2) %solid grey lines
            patchline([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
%             plot([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
            patchline([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
%             plot([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
            patchline([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
        end
        
        cd('G:\projectsBackup\rhys\prepProd2\docs')
        xlswrite('timingBySequence_prodlengthCorrected', dataConcat)
        
    case 'timingBySequence'
        
        T=[];
        
        %collect all produced intervals, followed by all target intervals.
        T1=tapply(A,{'subj','seqID'},{A.timingInterval1,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        T2=tapply(A,{'subj','seqID'},{A.timingInterval2,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        T3=tapply(A,{'subj','seqID'},{A.timingInterval3,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        T4=tapply(A,{'subj','seqID'},{A.timingInterval4,'nanmean','name','timing'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        
        K1 = tapply(A,{'subj','seqID'},{A.targetIntervalDiff1,'nanmean','name','targetIntervalDiff'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        K2 = tapply(A,{'subj','seqID'},{A.targetIntervalDiff2,'nanmean','name','targetIntervalDiff'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        K3 = tapply(A,{'subj','seqID'},{A.targetIntervalDiff3,'nanmean','name','targetIntervalDiff'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        K4 = tapply(A,{'subj','seqID'},{A.targetIntervalDiff4,'nanmean','name','targetIntervalDiff'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.points>=1&A.jitterMask==0);
        
        %concatenate...
        timingsConcatenated = [T1.timing, T2.timing, T3.timing, T4.timing]; %produce arrays containing all timings and targets
        targetsConcatenated = [K1.targetIntervalDiff, K2.targetIntervalDiff, K3.targetIntervalDiff, K4.targetIntervalDiff];
        
        %calculate deviation
        timingsConcatenatedError = timingsConcatenated - targetsConcatenated;
        
        %mean across interval positions
        T.timingError = mean(timingsConcatenatedError, 2);
        T.subj = T1.subj;
        T.seqID = T1.seqID;
        
        %produce figure
        figure
        myboxplot([T.seqID], T.timingError, 'xtickoff')
        set(gca,'xticklabel',{'Sequence 1','Sequence 2', 'Sequence 3','Sequence 4'})
        
        ylabel('Timing Error (ms)','FontName', 'Calibri', 'FontSize', 12)
        set(gca,'FontSize', 12, 'FontName', 'calibri')
        
        %to plot connecting lines between participant data
        hold on
        dataConcat = [T.timingError(T.seqID ==1), T.timingError(T.seqID ==2), T.timingError(T.seqID ==3), T.timingError(T.seqID ==4)];
        
        for i = 1:length(subj) %plot lines that connect each participant's data point
%             plot([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'color', [.5, .5, .5], 'lineWidth', 1.2) %solid grey lines
            patchline([1, 2], [dataConcat(i,1), dataConcat(i,2)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
%             plot([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
            patchline([2, 3], [dataConcat(i,2), dataConcat(i,3)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
%             plot([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'color', [.5, .5, .5], 'lineWidth', 1.2)
            patchline([3, 4], [dataConcat(i,3), dataConcat(i,4)], 'edgecolor',[0, 0, 0], 'edgealpha', 0.2, 'lineWidth', 1.2)
        end
        
        cd('E:\projects\rhys\prepProd2\docs')
        xlswrite('timingBySequence', dataConcat)
        
    case 'timingBySequenceBySubj'
        
        T.timingInterval1 = A.timingInterval1(A.seqID==1&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
        T.timingInterval2 = A.timingInterval2(A.seqID==1&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
        T.timingInterval3 = A.timingInterval3(A.seqID==1&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
        T.timingInterval4 = A.timingInterval4(A.seqID==1&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
        
        T.timingsConcatenated = [T.timingInterval1, T.timingInterval2, T.timingInterval3, T.timingInterval4];
        
        T.targetIntervalDiff1 = A.targetIntervalDiff1(A.seqID==1&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
        T.targetIntervalDiff2 = A.targetIntervalDiff2(A.seqID==1&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
        T.targetIntervalDiff3 = A.targetIntervalDiff3(A.seqID==1&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
        T.targetIntervalDiff4 = A.targetIntervalDiff4(A.seqID==1&A.BN>=47&A.BN<=52&A.trialType==1&A.mode==2&A.trialType==1&A.points>=1);
        
        T.targetsConcatenated = [T.targetIntervalDiff1, T.targetIntervalDiff2, T.targetIntervalDiff3, T.targetIntervalDiff4];
        
        [~, T.longestIntervalIndex] = max(T.targetsConcatenated(:,1));
        
        T.longestInterval = T.timingsConcatenated(:, T.longestIntervalIndex);
        T.longestInterval = repmat(T.longestInterval, 1, 4);
        
        T.adjustedTimingsConcatenated = T.timingsConcatenated./ T.longestInterval;
        
    case 'RT_TrainingXfMRI'
        
        %Across participants
        
        training=tapply(A,{'subj'},{A.timing(:,1),'nanmean','name','initTime'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
        fMRI=tapply(A,{'subj'},{A.timing(:,1),'nanmean','name','initTime'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.BN<=520&A.points>=1&A.jitterMask==0);
        
        a=training.initTime;
        b=fMRI.initTime;
        [t p] = ttest(a,b,2,'paired')
        
        trainMean = mean(a);
        trainStd = std(a);
        fMRIMean = mean(b);
        fMRIStd = std(b);
        
        disp(['Training Mean = ' num2str(trainMean)])
        disp(['Standard deviation = ' num2str(trainStd)])
        disp(['fMRI Mean = ' num2str(fMRIMean)])
        disp(['Standard deviation = ' num2str(fMRIStd)])
        
%         figure
%         lineplot(Y.tempID,Y.initTime,'style_thickline','leg','auto');
%         ylabel('Initiation Time (ms)');
%         xlabel('Timing')
%         ylim([250 600])
%         
%         
%         if P < 0.05 && P >= 0.001
%             title('*')
%         elseif P < 0.001
%             title('**')
%         end
        
    case 'deviation_TrainingXfMRI'
        
        %Across participants
        
        training=tapply(A,{'subj'},{A.RTabs(:,1),'nanmean','name','absDev'},'subset',A.trialType==1&A.mode==2&A.BN>=29&A.BN<=40&A.points>=1&A.jitterMask==0);
        fMRI=tapply(A,{'subj'},{A.RTabs(:,1),'nanmean','name','absDev'},'subset',A.trialType==1&A.mode==2&A.BN>=47&A.BN<=520&A.points>=1&A.jitterMask==0);
        
        
        a=training.absDev;
        b=fMRI.absDev;
        [t p] = ttest(a,b,2,'paired')
        
        trainMean = mean(a);
        trainStd = std(a);
        fMRIMean = mean(b);
        fMRIStd = std(b);
        
        disp(['Training Mean = ' num2str(trainMean)])
        disp(['Standard deviation = ' num2str(trainStd)])
        disp(['fMRI Mean = ' num2str(fMRIMean)])
        disp(['Standard deviation = ' num2str(fMRIStd)])
        
        figure
%         lineplot(Y.tempID,Y.absDev,'style_thickline','leg','auto');
%         ylabel('Absolute Deviation (ms)');
%         xlabel('Timing')
%         ylim([100 450])
%         
%         if P < 0.05 && P >= 0.001
%             title('*')
%         elseif P < 0.001
%             title('**')
%         end
        
end