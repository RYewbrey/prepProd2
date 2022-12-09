function prepProd2dynamics_ana_run_RY

%%% Add before starting any scripts (comment in when pasting into command line):
% addpath(genpath('E:\projects\rhys\prepProd2\matlab\dynamics'));
% addpath(genpath('E:\projects\rhys\prepProd2\data\behavioural\forces')); %path to behavioural data containing force traces
% addpath(genpath('D:\projects\toolboxes\userfun')); %joern's util tools (open source)

subjName={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
    's11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23',...
    's24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34','s35','s36',...
    's37','s38','s39','s40','s41','s42','s43','s44','s45','s46','s47','s48','s49',...
    's50','s51','s52','s53','s54','s55','s56','s57','s58','s59','s60'}; %% chronological without missing subject numbers, for later vector references

baseDir = 'E:\projects\rhys\prepProd2\data\behavioural\forces\group';
cd(baseDir)

subj=[3,5,6,7,9,10,13,16,17,18,20,21,22,25,26,31,32,34,36,38,39,40,41,42]; %meet both criteria (interaction & error rate)
%8
T = struct([]);

for i = subj
    load(sprintf('%s_overlap_data.mat', subjName{i}))
    
    if i == subj(1)
        T = R;
    else
        T = addstruct(T, R);
    end
    
end


N = tapply(T,{'subj'},{'areaIntAvg','mean','name','areaIntAvg'},'subset',T.errorFinger == 0 & T.mode == 2 & T.BN >= 47); %standard analysis
% N = tapply(T,{'seqID','subj'},{'areaIntAvg','mean','name','areaIntAvg'},'subset',T.errorFinger == 0 & T.mode == 2 & T.BN >= 47); %standard analysis
% N = tapply(T,{'cond','subj','TNminiblock'},{'areaIntAvg','mean','name','areaIntAvg'},'subset',T.correctOrder == 1 & T.BN > 69); %split by number of exposures to each sequence

%%single measure per participants
results(:,1) = subj;
results(:,2) = N.areaIntAvg;

%%split by seqID
% results(:,1) = subj;
% results(:,2) = N.areaIntAvg(R.cond == 1);
% results(:,3) = N.areaIntAvg(R.cond == 2);
% results(:,4) = N.areaIntAvg(R.cond == 3);
% results(:,5) = N.areaIntAvg(R.cond == 4);

%%split by number of exposures
% results(:,1) = subj;
% loopCounter = 2;
% for i=1:length(unique(N.cond))
%     for j = 1:length(unique(N.TNminiblock))
%         results(:,loopCounter) = N.areaIntAvg(R.cond == i & N.TNminiblock == j);
%         loopCounter = loopCounter + 1;
%     end
% end


xlswrite ('overlap_data', results)
% xlswrite ('overlap_data_splitbyseqID', results)
% xlswrite ('overlap_data_splitbyreps', results)