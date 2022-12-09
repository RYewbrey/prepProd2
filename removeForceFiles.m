clearvars

subj_name={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
    's11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23',...
    's24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34','s35','s36',...
    's37','s38','s39','s40','s41','s42','s43','s44','s45','s46','s47','s48','s49',...
    's50','s51','s52','s53','s54','s55','s56','s57','s58','s59','s60'}; %% chronological without missing subject numbers, for later vector references

% subj=[3,5,6,7,9,10,13,16,17,18,20,21,22,25,26,31,32,34,36,38,39,40,41,42];
subj=2;

baseDir = 'E:\projects\rhys\prepProd2\data\behavioural\';
newDir = 'E:\projects\rhys\prepProd2\data\behavioural\forces\';

minBN = 0;
blockNo = 52;

for i = subj %loop through participants

    for j=minBN:blockNo %loop through blocks
        
        cd(fullfile(baseDir, subj_name{i}))
        fname=sprintf('exp_BN%02d.mat',(j)); %identify and load each block
        
        fileIn=fullfile(baseDir,subj_name{i},fname);
        load(fileIn);
        
        E = B; %E is to be a complete copy of B. The force files from B are then to be deleted (to help load times on ana scripts)
        
        if iscell(B.forces) %for each force-related variable, clear the variable loaded in 'B'
            B.forces = [];
        end
        
        if iscell(B.timeStamps)
            B.timeStamps = [];
        end
        
        if iscell(B.thresholdedForces)
            B.thresholdedForces = [];
        end
        
        if iscell(B.thresholdedSmoothedForces)
            B.thresholdedSmoothedForces = [];
        end
        
        save(fname, 'B') %save the new, compact B in the exact same place it was loaded
        
        cd(newDir)
        
        if ~exist([newDir subj_name{i}], 'dir') %produce directory for forces data, if it doesn't exist
            mkdir([newDir subj_name{i}])
        end
        save([newDir subj_name{i} '\' fname], 'E') %save into subject folder in forces subfolder
        
    end
end