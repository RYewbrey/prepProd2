function prepProd2_imanaBatch_RY

%vector containing subject IDs in ascending order, add more as needed
subj_name={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
    's11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23',...
    's24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34','s35','s36',...
    's37','s38','s39','s40','s41','s42','s43','s44','s45','s46','s47','s48','s49',...
    's50','s51','s52','s53','s54','s55','s56','s57','s58','s59','s60'}; %% chronological without missing subject numbers, for later vector references

%subjects to be analysed
subj=[2,3,4,5,6,8,9,10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Preprocessing:

%

for i=subj
    disp('Subj#'); disp(subj_name{i});
    
%     prepProd2_imana_RY('make_nii',i);

%     prepProd2_imana_RY('set_AC',i); %%TO BE DONE MANUALLY

    prepProd2_imana_RY('slice_timing',i);
    prepProd2_imana_RY('realign_unwarp',i);
    prepProd2_imana_RY('meanepi',i);

%     prepProd2_imana_RY('coreg',i); %%TO BE DONE MANUALLY
    
%     prepProd2_imana_RY('glm_set',i);
%     prepProd2_imana_RY('glm_estimate',i);
%     prepProd2_imana_RY('glm_contrast',i);
%     prepProd2_imana_RY('con_smooth', i);

%     prepProd2_imana_RY('segment',i);
%     prepProd2_imana_RY('make_mask',i);
%     prepProd2_imana_RY('MVA_search',i);
%     prepProd2_imana_RY('MVA_do_overallMov',i);
%     prepProd2_imana_RY('MVA_do_overallPrep',i);
%     prepProd2_imana_RY('MVA_do_spatOneout_Mov',i);
%     prepProd2_imana_RY('MVA_do_spatOneout_Prep',i);
%     prepProd2_imana_RY('MVA_do_tempOneout_Mov',i);
%     prepProd2_imana_RY('MVA_do_tempOneout_Prep',i);
%     prepProd2_imana_RY('MVA_do_Int_Mov',i);
%     prepProd2_imana_RY('MVA_do_Int_Prep',i);
%     prepProd2_imana_RY('MVA_group');
%     prepProd2_imana_RY('MVA_estimate');
%     prepProd2_imana_RY('MVA_zValue',i);
%     prepProd2_imana_RY('MVA_zValue_oneOut',i);
%     prepProd2_imana_RY('MVA_smooth',i);
%     prepProd2_imana_RY('MNI_normalization',i);

end

%     prepProd2_imana_RY('glm_contrastGroup');
%     prepProd2_imana_RY('glm_contrastEstimate');