function prepProd2_imanaBatch_RY

%vector containing subject IDs in ascending order, add more as needed
subjName={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
    's11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23',...
    's24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34','s35','s36',...
    's37','s38','s39','s40','s41','s42','s43','s44','s45','s46','s47','s48','s49',...
    's50','s51','s52','s53','s54','s55','s56','s57','s58','s59','s60'}; %% chronological without missing subject numbers, for later vector references

%subjects to be analysed
subj=[3,5,6,7,9,10,13,16,17,18,20,21,22,25,26,31,32,34,36,38,39,40,41,42]; %Ps that reached performance threshold
% subj=[3 4 5 6 7 9 10 11 12 13 15 16 17 18 20 21 22 25 26 31 32 33 34 36 37 38 39 40 41 42]; %all Ps that were scanned regardless of performance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Analysing: ' num2str(subj)])
disp(['Subj N = ' num2str(length(subj))])

% parfor i=1:numel(subj)
%     disp('Subj#'); disp(subjName{subj(i)});
%     %
%     %     prepProd2_imana_RY('make_nii',subj(i));
%     %
%     %     prepProd2_imana_RY('set_AC',subj(i)); %%TO BE DONE MANUALLY
%     %
%     %     prepProd2_imana_RY('slice_timing',subj(i));
%     %     prepProd2_imana_RY('realign_unwarp',subj(i));
%     %     prepProd2_imana_RY('meanepi',subj(i));
%     %
%     %     prepProd2_imana_RY('coreg',subj(i)); %%TO BE DONE MANUALLY
%     %
%     %     prepProd2_imana_RY('glm_set',subj(i));
%     %     prepProd2_imana_RY('glm_estimate',subj(i));
%     %     prepProd2_imana_RY('glm_contrast',subj(i));
%     %         prepProd2_imana_RY('con_smooth', subj(i));
%     %
%     %     prepProd2_imana_RY('segment',subj(i));
%     %     prepProd2_imana_RY('make_mask',subj(i));
%     %     prepProd2_imana_RY('MVA_search',subj(i));
%     %     prepProd2_imana_RY('MVA_do_overallMov',subj(i));
%     %     prepProd2_imana_RY('MVA_do_overallPrep',subj(i));
%     %     prepProd2_imana_RY('MVA_do_spatOneout_Mov',subj(i));
%     %     prepProd2_imana_RY('MVA_do_spatOneout_Prep',subj(i));
%     %     prepProd2_imana_RY('MVA_do_tempOneout_Mov',subj(i));
%     %     prepProd2_imana_RY('MVA_do_tempOneout_Prep',subj(i));
%     %     prepProd2_imana_RY('MVA_do_Int_Mov',subj(i));
%     %     prepProd2_imana_RY('MVA_do_Int_Prep',subj(i));
%     %     prepProd2_imana_RY('MVA_do_Int_Mov_2x2',subj(i));
%     %     prepProd2_imana_RY('MVA_do_Int_Prep_2x2',subj(i));
%     %     prepProd2_imana_RY('MVA_zValue',subj(i));
%     %     prepProd2_imana_RY('MVA_zValue_oneOut',subj(i));
%     %     prepProd2_imana_RY('MVA_smooth',subj(i));
%     %     prepProd2_imana_RY('MNI_normalization',subj(i));
%     %     prepProd2_imana_RY('MNI_normalization_perc&con',subj(i));
%     
% end

%%%%%%%%%%%%% WHOLE BRAIN GROUP %%%%%%%%%%%%
%     prepProd2_imana_RY('group_avg')
%     prepProd2_imana_RY('glm_contrastGroup');
%     prepProd2_imana_RY('glm_contrastEstimate');
%     prepProd2_imana_RY('MVA_group');
%     prepProd2_imana_RY('MVA_estimate');



%%%%%%%%%%%%%% SUBCORTICAL ANALYSIS %%%%%%%%%%%%%%
% parfor i=1:numel(subj)
%     disp(['RSA Subcortical Analysis: ' subjName{subj(i)}])
% %     prepProd2_imana_RY('subcortical_make_nii',subj(i))
% %     prepProd2_imana_RY('subcortical_make_structs',subj(i))
% %     prepProd2_imana_RY('subcortical_reslice_structs',subj(i))
% %     prepProd2_imana_RY('subcortical_funcmask_structs',subj(i))
% %     prepProd2_imana_RY('subcortical_make_ROIs',subj(i))
% %     prepProd2_imana_RY('subcortical_make_search',subj(i))
%     
% %     prepProd2_imana_RY('subcortical_run_search_RSA',subj(i),1)
%     prepProd2_imana_RY('subcortical_run_overallMov_LDA',subj(i),1)
%     prepProd2_imana_RY('subcortical_run_overallPrep_LDA',subj(i),1)
% %     prepProd2_imana_RY('subcortical_run_spatMov_LDA',subj(i))
% %     prepProd2_imana_RY('subcortical_run_spatPrep_LDA',subj(i))
% %     prepProd2_imana_RY('subcortical_run_tempMov_LDA',subj(i))
% %     prepProd2_imana_RY('subcortical_run_tempPrep_LDA',subj(i))
% %     prepProd2_imana_RY('subcortical_run_intMov_LDA',subj(i))
% %     prepProd2_imana_RY('subcortical_run_intPrep_LDA',subj(i))
%     prepProd2_imana_RY('subcortical_zValue_LDA',subj(i))
% %     prepProd2_imana_RY('subcortical_wholebrain_percent_signal', subj(i))
% %     prepProd2_imana_RY('subcortical_segment_percent_signal',subj(i))
% %     prepProd2_imana_RY('subcortical_segment_contrasts',subj(i))
%     prepProd2_imana_RY('subcortical_smooth',subj(i))
% %     prepProd2_imana_RY('subcortical_calc_dissimilarity_maps',subj(i))
% %     prepProd2_imana_RY('subcortical_normalise_RSA',subj(i))
%     prepProd2_imana_RY('subcortical_normalise_LDA',subj(i))
% %     prepProd2_imana_RY('subcortical_normalise_conperc',subj(i))
% 
% 
% 
% %     
% %     prepProd2_imana_RY('subcortical_voxel_counts',subj(i))
% end

% disp('Pre-whitening subcortical...')
% prepProd2_imana_RY('subcortical_preWhiten',1)
% 
% disp('calculating subcortical RDMs and formatting LDA...')
% prepProd2_imana_RY('subcortical_calculate')

%Sub Group functions
% prepProd2_imana_RY('subcortical_group_avg_RSA')
% prepProd2_imana_RY('subcortical_group_avg_LDA')
% prepProd2_imana_RY('subcortical_group_avg_perc')
% prepProd2_imana_RY('subcortical_group_randomeffects_RSA')
% prepProd2_imana_RY('subcortical_group_randomeffects_LDA')
% prepProd2_imana_RY('subcortical_group_randomeffects_con')
% prepProd2_imana_RY('subcortical_group_estimate_RSA')
% prepProd2_imana_RY('subcortical_group_estimate_LDA')
% prepProd2_imana_RY('subcortical_group_estimate_con')
% prepProd2_imana_RY('subcortical_normalise_anat_masks',subj(i))


parfor i=1:numel(subj)
    disp(['Cerebellum Analysis: ' subjName{subj(i)}])
    
%     prepProd2_imana_RY('cerebellum_make_nii',subj(i))
%     prepProd2_imana_RY('cerebellum_suit_normalise',subj(i))
%     prepProd2_imana_RY('cerebellum_make_mask',subj(i))
%     if i == 1
%         prepProd2_imana_RY('cerebellum_make_structs')
%     end
%     prepProd2_imana_RY('cerebellum_structs_to_anat', subj(i))
%     prepProd2_imana_RY('cerebellum_reslice_structs', subj(i))
%     prepProd2_imana_RY('cerebellum_funcmask_structs', subj(i))
%     prepProd2_imana_RY('cerebellum_make_ROIs', subj(i))
%     prepProd2_imana_RY('')
%     prepProd2_imana_RY('')
%     prepProd2_imana_RY('')

%     prepProd2_imana_RY('cerebellum_make_search',subj(i))
%     prepProd2_imana_RY('cerebellum_run_search_RSA',subj(i),1) %varargin 2 = 1 for blueBear
%     prepProd2_imana_RY('cerebellum_run_spatMov_LDA',subj(i))
%     prepProd2_imana_RY('cerebellum_run_spatPrep_LDA',subj(i))
%     prepProd2_imana_RY('cerebellum_run_tempMov_LDA',subj(i))
%     prepProd2_imana_RY('cerebellum_run_tempPrep_LDA',subj(i))
    prepProd2_imana_RY('cerebellum_run_intMov_LDA',subj(i))
    prepProd2_imana_RY('cerebellum_run_intPrep_LDA',subj(i))
    prepProd2_imana_RY('cerebellum_zValue_LDA',subj(i))
%     prepProd2_imana_RY('cerebellum_reslice_contrast',subj(i))
    prepProd2_imana_RY('cerebellum_smooth',subj(i))
%     prepProd2_imana_RY('cerebellum_normalise_contrast',subj(i))
%     prepProd2_imana_RY('cerebellum_calc_dissimilarity_maps',subj(i))
%     prepProd2_imana_RY('cerebellum_normalise_RSA',subj(i))
    prepProd2_imana_RY('cerebellum_normalise_LDA',subj(i))
end

% disp('Pre-whitening cerebellum...')
% prepProd2_imana_RY('cerebellum_preWhiten',1)
% 
% disp('calculating cerebellum RDMs and formatting LDA...')
% prepProd2_imana_RY('cerebellum_calculate')

%CB Group functions
% prepProd2_imana_RY('cerebellum_group_avg_RSA')
prepProd2_imana_RY('cerebellum_group_avg_LDA')
% prepProd2_imana_RY('cerebellum_group_avg_psc')
prepProd2_imana_RY('cerebellum_group_randomeffects')
prepProd2_imana_RY('cerebellum_group_estimate')

% parfor i=1:numel(subj)
%     disp(['RSA Cortical Analysis: ' subjName{subj(i)}])
%     
% %     prepProd2_imana_RY('cortical_make_search',subj(i))
% %     prepProd2_imana_RY('cortical_run_search',subj(i),1) %set varargin 2 to 1 for blueBear
% %     prepProd2_imana_RY('cortical_run_search_integrated',subj(i),1) %set varargin 2 to 1 for blueBear
% %     prepProd2_imana_RY('cortical_smooth',subj(i))
% %     prepProd2_imana_RY('cortical_calc_dissimilarity_maps',subj(i))
% %     prepProd2_imana_RY('cortical_calc_dissimilarity_maps_integrated',subj(i))
% %     prepProd2_imana_RY('cortical_normalise',subj(i))
% 
% end



% prepProd2_imana_RY('')