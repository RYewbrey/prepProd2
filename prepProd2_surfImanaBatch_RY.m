function prepProd2_surfImanaBatch_RY

%vector containing subject IDs in ascending order, add more as needed
subjName={'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10',...
    's11','s12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23',...
    's24','s25','s26','s27','s28','s29','s30','s31','s32','s33','s34','s35','s36',...
    's37','s38','s39','s40','s41','s42'}; %chonological without missing subject numbers, for later vector references

%subjects to be analysed
% subj= [3 5 6 7 9 10 13 15 16 17 18 20 21 22 25 26 31 32 33 34 35 36 37 38 39 40 41 42]; %all
subj= [3 5 6 7 9 10 13 16 17 18 20 21 22 25 26 31 32 34 36 38 39 40 41 42]; %subjs for analysis & writeup

% subj= [20];

disp(['Group to be analysed: ' num2str(subj)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PreProcessing:

% Perform surface reconstruction first using freesurfer; see freesurfer_notes file on office mac
% Open matlab through the terminal, after initialising freesurfer and
% CDing to $subjects_dir

for i=subj
    disp('Subj'); disp(subjName{i});
    
%     prepProd2_surfImana('surf_xhemireg',i)
%     prepProd2_surfImana('surf_map_ico',i,2) %varargin{2} refers to atlas (2=x)
%     prepProd2_surfImana('surf_make_caret',i,2)
% 
%     prepProd2_surfImana('surf_map_con',i)
%     for h=1:2 %define_search required to loop across both hemispheres with varargin{2}
%     prepProd2_surfImana('surf_define_search',i,h) %requires maskbrain.nii from volumetric analysis
%     end
%     
%     prepProd2_surfImana('MVA_search_surf',i,2) %VOLUMETRIC, not used in later MVA
      
%     prepProd2_surfImana('MVA_do_overall_Mov_surf',i)
%     prepProd2_surfImana('MVA_do_overall_Prep_surf',i)
%     prepProd2_surfImana('MVA_do_spatOneout_Mov_surf',i)
%     prepProd2_surfImana('MVA_do_spatOneout_Prep_surf',i)
%     prepProd2_surfImana('MVA_do_tempOneout_Mov_surf',i)
%     prepProd2_surfImana('MVA_do_tempOneout_Prep_surf',i)
%     prepProd2_surfImana('MVA_do_Int_Mov_surf',i)
%     prepProd2_surfImana('MVA_do_Int_Prep_surf',i)
      
%     prepProd2_surfImana('surf_map_acc', i)
    
end
%     prepProd2_surfImana('surf_avrgcoord')
%     prepProd2_surfImana('surf_makeGroup')
%     prepProd2_surfImana('surf_zacc')
%     prepProd2_surfImana('surf_groupSmooth')
%     prepProd2_surfImana('surf_groupTtest')
%     prepProd2_surfImana('surf_groupTtest_unsmoothed') %produce group statistical maps on unsmoothed data (prompted by Eva Berlot)
    
        %RoI
for i=subj
%     prepProd2_surfImana('ROI_define',i)
%     prepProd2_surfImana('percent_signal',i)
end
% prepProd2_surfImana('ROI_make_BAs_paint')
% for i=subj
%     prepProd2_surfImana('ROI_BA_define',i)
% end
% prepProd2_surfImana('surf_stat_mask')
% prepProd2_surfImana('surf_stat_list')
% prepProd2_surfImana('surf_ROI_MVA')
% prepProd2_surfImana('surf_ROI_stat_RY',2,2)
% 
prepProd2_surfImana('result_crosssection_RY')
% prepProd2_surfImana('result_crosssectionContrasts_RY')