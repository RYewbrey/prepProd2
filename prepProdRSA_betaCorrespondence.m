function betas = prepProdRSA_betaCorrespondence()
%
%  betaCorrespondence.m is a simple function which should combine
%  three things: preBeta:	a string which is at the start of each file
%  containing a beta image, betas:	a struct indexed by (session,
%  condition) containing a sting unique to each beta image, postBeta:	a
%  string which is at the end of each file containing a beta image, not
%  containing the file .suffix
% 
%  use "[[subjectName]]" as a placeholder for the subject's name as found
%  in userOptions.subjectNames if necessary For example, in an experment
%  where the data from subject1 (subject1 name)  is saved in the format:
%  subject1Name_session1_condition1_experiment1.img and similarly for the
%  other conditions, one could use this function to define a general
%  mapping from experimental conditions to the path where the brain
%  responses are stored. If the paths are defined for a general subject,
%  the term [[subjectName]] would be iteratively replaced by the subject
%  names as defined by userOptions.subjectNames.
% 
%  note that this function could be replaced by an explicit mapping from
%  experimental conditions and sessions to data paths.
% 
%  Cai Wingfield 1-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

% nconditions = 8;
% 
% for condI=1:nconditions
%     if condI<10
%         fill = '000';
%     else
%         fill = '00';
%     end
%     betas(1,condI).identifier = sprintf('spmT_%s%d.img',fill,condI);
% end

betas(1,1).identifier = 'beta_0001.nii';
betas(1,2).identifier = 'beta_0003.nii';
betas(1,3).identifier = 'beta_0005.nii';
betas(1,4).identifier = 'beta_0007.nii';
betas(1,5).identifier = 'beta_0009.nii';
betas(1,6).identifier = 'beta_0011.nii';
betas(1,7).identifier = 'beta_0013.nii';
betas(1,8).identifier = 'beta_0015.nii';


end%function
