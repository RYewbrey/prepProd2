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

nconditions = 4*2; %4 sequences * prepProd
nrruns = 6; %nruns (referred to as nSessions in RSA toolbox)

runBSL=[0 0 0 0 0 0]; %rest baseline to attatch at the end of the vectors
prep      =[0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0];
prep=[repmat(prep,1,nrruns) runBSL];

prod      =[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0];
prod=[repmat(prod,1,nrruns) runBSL];

[~,prepCol] = find(prep>0); prepCol = reshape(prepCol,4,[]);
[~,prodCol] = find(prod>0); prodCol = reshape(prodCol,4,[]);
col = [prepCol; prodCol]';

for i=1:nrruns%for runs
    for j=1:nconditions%for conditions
        if col(i,j)<10
            fill = '000';
        elseif col(i,j)<100
            fill = '00';
        else
            fill = '0';
        end
        betas(i,j).identifier = sprintf('beta_%s%d.nii',fill,col(i,j)); %runs x conditions matrix
    end%for runs
end%for conditions
end%function
