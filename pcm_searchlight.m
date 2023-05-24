function R=pcm_searchlight(S,P,PAll,out,mva_func,SPMfiles,varargin)
% function pcm_searchlight(S,P,out,mva_func,varargin)
% Main analysis engine for SPM-based analysis on
% a local multivariate analysis for PCM function
%
% Adapted from lmva_spm by Rhys Yewbrey, 03/2023
% 
% The analysis function needs to take as the first input argument a PxN
% (voxel x trials) data matrix, and then can take an number of additional
% parameters as input 
% INPUT:
%   S:          Search light definition. This can either be a
%       - file name of a mat file, containing the variables
%               LI,voxmin,voxmax,and (optional) vox (see below)
%       - cell array that contains:
%               function handle to voxelselection function, followed by input parameters
%       - structure that contains (vox is optional for volume base write-out)
%               LI:         Px1 cell array with linear indices
%               voxmin:     Px3 Minimal voxelcoordinate in x,y,z direction
%               voxmax:     Px3 maximal voxelcoordinate in x,y,z direction
%               vox:        Px3 matrix of I,J,K voxle coordinates of centers of search-lights (empty
%                           for surface based search lights)
%
%   P:          Input images
%   out:        output images
%   mva_func:   multivariate analysis function, should take X (data, P x N)
%               as the first argument 
% VARARGIN:
%   params:     extra parameters passed to the mva-function
%   subset:     do the analysis only on subset of nodes 
%   normalise:  voxel-wise normalization of the beta-weights with image (usually ResMS.img) 
% JD June 2010
% v1.1: normalization option implemented 


params={};  % Extra parameter pased to mva_func
Vin = cell(1,length(P)); %preallocate
if (iscell(P)) %P must be a 1 x SN cell containing beta filenames
    for i = 1:length(P)%for participant
        Vin{1,i} = spm_vol(char(P{i})); %to be read by spm_vol
    end%for participant
else
    error('Requires filenames for each participant as cells.')
end%if cell

subset=[]; 
normalise=[]; 
isNP=1;                 % Assumes that the MVA function expects an P x N matrix 

vararginoptions(varargin,{'params','subset','normalise','isNP'});

% Check and resolve the searchlight input:
if (ischar(S))
    S=load(S);
elseif (isstruct(S))
elseif (iscell(S))
    % resolve function call.
end

if(~isempty(normalise)) 
    if (ischar(normalise)) 
        Vnorm=spm_vol(normalise); 
    else 
        Vnorm=normalise; 
    end
end

if (~iscell(S.LI) || size(S.LI,2)~=1)
    error('IJK must be a Px1 cell array');
end

if (~iscell(out))
    error('outfiles must be given in a cell array'); 
end 
num_out=length(out);
num_cent=size(S.LI,1);

% Check the output type 
if (isfield(S,'voxel'))
    type=1;                           
    voxel=S.voxel; 
elseif (isfield(S,'vox'))
    if (size(S.vox,2)~=3 || size(S.vox,1)~=num_cent)
        error('Vox needs to be a Px3 matrix');
    end
    type=1;
    voxel=surfing_subs2inds(Vin{1}(1).dim,S.vox(:,1:3));
else 
    type=2;           % Surface-based searchS.LIght: Write as metric files
end


% Now split the search volume in multiple Blocks and save resulting
% structures as temp files to preserve memory
isdone=false;
if (isempty(subset))
    subset=true(num_cent,1); 
end
iscalc=isnan(S.voxmin(:,1)) | ~subset;                        % Set all NaN to calculated
idealblock=7e5;                                              % Ideal block size
blocksize=max(S.voxmax-S.voxmin);       % Minimal block size 
ratio=idealblock./(prod(blocksize)*length(P)); 
if (ratio>1)
    blocksize=blocksize*ratio^(1/3); 
end 
b=1;

while (~isdone)
    
    % Find the candidate search lights
    % This starts in the x-direction and finds the slice of search lights
    % that are within that slice, and similar for y and z.
    % then it moves the block tightly to this location
    iscand=~iscalc;
    
    for i=1:3
        IJKc1(i)=min(S.voxmin(iscand,i));
        IJKc2(i)=IJKc1(i)+blocksize(i);     
        iscand=iscand & S.voxmin(:,i)>=IJKc1(i) & S.voxmax(:,i)<=IJKc2(i);
    end
    IJKc1=min(S.voxmin(iscand,:),[],1);
    IJKc2=IJKc1+blocksize;
    iscand=~iscalc & S.voxmin(:,1)>=IJKc1(1)  & S.voxmin(:,2)>=IJKc1(2) & S.voxmin(:,3)>=IJKc1(3) & ...
        S.voxmax(:,1)<=IJKc2(1)  & S.voxmax(:,2)<=IJKc2(2) & S.voxmax(:,3)<=IJKc2(3);
    j=find(iscand);
    if (isempty(j))
        break;
    end
    
    
    % Save the substructure as a tempory file
    T.LI={S.LI{j}}';
    T.j=j;
    if (isfield(S,'vox'))
        T.vox=S.vox(j,:);
    end
    fprintf('block %d Corner: %d %d %d length:%d  \n',b,IJKc1,length(j));
    save(sprintf('temp_%2.2d.mat',b),'T');
    iscalc(j)=true;
    isdone=all(iscalc);
    b=b+1;
    
end
numblocks=b-1;

% Free the memory
clear S T;

% InitiaS.LIze output arrays
R=zeros(num_out,num_cent)*NaN;

k=1;
spm_progress_bar('Init',numblocks,'MVA','number of Blocks');
tic;
for b=1:numblocks
    load(sprintf('temp_%2.2d.mat',b),'T');
    
    linVox=unique(cat(2,T.LI{:})');
    
    for subj=1:length(Vin)
        [I,J,K]=ind2sub(Vin{subj}(1).dim,linVox);
        
        X{subj}=sparse(double(max(linVox)),length(Vin{subj}));
        N{subj}=sparse(double(max(linVox)),1);
        
        for i=1:length(Vin{subj})
            X{subj}(linVox,i)=spm_sample_vol(Vin{subj}(i),double(I),double(J),double(K),0);
        end
    end
    if (~isempty(normalise)) 
        N(linVox,1)=spm_sample_vol(Vnorm,double(I),double(J),double(K),0);
        X(linVox,:)=bsxfun(@rdivide,X(linVox,:),sqrt(N(linVox,:))); 
    end
    clear I J K;
    for i=1:size(T.LI,1)
        if (isNP)
            for subj=1:length(X)
                
                Xdata{subj} = full(X{subj}(T.LI{i},:)');
            end
            R(:,T.j(i))=feval(mva_func,Xdata,params{:});
        else
            for subj=1:length(X)
                Xdata{subj} = full(X{subj}(T.LI{i},:));
            end
            R(:,T.j(i))=feval(mva_func,Xdata,params{:});
        end
        if (mod(k,50)==0)
            fprintf('%d done: %f\n',k,toc);
            tic;
        end
        k=k+1;
    end
    spm_progress_bar('Set',b);
    b=b+1;
end
fprintf('%d total: %f\n',k,toc);

% Clean up the search lights
delete('temp*.mat');



if (isempty(out))
    return;
end

if (type==1)
    % Write as a volume based image
    
    for i=1:size(R,1)
        Z=NaN(Vin(1).dim);
        [d, f, t m]=spm_fileparts(out{i});
        Vo      = struct(...
            'fname',    fullfile(d,[f t]),...
            'dim',      Vin(1).dim,...
            'dt',       [spm_type('float32') spm_platform('bigend')],...
            'mat',      Vin(1).mat,...
            'n',        [str2num(m) 1],...
            'descrip',  'multivariate analysis');
        Z(voxel)=R(i,:);            % Project back into voxel space 
        spm_write_vol(Vo,Z);
    end
else 
    for i=1:size(R)
        [d, f, t, m]=spm_fileparts(out{i});
        metricNames{i}= fullfile(d, [f t]);
        metricColumn{i}= m(2:end);
        if (~strcmp(t,'.metric'))
            error('file type must me a metric file'); 
        end
    end
    %make the struct
    M=struct();
    for i=1:size(R)     % Loop over all rows in the MVA-result
        j=1;
        endwhile=  numel(M);
        while j<= endwhile                  % Loop over all existing metric files 
            if ~isfield(M, 'save_names')    % If there is no one yet, initialize the first 
                M.save_names= metricNames{1};
                M.column_name= {metricColumn{1}}; 
                M.data=R(1,:)';
                M.num_rows=num_cent;
                M.num_cols= 1;
                minmax_colormapping=[-1 1];
                M.column_color_mapping=minmax_colormapping;
                M.encoding={'BINARY'};
                M.index=(0:(num_cent-1))';
                j=j+1;
            elseif strcmp(M(j).save_names, metricNames{i})      % If a metric file with this name is already existing: add column
                M(j).column_name= [M(j).column_name, {metricColumn{i}}];
                M(j).data=[M(j).data R(i,:)'];
                M(j).num_cols= M(j).num_cols+1;
                M(j).column_color_mapping=repmat(minmax_colormapping,M(j).num_cols,1);
                j=numel(M)+1;
            elseif j==numel(M)          % If not, generate additional metric file 
                M(j+1).save_names= metricNames{i};
                M(j+1).column_name= {metricColumn{i}};
                M(j+1).data=R(i,:)';
                M(j+1).num_rows=num_cent;
                M(j+1).num_cols= 1;
                minmax_colormapping=[-1 1];
                M(j+1).column_color_mapping=minmax_colormapping;
                M(j+1).encoding={'BINARY'};
                M(j+1).index=(0:(num_cent-1))';
                j=j+1;
            else
                j=j+1;
            end
        end
    end
    for i=1:numel(M)
        caret_savemetric(M(i).save_names, rmfield(M(i),'save_names') );   %flexible
    end
    
end
spm_progress_bar('Clear');

