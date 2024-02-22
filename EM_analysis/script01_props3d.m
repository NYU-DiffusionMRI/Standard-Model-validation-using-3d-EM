function script01_props3d(samplei) 
% use regionprops3 on the segmented EM samples
% samplei = the index of the EM sample to analyze
%
% This script needs the original DeepAcson segmentations in folder
% White_Matter_EM/*/
% These can be downloaded from:
% https://etsin.fairdata.fi/dataset/f8ccc23a-1f1a-4c98-86b7-b63652a809c3
% The results of this script will be saved in folder
% White_matter_EM_Prop/*/
%
% 01-Jun-2021 by Ricardo Coronado-Leija
if(nargin == 0)
close all
clc
rng(0)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Begin\n')
rootIn     = 'White_matter_EM/';
rootOut    = 'White_matter_EM_Prop/'; 
samples    = {'Sham_25_contra','Sham_25_ipsi','Sham_49_contra','Sham_49_ipsi',...
              'TBI_24_contra','TBI_24_ipsi','TBI_28_contra','TBI_28_ipsi',...
              'TBI_2_contra','TBI_2_ipsi'};
idsample   = [25 25 49 49 24 24 28 28 2 2 ]; 
sidesample = {'contra','ipsi','contra','ipsi','contra','ipsi',...
              'contra','ipsi','contra','ipsi'};
nsamples   = length(samples);
res        = {'HM','LM'};
ir         = 2; % 1 = HM, 2 = LM
voxsizes   = [15 15 50; 50 50 50]*1e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select the samples to loop
if(nargin > 0)
fidx = samplei;     
else
fidx = 1:nsamples;
end 
%
for f = fidx   
sample    = [num2str(idsample(f)) '_' sidesample{f}];
% create output folder if needed 
folderOut = [rootOut '/' samples{f} '/'];  
if(~isfolder(folderOut)); mkdir(folderOut); end
tic
fprintf('Reading data ... ')
% % raw
% iname    = [rootIn '/' samples{f} '/' res{ir} '_' sample '.h5'];
% idat     = squeeze(h5read(iname,'/raw'));
% myelin
mname    = [rootIn '/' samples{f} '/' res{ir} '_' sample '_myelin.mat'];
mdat     = load(mname);
% axons
aname    = [rootIn '/' samples{f} '/' res{ir} '_' sample '_myelinated_axons.mat'];
adat     = load(aname);
% nucleus
if(ir == 2)
nname    = [rootIn '/' samples{f} '/' res{ir} '_' sample '_nucleus.mat'];
ndat     = load(nname);    
end 
toc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
voxs      = voxsizes(ir,:);  % voxel size micras
if(ir == 1)
% resample High Resolution Data to be of isotropic resolution 
vref     = min(voxs);
nxi      = voxs(1)/vref;
nyi      = voxs(2)/vref;
nzi      = voxs(3)/vref;
% data     = imresize3(idat,'Scale',[nxi nyi nzi],'Method','linear');  
% clear('idat');
labels   = imresize3(adat.myelinated_axons,'Scale',[nxi nyi nzi],'Method','nearest') ;
clear('adat');
myelins  = imresize3(mdat.myelin,'Scale',[nxi nyi nzi],'Method','nearest') > 0.5;   
clear('mdat');
vox      = [vref vref vref];
isovox   = mean(vox);    
% ncells   = 0;
else
% set global name labels and combine with nucleus if needed    
tic
fprintf('Combining nucleus ... ')    
naxons       = uint32( max(adat.final_lbl(:)) );  % axons
ncells       = uint32( max(ndat.final_lbl(:)) );  % cells
% data         = idat;                              % data
% clear('idat');
labels       = uint32( adat.final_lbl );          % labels axons
clear('adat');
nucleus      = uint32( ndat.final_lbl ) + naxons; % new labels cell
clear('ndat');
mask         = nucleus > naxons;                  % mask of cells
labels(mask) = nucleus(mask);                     % adding cells to axon labels 
myelins      = mdat.myelin;
clear('mdat');
vox          = voxsizes(ir,:);
isovox       = mean(vox); 
clear('nucleus','mask');   
fprintf('axons = %d, nucleus = %d ..',naxons,ncells);
toc
end % if-else 
[nx,ny,nz] = size(myelins);
clear('myelins')
% compute number of objects and their volumes
nobj = max(labels(:));
fprintf('%s_%s: nx:%d, ny:%d, nz:%d, res:%.3f um, axons = %d\n',...
        samples{f},res{ir},nx,ny,nz,isovox,nobj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volume   Centroid   BoundingBox   PrincipalAxisLength   Orientation
% EigenVectors   EigenValues   SurfaceArea
%
% EigenValues  = (Lengths/4)^2;
% EigenVectors = eul2rotm(Orient*pi/180)
% Orientation  = rotm2euler(Vecs)*180/pi;

% Compute and save statistics
% Computation of solidity takes very very long time
fprintf('Compute and save statistics ... ')
tic
labels_stats = regionprops3(labels,'BoundingBox','Centroid',...
               'EigenValues','EigenVectors','Orientation','PrincipalAxisLength', ...
               'SurfaceArea','Volume');%,'Solidity');
sname        = [folderOut '/' res{ir} '_' sample '_myelinated_axons_stats.mat'];
save(sname,'labels_stats','-v7.3');
toc
clear('labels_stats')

% Does not take much time and it takes a lot of storage space
% But this one will be saved anyways
% This will help for further processing
fprintf('Get and save voxel idx list ... ')
tic
labels_idx_list = regionprops3(labels,'VoxelIdxList');
iname           = [folderOut '/' res{ir} '_' sample '_myelinated_axons_idx_list.mat'];
save(iname,'labels_idx_list','-v7.3');               
toc
clear('labels_idx_list')

clear('labels');
end % f

fprintf('End\n')
end % main