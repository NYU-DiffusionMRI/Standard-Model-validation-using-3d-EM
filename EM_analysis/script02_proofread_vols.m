function script02_proofread_vols(samplei) 
% merge axons, myelin & cells into one volume and remove small spurious 
% objects from the segmented volumes for better visualization
% samplei = the index of the sample to analyze
%
% This script needs the original DeepAcson segmentations in folder
% White_Matter_EM/*/
% These can be downloaded from:
% https://etsin.fairdata.fi/dataset/f8ccc23a-1f1a-4c98-86b7-b63652a809c3
% It also needs the functions in: 
% RefineSeg.m
% The results of this script will be saved in folder
% White_matter_EM_FullSubstrates/*/
%
% 02-Feb-2022 by Ricardo Coronado Leija
if(nargin == 0)
close all
clc
rng(0)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Begin\n')
rootIn     = 'White_matter_EM/';
rootOut    = 'White_matter_EM_FullSubstrates/'; 
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
if(ir == 1)
% raw
iname    = [rootIn '/' samples{f} '/' res{ir} '_' sample '.h5'];
idat     = squeeze(h5read(iname,'/raw'));
end % only for HM
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
% Binary volumes
if(ir == 1)
% resample High Resolution Data to be of isotropic resolution 
vref     = min(voxs);
nxi      = voxs(1)/vref;
nyi      = voxs(2)/vref;
nzi      = voxs(3)/vref;
data     = imresize3(idat                     ,'Scale',[nxi nyi nzi],'Method','linear')  ;
clear('idat');
axons    = imresize3(adat.myelinated_axons > 0,'Scale',[nxi nyi nzi],'Method','nearest') ;
clear('adat');
myelins  = imresize3(mdat.myelin           > 0,'Scale',[nxi nyi nzi],'Method','nearest') ;
clear('mdat');
nucleus  = [];
vox      = [vref vref vref];
isovox   = mean(vox);    
else
% set global names
tic
data     = [];
axons    = adat.final_lbl > 0;
clear('adat')
myelins  = mdat.myelin    > 0;
clear('mdat');
nucleus  = ndat.final_lbl > 0;
clear('ndat');
vox      = voxsizes(ir,:);
isovox   = mean(vox); 
toc
end % if-else 
[nx,ny,nz] = size(myelins);
% compute number of objects and their volumes
fprintf('%s_%s: nx:%d, ny:%d, nz:%d, res:%.3f um\n',samples{f},res{ir},nx,ny,nz,isovox)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% voxel size ( in um )
vs      = isovox; 
% to dilate ias (in um) 
rs      = 0.350;  
% to remove small objects
nvoxels = 100;   % first removal of isolated tiny objects before dilate (in pixels)
dm      = 0.300; % diameter small spherical objects (in um) 

% save resized data for HM
if(ir == 1)
dname = [folderOut '/' res{ir} '_' sample '_data.h5'];
h5create(dname,'/raw',[nx ny nz]);
h5write(dname ,'/raw',data)
clear('data')
end % ir

% set original substrate
tic 
fprintf('Original Substrate ... ')
SubstrateOriginal = RefineSeg.SetSubstrate(axons,myelins,nucleus);
toc
% save
tic
soname = [folderOut '/' res{ir} '_' sample '_substrate_original.mat'];
save(soname,'SubstrateOriginal','-v7.3'); 
toc
clear('SubstrateOriginal')

% set cleaned substrate
tic 
fprintf('Clean Substrate ... ')
SubstrateCleaned = RefineSeg.SetCleanSubstrate(axons,myelins,nucleus,vs,rs,dm,nvoxels);
toc
clear('axons','myelins','nucleus')
% save
tic
scname = [folderOut '/' res{ir} '_' sample '_substrate_cleaned.mat'];
save(scname,'SubstrateCleaned','-v7.3'); 
toc
clear('SubstrateCleaned')

end % f

fprintf('End\n')
end % main