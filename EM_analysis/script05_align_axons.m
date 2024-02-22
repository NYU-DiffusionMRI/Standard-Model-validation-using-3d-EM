function script05_align_axons(samplei) 
% This script aligns the axons to the z axis, it saves the transformations 
% Original->Aligned and Aligned->Original. Instead of rotate coordinates, 
% the final rotated volume is "interpolated" using the unrotated volume 
% as reference so the original can be recovered in a similar way.
% samplei = the index of the EM sample to analyze
%
% This script needs the original volumes and the results from script:
% script01_props3d.m
% It also needs the functions in: 
% RefineSeg.m
% and the information in the text file
% 000_animal_exp_axons_cells_nx_ny_nz_cclims_cglims.txt
% The results of this script will be saved in folder
% White_matter_EM_Aligned_Axons/*/
% 28-Aug-2022 by Ricardo Coronado-Leija
if(nargin == 0)
close all
clc
rng(0)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Begin\n')
root       = '';  
rootInfo   = '';
res        = {'HM','LM'};
ir         = 2; % 1 = HM, 2 = LM
rootOri    = [root 'White_matter_EM/'];
rootIn     = [root 'White_matter_EM_Prop/']; 
rootOut    = [root 'White_matter_EM_Aligned_Axons_' res{ir} '/']; 
Info       = '000_animal_exp_axons_cells_nx_ny_nz_cclims_cglims.txt';
samples    = {'Sham_25_contra','Sham_25_ipsi','Sham_49_contra','Sham_49_ipsi',...
              'TBI_24_contra','TBI_24_ipsi','TBI_28_contra','TBI_28_ipsi',...
              'TBI_2_contra','TBI_2_ipsi'};
idsample   = [25 25 49 49 24 24 28 28 2 2 ]; 
sidesample = {'contra','ipsi','contra','ipsi','contra','ipsi',...
              'contra','ipsi','contra','ipsi'};
nsamples   = length(samples);
voxsizes   = [15 15 50; 50 50 50]*1e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select the samples to loop
if(nargin > 0)
fidx = samplei;     
else
fidx = 1:nsamples;
end 
% plot flag
fPlot = rand < 0;
% matrix to rotate  z -> x
RotXtoZ = [0 0 -1; 0 1 0; 1 0 0];
%
for f = fidx
% ====================== % 
% ===== Read Vols ====== %
% ====================== %
% create output folder if needed    
sample    = [num2str(idsample(f)) '_' sidesample{f}];
folderOri = [rootOri '/' samples{f} '/'];
folderIn  = [rootIn  '/' samples{f} '/'];
folderOut = [rootOut '/' samples{f} '/'];  
if(~isfolder(folderOut)); mkdir(folderOut); end
tic
% Read myelin 
fprintf('Read myelin ... ')
mname           = [folderOri '/' res{ir} '_' sample '_myelin.mat'];
mdat            = load(mname).myelin;
% Read statistics 
fprintf('Read statistics ... ')
sname           = [folderIn  '/' res{ir} '_' sample '_myelinated_axons_stats.mat'];
labels_stats    = load(sname).labels_stats;
% Read indices
fprintf('Read voxel idx list ... ')
iname           = [folderIn  '/' res{ir} '_' sample '_myelinated_axons_idx_list.mat'];
labels_idx_list = load(iname).labels_idx_list;
toc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% ====================== % 
% ===== Resize Vol ===== %
% ====================== %
tic
voxs    = voxsizes(ir,:); % voxel size micras
if(ir == 1)
% resample High Resolution Data to be of isotropic resolution 
vref    = min(voxs);
nxi     = voxs(1)/vref;
nyi     = voxs(2)/vref;
nzi     = voxs(3)/vref;
myelins = imresize3(mdat,'Scale',[nxi nyi nzi],'Method','nearest') > 0.5;   
clear('mdat');
vox     = [vref vref vref];
isovox  = mean(vox);    
ncells  = 0;
else
myelins = mdat;
clear('mdat');
vox     = voxs;
isovox  = mean(vox); 
cglims  = load([rootInfo Info]); 
ncells  = cglims(4,f); % num cells
end % depending on volumes    
[nx,ny,nz] = size(myelins); % size volume
sz         = [nx ny nz];
nobjects   = size(labels_idx_list,1) - ncells; % number of objects (no cells)
toc    
clear('myelins')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Volume   Centroid   BoundingBox   PrincipalAxisLength   Orientation
% EigenVectors   EigenValues   SurfaceArea
%
% EigenValues  = (Lengths/4)^2;
% EigenVectors = eul2rotm(Orient*pi/180)
% Orientation  = rotm2euler(Vecs)*180/pi;
%
tic
fprintf('Align fibers %s ... ',samples{f}) 
for i = nobjects:-1:1 % so allocation ocurrs first
% tic
% indices of the object    
indj        = labels_idx_list(i,1).VoxelIdxList{:};
% check emptiness
if(~isempty(indj))
% ====================== % 
% ===== Props Axon ===== %
% ====================== %    
% length of major axes of the ellipsoid that has the same normalized second
% central moments as the object (not necesarily equal to length of object)
% Length      = labels_stats(i,4).PrincipalAxisLength*isovox;
% eigenvectors (rotation matrix from x to orientation of the object)
EigVec      = labels_stats(i,6).EigenVectors{:};
% Bounding Box (integers)
BoundingBox = round( labels_stats(i,3).BoundingBox );
% Rotation Matrices
% (first to x, inverting eigenvectors and then from x->z)
MatAlign    = RotXtoZ*(EigVec');
MatRestore  = EigVec*(RotXtoZ');
% ====================== % 
% ===== Align Axon ===== %
% ====================== %
% get coordinates of the axon
[x,y,z]    = ind2sub(sz,indj);
subj       = [x y z]; 
% [subjnewI,indjnew,subjori,sznew] = RefineSeg.GetAlignedAxon(subj,MatAlign);
[~,indjnew,newori,sznew] = RefineSeg.GetAlignedAxon(subj,MatAlign);
% For later only need: subjori and sznew (indjnew is just for testing)
% ======================== % 
% === Original Volume ==== %
% ======================== %
szold  = BoundingBox([5 4 6]);
oldori = BoundingBox([2 1 3]) - 1;
fiber1 = RefineSeg.SetAxonOrigin(indj,sz,szold,oldori);
% ======================== % 
% === Rotated Volume ===== %
% ======================== %
finalindices = RefineSeg.GetAlignedAxonBetter(fiber1,sznew,oldori,newori,MatRestore);
% get struct with the axon indices (large)
fb(i).VoxelIdxList   = finalindices;
fb(i).SizeOriginal   = sz;
fb(i).SizeOld        = szold;
fb(i).SizeNew        = sznew;
fb(i).OriginOld      = oldori;
fb(i).OriginNew      = newori;
fb(i).MatAlign       = MatAlign;
fb(i).MatRestore     = MatRestore;
fb(i).VoxelSize      = isovox;
% get struct without the axon indices (small)
tr(i).SizeOriginal   = sz;
tr(i).SizeOld        = szold;
tr(i).SizeNew        = sznew;
tr(i).OriginOld      = oldori;
tr(i).OriginNew      = newori;    
tr(i).MatAlign       = MatAlign;
tr(i).MatRestore     = MatRestore;
tr(i).VoxelSize      = isovox;
% ----------------------------------------------------------------------- %
if(fPlot)
% ======================== % 
% ===== Set New Axon ===== %
% ======================== %
fiber2               = zeros(sznew,'logical');
fiber2(finalindices) = 1;
fiber3               = zeros(sznew,'logical');
fiber3(indjnew)      = 1;
% ======================== % 
% ===== Restore Axon ===== %
% ======================== %
oldindices         = RefineSeg.GetAlignedAxonBetter(fiber2,szold,newori,oldori,MatAlign);
fiber4             = zeros(szold,'logical');
fiber4(oldindices) = 1;
% ======================== % 
% ==== Compare Axons ===== %
% ======================== %
figure; sliceViewer(fiber1 + fiber4); % Original vs Restored
figure; sliceViewer(fiber2 + fiber3); % Simple Rotation vs Better Rotation 
fprintf('%u ',i);
end % length > 15
% ----------------------------------------------------------------------- %
fprintf('%u ',i);
else
% ======================== % 
% ====== Empty Axons ===== %
% ======================== %
% get struct with the axon indices (large)
fb(i).VoxelIdxList   = 1;
fb(i).SizeOriginal   = sz;
fb(i).SizeOld        = [1 1 1];
fb(i).SizeNew        = [1 1 1]; 
fb(i).OriginOld      = [1 1 1];
fb(i).OriginNew      = [1 1 1];
fb(i).MatAlign       = MatAlign;
fb(i).MatRestore     = MatRestore;
fb(i).VoxelSize      = isovox;
% get struct without the axon indices (small)
tr(i).SizeOriginal   = sz;
tr(i).SizeOld        = [1 1 1]; 
tr(i).SizeNew        = [1 1 1]; 
tr(i).OriginOld      = [1 1 1];
tr(i).OriginNew      = [1 1 1];   
tr(i).MatAlign       = MatAlign;
tr(i).MatRestore     = MatRestore;
tr(i).VoxelSize      = isovox;
end % no empty    
clear('fiber1','subj','indj','indjnew','x','y','z','finalindices')
labels_idx_list(i,:) = [];
labels_stats(i,:)    = [];
% pause(0.01);
end % i
fprintf('\n');
clear('labels_idx_list','labels_stats');
% ======================== % 
% ====== Structures ====== %
% ======================== %
% transform structure to table       
aligned_fibers_indices     = struct2table(fb);          
aligned_fibers_parameters  = struct2table(tr);       
clear('fb','tr')
% save
toc
% ======================== % 
% ====== Save Axons ====== %
% ======================== %
tic
fprintf('Save structures ... ')
fbname = [folderOut '/' res{ir} '_' sample '_my_alig_axons_ind_list.mat'];
save(fbname,'aligned_fibers_indices','-v7.3');
clear('aligned_fibers_indices')
%
fbname = [folderOut '/' res{ir} '_' sample '_my_alig_axons_pars.mat'];
save(fbname,'aligned_fibers_parameters');
clear('aligned_fibers_parameters');
toc
end % f
fprintf('End\n');
end % main

