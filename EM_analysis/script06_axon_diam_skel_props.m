function script06_axon_diam_skel_props(samplei) 
% compute skeletons and cross seccional areas on aligned axons
% and also check several properties of the axons
% samplei = the index of the EM sample to analyze
%
% This script needs the original volumes and the results from scripts:
% script01_props3d.m
% script05_align_axons.m
% It also needs the functions in: 
% RefineSeg.m
% The results of this script will be saved in folder
% White_matter_EM_Diam_Skel_Props/*/
% 29-Aug-2022 by Ricardo Coronado-Leija
if(nargin == 0)
close all
clc
rng(0)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Begin\n')
root       = ''; 
ir         = 2; % 1 = HM, 2 = LM
res        = {'HM','LM'};
rootIn     = [root 'White_matter_EM_Aligned_Axons_' res{ir} '/']; 
rootOut    = [root 'White_matter_EM_Diam_Skel_Props_' res{ir} '/']; 
rootInfo   = [root 'White_matter_EM_Prop/'];
samples    = {'Sham_25_contra','Sham_25_ipsi','Sham_49_contra','Sham_49_ipsi',...
              'TBI_24_contra','TBI_24_ipsi','TBI_28_contra','TBI_28_ipsi',...
              'TBI_2_contra','TBI_2_ipsi'};
idsample   = [25 25 49 49 24 24 28 28 2 2 ]; 
sidesample = {'contra','ipsi','contra','ipsi','contra','ipsi',...
              'contra','ipsi','contra','ipsi'};
nsamples   = length(samples);
voxsizes   = [15 50]*1e-3;
pminv      = 0.10; % min proportion volume (broken axons)
pmina      = 0.05; % min proportion area (broken cross sections)
rmin       = 2;    % pixels (for closing operation)
lminum     = 0.5;  % um 
lmin       = round( lminum/voxsizes(ir) ); % lmin pixels
Lmin       = 2*lmin+1; % only axons larger than this length (1um) will survive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select the samples to loop
if(nargin > 0)
fidx = samplei;     
else
fidx = 1:nsamples;
end 
% % plot flag
% fPlot = rand < 0;
% % matrix to rotate  z -> x
% RotXtoZ = [0 0 -1; 0 1 0; 1 0 0];
%
for f = fidx
sample     = [num2str(idsample(f)) '_' sidesample{f}];    
% ====================== % 
% ===== Read Vols ====== %
% ====================== %
% create output folder if needed      
folderIn   = [rootIn   '/' samples{f} '/'];
folderInfo = [rootInfo '/' samples{f} '/'];
folderOut  = [rootOut  '/' samples{f} '/'];  
if(~isfolder(folderOut)); mkdir(folderOut); end
tic
fprintf('Read statistics ... ')
sname        = [folderInfo '/' res{ir} '_' sample '_myelinated_axons_stats.mat'];
labels_stats = load(sname).labels_stats;
%
fprintf('Read aligned axons indices ... ')
fbname = [folderIn '/' res{ir} '_' sample '_my_alig_axons_ind_list.mat'];
fb     = load(fbname).aligned_fibers_indices;
nobj   = size(fb,1);
voxs   = voxsizes(ir); % voxel size micras isotropic
toc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
tic
fprintf('Compute properties %s ... ',samples{f}) 
% for i = [1 10 50 100 150 200 216 250 279 293 300 350] %
for i = nobj:-1:1 % 216
% ======================== % 
% ======= Set Axon ======= %
% ======================== %    
% Size of axon bounding box    
AxSize  = fb(i,4).SizeNew;
% axons longer than 1um in z
% empty axons have size (1,1,1) so this also excludes them
if(AxSize(3) >  Lmin )
% Indices axon
AxIdx   = fb(i,1).VoxelIdxList{:};
% Volume
fiber   = zeros(AxSize,'logical');
% Set fiber
fiber(AxIdx) = true;
% RefineSeg.PlotAxon(fiber,i,sample,voxs(1));
% ======================== % 
% = Volumes and Box Size = %
% ======================== %   
% Check 00: Change volume beetween original and aligned
volori  = labels_stats(i,1).Volume; 
volnew  = length(AxIdx);            
volnewp = volnew/volori;    
clear('AxIdx')
% Check 01: Check number of objects composing the aligned fiber
[fiberu,nobj,vox1st,voxrest,rest1stp] = RefineSeg.KeepLargestComponentsMore(fiber,pminv);
clear('fiber');
% proportion with original volume
vol1stp  = vox1st/volori;
volkept  = (vox1st+voxrest)/volori; 
% Check 02: checking length box, max(xy)/z ratio
xrange0  = find(sum(sum(fiberu,2),3));
yrange0  = find(sum(sum(fiberu,1),3));
zrange0  = find(sum(sum(fiberu,1),2));
xlength  = length(xrange0);
ylength  = length(yrange0);
zlength  = length(zrange0);
xyzratio = max([xlength ylength])/zlength;
clear('xrange0','yrange0');
% ======================== % 
% == Skeleton CS Areas === %
% ======================== %  
Axon      = RefineSeg.AnalyzeCrossSections(fiberu,zrange0,lmin,rmin,pmina);
volfinal  = length(Axon.AxonVoxelIdxList); 
volfinalp = volfinal/volori;
clear('fiberu');
% ======================== % 
% ====== Save Axon ======= %
% ======================== %  
Axon.propvolorifinal    = volfinalp; % proportion final axon / original axon
% Check 00
Axon.volumeoriginal     = volori;   % Original Axon Volume
Axon.volumealigned      = volnew;   % Original Axon Volume
Axon.propvolorialigned  = volnewp;  % Original Axon Volume
% Check 01
Axon.numobjectsaligned  = nobj;      % number of objects in the original axon
Axon.volumelargest1st   = vox1st;    % volume largest object
Axon.volumelargestrest  = voxrest;   % volume sum other largest objects (kept but smaller than largest)
Axon.propvolorikept     = volkept;   % proportion (total kept objects volume)/original axon volume
Axon.propvolori1st      = vol1stp;   % proportion largest/original axon volume
Axon.propvol1strest     = rest1stp;  % proportion largest/ other object volumes
% Check 02
Axon.largestrangez      = uint32(zrange0);  % range of axon(largest object) in z-axis 
Axon.largestboundingbox = [xlength ylength zlength]; % bounding box size axon (largest object)
Axon.largestxyzratio    = xyzratio; % ratio between max(xy)/z lengths of bounding box axon (largest object)
% Legacy Variables
Axon.SizeOriginal   = fb(i,2).SizeOriginal;
Axon.SizeOld        = fb(i,3).SizeOld;
Axon.SizeNew        = fb(i,4).SizeNew;
Axon.OriginOld      = fb(i,5).OriginOld;
Axon.OriginNew      = fb(i,6).OriginNew;
Axon.MatAlign       = fb(i,7).MatAlign{:};
Axon.MatRestore     = fb(i,8).MatRestore{:};
Axon.VoxelSize      = voxs;%fb(i,9).VoxelSize;
% set into bigger structure   
Axons(i) = Axon;
% ----------------------------------------------------------------------- %
fprintf('%u ',i);
labels_stats(i,:) = [];
fb(i,:)           = [];
clear('Axon','zrange0');
end % no empty    
end % i
toc
fprintf('\n');
% ======================== % 
% ====== Structures ====== %
% ======================== %
tic
% transform structure to table  
Props                  = rmfield(Axons,{'AxonVoxelIdxList','skeleton','csareas','csareasalign'});
axon_skel_areacs_props = struct2table(Axons);  
clear('Axons');
skeletons              = axon_skel_areacs_props(:,3);     
areas_cross_sections   = axon_skel_areacs_props(:,[8 9]); 
axon_props             = struct2table(Props); 
clear('Props')
% save
toc
% ======================== % 
% ====== Save Axons ====== %
% ======================== %
tic
fprintf('Save structure ... ')
fbname = [folderOut '/' res{ir} '_' sample '_axons_skel_areacs_props.mat'];
save(fbname,'axon_skel_areacs_props','-v7.3');
clear('axon_skel_areacs_props')
%
fbname = [folderOut '/' res{ir} '_' sample '_axons_skeletons.mat'];
save(fbname,'skeletons');
clear('skeletons');
%
fbname = [folderOut '/' res{ir} '_' sample '_axons_areas_cross_sections.mat'];
save(fbname,'areas_cross_sections');
clear('areas_cross_sections');
%
fbname = [folderOut '/' res{ir} '_' sample '_axons_props.mat'];
save(fbname,'axon_props');
clear('axon_props');
%
clear('labels_idx_list','labels_stats');
toc
end % f
fprintf('End\n');
end % main