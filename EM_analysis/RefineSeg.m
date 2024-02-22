classdef RefineSeg
% set of function for the analysis of segmented 3d axons from EM
% 31-Jan-2022 by Ricardo Coronado-Leija
methods (Static)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotAxon(fiber,flabel,samplesN,voxsize)
[nx,ny,nz] = size(fiber);
tstr = [samplesN ' Axon: ' num2str( flabel ) ' Length Z: ' num2str(nz*voxsize,'%.2f') ... 
    ' Length X: ' num2str(nx*voxsize,'%.2f') ' Length Y: ' num2str(ny*voxsize,'%.2f')];
fprintf('%s\n',tstr)
% figure('Name','axon','Visible','off')
figure
subplot(1,3,1)
h = RefineSeg.plotaxon(fiber); 
xlabel('x')
ylabel('y')
zlabel('z')
view([90 0])
set(h,'edgealpha',0,'facealpha',0.6,'facecolor',[0.3 0.3 0.3])
axis equal off
hold on
material dull; lightangle(-45,30)
subplot(1,3,2)
h = RefineSeg.plotaxon(fiber); 
xlabel('x')
ylabel('y')
zlabel('z')
view([0 0])
set(h,'edgealpha',0,'facealpha',0.6,'facecolor',[0.3 0.3 0.3])
axis equal off
hold on
material dull; lightangle(-45,30)
subplot(1,3,3)
h = RefineSeg.plotaxon(fiber); 
xlabel('x')
ylabel('y')
zlabel('z')
view([0 90])
set(h,'edgealpha',0,'facealpha',0.6,'facecolor',[0.3 0.3 0.3])
axis equal off
hold on
material dull; lightangle(-45,30)
sgtitle(tstr,'fontsize',30,'fontweight','bold')
% close('axon')
end
%
function PlotAxons(fiber,tstr,fslope)
if(nargin<3)
fslope = false;    
end
[nx,ny,nz] = size(fiber);    
fprintf('%s\n',tstr)
% plot axons
RefineSeg.plotaxon(fiber,fslope); 
xlim([1 ny])
ylim([1 nx])
zlim([1 nz])
% view([0 0])
% axis equal 
axis off
hold on
material dull; 
% light('Position',[1 0 0],'style','local'); 
end
%
function plotaxon(BW,fslope)
fv = isosurface(BW,0);
% p2 = patch(fv,'edgecolor','none','edgealpha',0,'facealpha',0.1,'facecolor',[0.2 0.2 0.2]); % MC simulation
p2 = patch(fv,'edgecolor','none','edgealpha',0,'facealpha',0.9,'facecolor',[255 192 0]/255);
isonormals(BW,p2)
view(3) 
daspect([1 1 1])
% axis tight
axis equal 
camlight 
if(fslope)
camlight(-80,-10) % 5 axons
else
view(0,0) % 12 axons
end
lighting gouraud
clear('fv','p2')
end
%
% %
% function PlotAxons(fiber,tstr)
% [nx,ny,nz] = size(fiber);    
% fprintf('%s\n',tstr)
% % figure
% h = RefineSeg.plotaxon(fiber); 
% % xlabel('x')
% % ylabel('y')
% % zlabel('z')
% xlim([1 ny])
% ylim([1 nx])
% zlim([1 nz])
% view([0 0])
% axis equal 
% axis off
% set(h,'edgealpha',0,'facealpha',0.4,'facecolor',[0.3 0.3 0.3])
% hold on
% material dull; light('Position',[1 0 0],'style','local'); %lightangle(-45,30)
% % title(tstr,'fontsize',30,'fontweight','bold')
% end
% %
% function h = plotaxon(BW)
%     [f,v] = isosurface(BW,0);
%     TR    = triangulation(f,v);
%     h     = trisurf(TR);
%     clear('f','v','TR')
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Volume = SetAxonOrigin(AxonIndices,SizeVol,SizeAxon,Origin)
% get original coordinates
[x,y,z]       = ind2sub(SizeVol,AxonIndices);
subj          = [x y z]; 
clear('x','y','z')
% box for the original isolated axon
Volume        = zeros(SizeAxon,'logical');
% set origin at  (1,1,1) 
subj1         = subj - Origin;
clear('subj')
% put original axon in its box
indj1         = sub2ind(SizeAxon,subj1(:,1),subj1(:,2),subj1(:,3)); 
clear('subj1')
Volume(indj1) = true; 
clear('AxonIndices','indj1')
end % SetAxonOrigin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AlignedCoordinates,AlignedIndices,Origin,BoxSize] = ...
        GetAlignedAxon(OriginalCoordinates,RotMat)
% Input:
% OriginalCoordinates = (N x 3)
% RotMat              = (3 x 3)
% Output:
% AlignedCoordinates  = (N x 3)
% AlignedIndices      = (N x 1)
% Origin              = (1 x 3)
% BoxSize             = (1 x 3)

% Align axon to z  
subjrot  = RotMat * OriginalCoordinates' ;
clear('OriginalCoordinates')
% get minimum values
subjori  = min(subjrot,[],2) - 1;
% get integer coordinates starting in 1,1,1
subjnew  = (subjrot - subjori); 
clear('subjrot') 
subjnewI = round(subjnew); 
clear('subjnew') 
% get max values (size)
sznew    = ( max(subjnewI,[],2) )';
% get indices with respect  to new volume
indjnew  = sub2ind(sznew,subjnewI(1,:),subjnewI(2,:),subjnewI(3,:))'; 
% Outs
AlignedCoordinates = subjnewI';
AlignedIndices     = indjnew;
Origin             = subjori';
BoxSize            = sznew;
clear('subjnewI','indjnew','subjori');
end % GetAlignedAxon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AlignedIndices,AlignedValues] = ...
        GetAlignedAxonBetter(OriginalVolume,BoxSizeNew,OriginOld,OriginNew,RotMat)
% In:
% OriginalVolume (3D volume)
% BoxSizeNew (1 x 3)
% OriginOld (1 x 3)
% OriginNew (1 x 3)
% RotMat (3 x 3)
%
% Out:
% AlignedIndices (N x 1)

% Box size original volume
BoxSizeOld         = size(OriginalVolume);
if(length(BoxSizeOld) < 3)
Temp = [1 1 1];
Temp(1:length(BoxSizeOld)) = BoxSizeOld;
BoxSizeOld = Temp;
end % 3d
% get coordinates and indices of the full rotated volume
[X,Y,Z]            = meshgrid(1:BoxSizeNew(1),1:BoxSizeNew(2),1:BoxSizeNew(3));
subjnewJ           = [X(:) Y(:) Z(:)]';
clear('X','Y','Z');
indjnewJ           = sub2ind(BoxSizeNew,subjnewJ(1,:),subjnewJ(2,:),subjnewJ(3,:)); 
% Remove translation of the volume coordinates
subjtra            = subjnewJ + OriginNew'; 
clear('subjnewJ')
% Remove rotation of the volume coordinates
subjold            = RotMat * subjtra ; 
clear('subjtra')
% Move volume coordinates to origin
subj2              = round( subjold - OriginOld' );
clear('subjold')
% get indices of voxels inside bounds (everything else is zero)
flagvoxin           = (1 <= subj2(1,:)) & (subj2(1,:) <= BoxSizeOld(1)) & ...
                      (1 <= subj2(2,:)) & (subj2(2,:) <= BoxSizeOld(2)) & ...
                      (1 <= subj2(3,:)) & (subj2(3,:) <= BoxSizeOld(3));
subj2(:,~flagvoxin) = 1; % for the next step not to fail, but these voxels do not matter (will be set to zero)          
indj2               = sub2ind(BoxSizeOld,subj2(1,:),subj2(2,:),subj2(3,:)); 
clear('subj2')
% get values from original (not rotated) volume (same as NN interpolation)
valsin              = squeeze( OriginalVolume( indj2 ) );
clear('OriginalVolume','indj2')
valsinN             = valsin(:)';
clear('valsin')
% get final indices of the rotated axon
if(nargout == 1)
finalflags          = flagvoxin & valsinN;    
clear('flagvoxin','valsinN')
AlignedIndices      = indjnewJ(finalflags)';
clear('finalflags')
else
AlignedIndices      = indjnewJ(flagvoxin)';   
AlignedValues       = valsinN(flagvoxin)';
clear('flagvoxin','valsinN')
end % outputs
clear('indjnewJ')
end % GetAlignedAxonBetter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AlignedIndices,AlignedValues1,AlignedValues2] = ...
        GetAlignedAxonBetterTwoVols(OriginalVolume1,OriginalVolume2,BoxSizeNew,OriginOld,OriginNew,RotMat)
% In:
% OriginalVolume1 (3D volume)
% OriginalVolume2 (3D volume) - same size as 1
% BoxSizeNew (1 x 3)
% OriginOld (1 x 3)
% OriginNew (1 x 3)
% RotMat (3 x 3)
%
% Out:
% AlignedIndices (N x 1)

% Box size original volume
BoxSizeOld         = size(OriginalVolume1);
if(length(BoxSizeOld) < 3)
Temp = [1 1 1];
Temp(1:length(BoxSizeOld)) = BoxSizeOld;
BoxSizeOld = Temp;
end % 3d
% 1. first the indices
% get coordinates and indices of the full rotated volume
[X,Y,Z]            = meshgrid(1:BoxSizeNew(1),1:BoxSizeNew(2),1:BoxSizeNew(3));
subjnewJ           = [X(:) Y(:) Z(:)]';
clear('X','Y','Z');
indjnewJ           = sub2ind(BoxSizeNew,subjnewJ(1,:),subjnewJ(2,:),subjnewJ(3,:)); 
% Remove translation of the volume coordinates
subjtra            = subjnewJ + OriginNew'; 
clear('subjnewJ')
% Remove rotation of the volume coordinates
subjold            = RotMat * subjtra ; 
clear('subjtra')
% Move volume coordinates to origin
subj2              = round( subjold - OriginOld' );
clear('subjold')
% get indices of voxels inside bounds (everything else is zero)
flagvoxin           = (1 <= subj2(1,:)) & (subj2(1,:) <= BoxSizeOld(1)) & ...
                      (1 <= subj2(2,:)) & (subj2(2,:) <= BoxSizeOld(2)) & ...
                      (1 <= subj2(3,:)) & (subj2(3,:) <= BoxSizeOld(3));
% get final indices of the rotated axon
AlignedIndices      = indjnewJ(flagvoxin)';  
clear('indjnewJ')
% 2. then the values
subj2(:,~flagvoxin) = 1; % for the next step not to fail, but these voxels do not matter (will be set to zero)          
indj2               = sub2ind(BoxSizeOld,subj2(1,:),subj2(2,:),subj2(3,:)); 
clear('subj2')
% get values from original (not rotated) volume (same as NN interpolation)
valsin1             = squeeze( OriginalVolume1( indj2 ) );
clear('OriginalVolume1')
valsin2             = squeeze( OriginalVolume2( indj2 ) );
clear('OriginalVolume2','indj2')
valsinN1            = valsin1(:)';
valsinN2            = valsin2(:)';
clear('valsin1','valsin2')
% get final values of the rotated axon
AlignedValues1      = valsinN1(flagvoxin)';
AlignedValues2      = valsinN2(flagvoxin)';
clear('flagvoxin','valsinN1','valsinN2')
end % GetAlignedAxonBetterTwoVols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ObjectFill = CrossSectionHoleFill(Object,zrange)
Size       = size(Object);    
ObjectFill = zeros(Size,'logical');
for k = zrange'
ObjectFill(:,:,k) = imfill(Object(:,:,k),'holes');
end % k
clear('Object','zrange')
end % CrossSectionHoleFill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WSmask = WatershedMask(Axon,WSlabels,pthr)
NewSize = size(Axon);
WSmask  = zeros(NewSize,'logical');
nvoxax  = sum(Axon(:));   
IL      = uint16(Axon).*WSlabels; 
clear('Axon')
ILu     = unique(IL(:));
% keep considerable percentage
for k = 1:length(ILu)
if(ILu(k) > 0)
ILfax   = IL == ILu(k);
pIL     = sum(ILfax(:))/nvoxax;
clear('ILfax')
if(pIL >= pthr)
ILfws   = WSlabels == ILu(k);    
WSmask(ILfws) = true;
clear('ILfws')
end % if percentage
end % if no zero
end % k
clear('WSlabels','IL')
end % WatershedMask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FiberKept = KeepLargestComponents(Fiber,pthr)
FiberKept = false(size(Fiber));
CC        = bwconncomp(Fiber);
clear('Fiber')
numVoxels = cellfun(@numel,CC.PixelIdxList);
perVoxels = numVoxels/sum(numVoxels);
for k = 1:CC.NumObjects
if( perVoxels(k) >= pthr )    
FiberKept(CC.PixelIdxList{k}) = true;
CC.PixelIdxList{k} = [];
end % if 
end % k
clear('CC')
end % KeepLargestComponents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FiberKept,numObj,vox1st,voxrest,prest1st] = ...
        KeepLargestComponentsMore(Fiber,pthr)
FiberKept  = false(size(Fiber));
CC         = bwconncomp(Fiber);
clear('Fiber')
% num voxels (percentage) per object
numVoxels  = cellfun(@numel,CC.PixelIdxList);
sumVoxels  = sum(numVoxels);
% Sort
[vs,is]    = sort(numVoxels,'descend');
% num objects
numObj     = CC.NumObjects;
% -------------------------------------- %
% num voxels bigger
vox1st     = vs(1);
% num voxels all kept
if(numObj > 1)
idxrest    = vs/sumVoxels >= pthr ;
idxrest(1) = 0; % bigger
voxrest    = sum( vs( idxrest ) );
else
voxrest    = 0.0;
end % if-else sectioned 
% proportion volumes bigger and all the rest
prest1st   = voxrest/vox1st; % proportion
% -------------------------------------- %
% keep few 
perVoxels = numVoxels/sumVoxels;
if(pthr <= 0 || 1 <= pthr)
% keep only the largest object
FiberKept(CC.PixelIdxList{is(1)}) = true;
else
% keep objects larger than pthr 
for k = 1:numObj
if( perVoxels(k) >= pthr )    
FiberKept(CC.PixelIdxList{k}) = true;
CC.PixelIdxList{k} = [];
end % if 
end % k
end % if-else
% -------------------------------------- %
clear('CC')
end % KeepLargestComponentsBetter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Axon = AnalyzeCrossSections(fiber,zrange0,minl,minr,minprop)
% start with closing operation
[nx,ny,nz] = size(fiber);
SE         = strel('sphere',minr);
fibercl   = imclose(fiber,SE); 
clear('fiber','SE')
% fibercl1   = imresize3(fibercl,4,'nearest') > 0;
% fibercl2   = imresize3(fibercl,4,'linear') > 0;
% zrange0    = find(sum(sum(fibercl,1),2));
% [nx,ny,nz] = size(fibercl);
% clear('fibercl1')
% ----------------------------------------------------------------------- %
fibern     = zeros(nx,ny,nz,'uint8');  % fill holes per section / break into objects  
ncs        = length(zrange0);  % number of cross sections
skels      = zeros(nz,3,4);    % skeletons (up to three per cross-section + average)
skelsn     = zeros(nz,3,4);    % skeletons normalized
nskelwin   = zeros(nz,1);      % skeleton winner (down->up)
nocso      = zeros(nz,1);      % number of objects per cross section original
nocsn      = zeros(nz,1);      % number of objects per cross section after corrections 
nvcso      = zeros(nz,1);      % total voxels per cross section original
nvcsn      = zeros(nz,1);      % total voxels per cross section (up to 3 objects)  
% ======================================================================= %
% fill cross sections and checking them for multiple objects 
for z = zrange0'
% extract, fill and decompose cross section    
fiber2d    = squeeze(squeeze(fibercl(:,:,z)));  % get cross section
fiber2df   = imfill(fiber2d,'holes');           % fill holes 
clear('fiber2d')
cc2d       = bwconncomp(fiber2df);              % get objects
nobj       = cc2d.NumObjects;                   % number of objects
nocso(z)   = nobj;
% ----------------------------------------------------------------------- %
% figure(1); subplot(1,2,1); imshow(fiber2df);
clear('fiber2df')
% remove smaller objects
if(nobj > 1)
nv         = cellfun(@numel,cc2d.PixelIdxList);
[nvs,ivs]  = sort(nv,'descend');
minvox     = round(minprop*sum(nvs)); % remove objects smaller than certain proportion of bigger object
% remove very small objects
nobjn = nobj;    
for i = 2:nobjn
if(nvs(i) < minvox)     
nobjn = nobjn - 1;
end % if very small
end % i
% removing medium size objects if there are more than three 
if(nobjn > 3)
nobjn = 3;
end % if more than 3
else
% only one object
nvs   = numel(cc2d.PixelIdxList{1});  
ivs   = 1;
nobjn = 1;
end % if-else more than 1
% ----------------------------------------------------------------------- %
% compute center of masses for each object
nocsn(z)        = nobjn;
fiber2dn        = zeros(nx,ny,'uint8');
for i = 1:nobjn
[I,J]           = ind2sub([nx,ny],cc2d.PixelIdxList{ivs(i)});
cm              = [mean(I),mean(J),z];
cmn             = normr(cm); 
skels(z,:,i+1)  = cm;
skelsn(z,:,i+1) = cmn;
fiber2dn(cc2d.PixelIdxList{ivs(i)}) = i;
clear('I','J','cm','cmn')
end % i
% compute center of mass considering all objects
[I,J]         = ind2sub([nx,ny],find(fiber2dn));
cm            = [mean(I),mean(J),z];
cmn           = normr(cm); 
skels(z,:,1)  = cm;
skelsn(z,:,1) = cmn;
clear('I','J','cm','cmn')
% figure(1); subplot(1,2,2); imshow(fiber2dn,[0 3]);
% ----------------------------------------------------------------------- %
% set new cross-section (splitted)
fibern(:,:,z) = fiber2dn;
% volumes
nvcso(z)      = sum(nvs);
nvcsn(z)      = sum(nvs(1:nobjn));
clear('cc2d','fiber2dn')
end % z
clear('fibercl')
% ======================================================================= %
% some flags indicating certain characteristics
skone    = nocsn == 1;       % cross-sections with only one object 
sknozero = nocsn >  0;       % cross-sections with at least one object
skmany   = nocsn >  1;       % cross-sections with more than one object
% check lengths of the parts of the axon with only one skeleton
[~,nobjl1,~,locobj1] = RefineSeg.CheckLengthObjects(skone,minl);
% check number of cross sections with more than 1 object
nobjm  = sum( skmany );
pobjl1 = sum(locobj1)/ncs;
% ======================================================================= %
% ======================================= %
% ======= Select only one skeleton ====== %
% ======================================= %
% fskel (0: bad: nothing done, 1: good: only one, 2: bad: but corrected)
% fskel = -1;
skel  = zeros(nz,3);
if( nobjm == 0 )
% only one skeleton
fskel = 1;
% fprintf('One Skeleton: Nothing\n')
nskelwin(sknozero) = 2;  % 1st (only one object)
skel               = squeeze( skels(:,:,2) );
else
% more than one skeleton
% fprintf('Multiple Skeletons\n')
if(nobjl1 > 0)
% keep only one skeleton
fskel = 2;
% fprintf('Good: Keep one\n')    
nskelwin(skone) = 2; % 1st (only one object)
skel(skone,:)   = squeeze( skels(skone,:,2) );
% --- down-up --- % 
cmn   = [];
for k = 1:nz
% center of mass of this cross section if there is only one 
if(locobj1(k))     
cmn   = skelsn(k,:,2);
end % if one
% is this region one containing multiple objects
if(skmany(k))
% correct only if there is a cm of reference    
if(~isempty(cmn))
cmsn  = squeeze( skelsn(k,1:3,1:(nocsn(k)+1)) )';
pps   = abs( sum( cmsn .* cmn , 2) );
[~,i] = max(pps);
cmn         = cmsn(i,:);
nskelwin(k) = i;
skel(k,:)   = squeeze( skels(k,:,i) );
end % if non empty
end % if many
end % k 
% --- up-down --- % 
cmn   = [];
for k = nz:-1:1
% center of mass of this cross section if there is only one 
if(locobj1(k))    
cmn   = skelsn(k,:,2);
end % if one
% is this region one containing multiple objects
if(skmany(k))
% correct only if there is a cm of reference    
if(~isempty(cmn))
cmsn  = squeeze( skelsn(k,1:3,1:(nocsn(k)+1)) )';
pps   = abs( sum( cmsn .* cmn , 2) );
[~,i] = max(pps);
%
if(nskelwin(k) == 0)
% set if zero
cmn         = cmsn(i,:);
nskelwin(k) = i;
skel(k,:)   = squeeze( skels(k,:,i) );
else
% compare with previous one
if(nskelwin(k) ~= i)
% if different set to average   
cmn         = cmsn(1,:);
nskelwin(k) = 1;   
skel(k,:)   = squeeze( skels(k,:,1) );
end % different 
end % if-else set or not set
%
end % if non empty
end % if many
clear('cmsn','pps')
end % k 
else
% keep everything    
fskel              = 0;
% fprintf('Not long enough regions with one skeleton: keeping all with no modification\n')        
nskelwin(sknozero) = 1; % all objects
skel               = squeeze( skels(:,:,1) );
end
end % if-else one vs many
clear('skone','sknozero','skmany')
% ======================================= %
% ======================================= %
% ======================================= %
% ======================================================================= %
% % Plot original/unique_regions/selected
% propmobj = 100*sum(nocsn > 1)/ncs;
% % all
% figure; 
% subplot(1,3,1); hold on; 
% scatter3(skels(:,1,1),skels(:,2,1),skels(:,3,1),10,'r','filled'); 
% scatter3(skels(:,1,2),skels(:,2,2),skels(:,3,2),10,'g','filled'); 
% scatter3(skels(:,1,3),skels(:,2,3),skels(:,3,3),10,'b','filled'); 
% scatter3(skels(:,1,4),skels(:,2,4),skels(:,3,4),10,'k','filled'); 
% axis('equal');
% title(num2str(propmobj))
% view([0 0])
% % just regions with one
% subplot(1,3,2); hold on; 
% scatter3(skels(locobj1,1,1),skels(locobj1,2,1),skels(locobj1,3,1),10,'r','filled'); 
% scatter3(skels(locobj1,1,2),skels(locobj1,2,2),skels(locobj1,3,2),10,'g','filled'); 
% scatter3(skels(locobj1,1,3),skels(locobj1,2,3),skels(locobj1,3,3),10,'b','filled'); 
% scatter3(skels(locobj1,1,4),skels(locobj1,2,4),skels(locobj1,3,4),10,'k','filled'); 
% axis('equal');
% title(num2str(propmobj))
% view([0 0])
% % selected
% subplot(1,3,3); hold on; 
% scatter3(skel(:,1),skel(:,2),skel(:,3),10,'k','filled'); 
% axis('equal');
% title(num2str(propmobj))
% view([0 0])
% ======================================================================= %
clear('skels','skelsn')
% clear('locobj1')
% set fiber / num obj / volumes 
fibersel = zeros(nx,ny,nz,'logical');
nocss    = zeros(nz,1); % number of objects
nvcss    = zeros(nz,1); % number of voxels
for z = zrange0'
fiber2d  = fibern(:,:,z);    
if(nskelwin(z) == 1)
% all objects kept    
fibersel(:,:,z) = fiber2d > 0;
nvcss(z)        = nvcsn(z); 
nocss(z)        = nocsn(z);
else
fibersel(:,:,z) = fiber2d == (nskelwin(z)-1);
nvcss(z)        = sum(fibersel(:,:,z),'all'); 
nocss(z)        = 1;
end
clear('fiber2d')
end % z
clear('fibern','nvcsn','nocsn','nskelwin')
% ======================================================================= %
% check lengths of the parts of the axon with narrow area
fewvox = nvcss < 2*minr+1; % cross-sections with narrow areas (or zero)
[~,nobjln,~,locobjn] = RefineSeg.CheckLengthObjects(fewvox,minl);
nvoxn  = sum(fewvox);
pobjln = sum(locobjn)/ncs;
clear('locobjn')
% check lengths of the parts of the axon with more than one object
locmany = nocss > 1;
[~,nobjlm,~,locobjm] = RefineSeg.CheckLengthObjects(locmany,minl);
nvoxm  = sum(locmany);
pobjlm = sum(locobjm)/ncs;
clear('locobjm')
% check lengths of the parts of the axon with only one skeleton
skone    = nocss == 1;       % cross-sections with only one object 
[~,nobjl2,~,locobj2] = RefineSeg.CheckLengthObjects(skone,minl);
pobjl2 = sum(locobj2)/ncs;
clear('skone')
% ======================================================================= %
% connected components
CC   = bwconncomp(fibersel);
nv   = cellfun(@numel,CC.PixelIdxList);
snv  = sum(nv);
nvs  = sort(nv,'descend');
nobj = CC.NumObjects;
if(nobj == 1)
p1st = 1;
else
p1st = nvs(1)/snv;    
end
% ======================================================================= %
% Compute axon euclidean and geodesic lengths
[leuclidean,lgeodesic,sinuosity] = RefineSeg.FiberLengths(skel);
% Compute cross section areas along axon
[csareasalign,csareas] = RefineSeg.CrossSectionAreas(fibersel,skel); 
% figure; plot(1:nz,csareas,'-k',1:nz,nvcss,'-r'); % they should be the same
% ======================================================================= %
VoxelIdxList = find(fibersel);
% ======================================================================= %
% Outs
Axon.AxonVoxelIdxList = VoxelIdxList; % axon indices
Axon.size             = [nx ny nz];   % size
Axon.skeleton         = skel;         % skeleton
Axon.volume           = snv;          % final volume axon
Axon.euclideanlen     = leuclidean;   % euclidean length    (exluding empty spaces)
Axon.geodesiclen      = lgeodesic;    % geodesic length     (exluding empty spaces)
Axon.sinuosity        = sinuosity;    % sinuosity     
Axon.csareas          = csareas;      % cross section areas (exluding empty spaces)
Axon.csareasalign     = csareasalign; % cross section areas (exluding empty spaces)
%
Axon.numobjects       = nobj;     % final axon objects
Axon.vollargest       = nvs(1);   % volume largest
Axon.proplargest      = p1st;     % proportion largest object
Axon.skflag           = fskel;    % flag about proc skeleton  (0: bad, nothing done, 1: good, only one, 2: regular, many but corrected)
% Axon.nvoxcrosssecorig = nvcso;    % number voxels per cross section original
% Axon.nobjcrosssecorig = nocso;    % number objects per cross section original
Axon.nvoxcrosssec     = uint32(nvcss); % number voxels per cross section corrected
Axon.nobjcrosssec     = uint32(nocss); % number objects per cross section corrected
Axon.loclongoneskel   = locobj1;  % location of long regions where there is only one skeleton (before final selection)
Axon.numlongoneskel   = nobjl1;   % number of long regions where there is only one skeleton (before final selection) 
Axon.proplongoneskel  = pobjl1;   % length proportion of long regions where there is only one skeleton (before final selection) 
Axon.loclongoneskel2  = locobj2;  % location of long regions where there is only one skeleton (before final selection)
Axon.numlongoneskel2  = nobjl2;   % number of long regions where there is only one skeleton (before final selection) 
Axon.proplongoneskel2 = pobjl2;   % length proportion of long regions where there is only one skeleton (before final selection) 
Axon.numnarrow        = nvoxn;    % number of cross sections with narrow areas (after final selection) 
Axon.locnarrow        = fewvox;   % locations of cross sections with narrow areas (after final selection) 
Axon.numlongnarrow    = nobjln;   % number of long regions where there is narrow area (after final selection) 
Axon.proplongnarrow   = pobjln;   % length proportion of long regions where there is narrow area (after final selection) 
Axon.nummanyskel      = nvoxm;    % number of cross sections with more than one skeleton (after final selection)
Axon.locmanyskel      = locmany;  % locations of cross sections with more than one skeleton (after final selection)
Axon.numlongmanyskel  = nobjlm;   % number of long regions where there is more than one skeleton (after final selection) 
Axon.proplongmanyskel = pobjlm;   % length proportion of long regions where there is more than one skeleton (after final selection) 
clear('fibersel','VoxelIdxList','skel','zrange0','CC')
clear('nvcso','nvcsn','nvcss','nocso','nocsn','nocss')
clear('locobj1','locobj2','nv','nvs','csareasalign','csareas')
clear('fewvox','locmany')
end % CheckCrossSections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nobj,nobjl,lenobj,locobj] = CheckLengthObjects(FlagsPerCrossSec,minlen)
% flag per cross sections 
nsf    = FlagsPerCrossSec > 0; % flag objects
locobj = false(size(nsf));
if(sum(nsf) == 0)
% empty
nobj   = 0;
nobjl  = 0;
lenobj = 0;
else
cc1d   = bwconncomp(nsf);
nobj   = cc1d.NumObjects;
lenobj = cellfun(@numel,cc1d.PixelIdxList);
nobjl  = 0;
for i = 1:nobj
if(lenobj(i) >= minlen) % larger than 1um
locobj(cc1d.PixelIdxList{i}) = 1;
nobjl  = nobjl + 1;
end % if
end % i
end % if-else
clear('FlagsPerCrossSec','nsf','cc1d')
end % CheckLengthObjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is wrong and is the function used for results White_matter_EM_Diam_Skel_Props
% Not much problem just use the skeletons to compute again but with New function
function [leuclidean,lgeodesic,sinuosity] = FiberLengths(skel)
% regions with no skeleton 
fskel   = sum(skel,2) == 0;
skelnan = skel;
skelnan(fskel,:) = nan;
% compute differences 
df   = diff(skelnan,1,1);
% compute differential distances
dfn  = sqrt(sum(df.*df,2));
% keep those not affected by nans (finite)
ffin = ~isnan(dfn);
% geodesic length 
lgeodesic  = sum(dfn(ffin));
% euclidean length 
leuclidean = sum(ffin); % <--- this is the part that its wrong 
% sinuosity
sinuosity  = lgeodesic/leuclidean;
clear('skel','skelnan','df','dfn','fskel','ffin')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [leuclidean,lgeodesic,sinuosity] = FiberLengthsNew(skel)
% regions with no skeleton 
fskel   = sum(skel,2) == 0;
skelnan = skel;
skelnan(fskel,:) = nan;
% compute differences 
df   = diff(skelnan,1,1);
% compute differential distances
dfn  = sqrt(sum(df.*df,2));
% keep those not affected by nans (finite)
ffin = ~isnan(dfn);
% geodesic length 
lgeodesic  = sum(dfn(ffin));
% euclidean length 
ske        = skelnan(~fskel,:);
dfe        = ske(end,:) - ske(1,:);
leuclidean = sqrt(sum(dfe.*dfe));
% sinuosity
sinuosity  = lgeodesic/leuclidean;
clear('skel','skelnan','df','dfn','fskel','ffin','ske','dfe')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [csareasaligned,csareas] = CrossSectionAreas(fiber,skel)
[nx,ny,nz] = size(fiber);    
% regions with no skeleton 
fskel   = sum(skel,2) == 0;
skelnan = skel;
skelnan(fskel,:) = nan;
clear('fskel')
% compute normals to cross sections
tg    = gradient(skelnan')'; 
tgn   = normr(tg);
tgnf  = sum(tg,2); 
clear('skelnan','tg')
% keep those not affected by nans (finite)
tfin  = ~isnan(tgnf);
nfin  = sum(tfin);
tgns  = tgn(tfin,:);
skels = skel(tfin,:);
axon  = fiber(:,:,tfin);
clear('skel','tgn','tgnf')
% rotation matrices
zax   = repmat([0 0 1],[nfin 1]);
rmat  = RefineSeg.FromToRotationMatrices(zax,tgns);
clear('zax')
% clear('tgns')
% base cross section
nx2   = ceil(nx/2);
ny2   = ceil(ny/2);
[X,Y] = meshgrid(-nx2:0.5:nx2,-ny2:0.5:ny2); 
xy    = [X(:) Y(:) zeros(numel(X),1)];
clear('X','Y')
% ======================================================================= %
% Compute cross-section areas
nv = zeros(nfin,1); % not aligned
na = zeros(nfin,1); % aligned
for i = 1:nfin
% % === HHLee method is slower and results are very similar === %  
% fiberi = bwalign(axon,skels(i,:),tgns(i,:),[nx ny 1]); 
% nv(i)  = sum(fiberi,'all');
% === my method is less elegant but is faster and results are similar === %
% the noisy peaks seem to be not because of the alignment but because of  
% small errors in segmentation which causes discontinuities in the skeletons 
% 1. rotate
rxy   = rmat(:,:,i) * xy';
% sum(rxy' * tgns(i,:)') % should be zero (each entry should be zero)
% 2. translate
rxyt  = round( rxy + skels(i,:)' );
% 3. get indices of voxels inside bounds (everything else is zero)
fvin  = (1 <= rxyt(1,:)) & (rxyt(1,:) <= nx) & ...
        (1 <= rxyt(2,:)) & (rxyt(2,:) <= ny) & ...
        (1 <= rxyt(3,:)) & (rxyt(3,:) <= nz);
indrt = sub2ind([nx,ny,nz],rxyt(1,fvin),rxyt(2,fvin),rxyt(3,fvin)); 
% 4. get values from original (not rotated) volume (same as NN interpolation)
valin = squeeze( fiber( indrt ) );
nv(i) = sum(valin)/4; % factor 4 because of the 0.5 spacing on x and y
na(i) = sum(axon(:,:,i),'all');
%
% % show cross section
% vals  = fvin; vals(fvin) = valin;
% figure(100); scatter3(xy(:,2) ,xy(:,1) ,xy(:,3) ,10,vals,'filled'); view([0 90]); title(num2str(i))
%
% % plot plane / normal vector / rotated plane
% figure; hold on
% scatter3(xy(:,1) ,xy(:,2) ,xy(:,3) ,10,'k','filled');
% scatter3(rxy(1,:),rxy(2,:),rxy(3,:),10,'r','filled');
% plot3(10*[0 tgns(i,1)],10*[0 tgns(i,2)],10*[0 tgns(i,3)],'-g','linewidth',5);
% plot rotates translated plane / skeletons
% figure; hold on 
% scatter3(skels(:,1),skels(:,2),skels(:,3),10,'r','filled');
% scatter3(rxyt(1,:),rxyt(2,:),rxyt(3,:),10,'k','filled');
clear('rxy','rxyt','indrt','valin','fvin')
end % i
clear('xy','fiber','axon','skels','rmat','tgns')
% ======================================================================= %
% Outputs
csareasaligned = zeros(nz,1);
csareas        = zeros(nz,1);
csareasaligned(tfin) = nv;
csareas(tfin)        = na;
% figure; plot(1:nz,csareas,'-k',1:nz,csareasaligned,'-r')
clear('na','nv','tfin')
end %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rot = FromToRotationMatrices(From,To)
% Rot = FromToRotationMatrices(From,To)
% Create Rotation Matrices Rot in a way that
% Rot(:,:,i)*From(i,:) = To(i,:)
%
% Weird results when the angle between From(i) and To(i) is close to 180.
% Maybe because the cross product fail to get the correct normal to the
% plane since the solutions are infinite if From(i) and To(i) are almost
% parallel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify inputs
sF = size(From);
sT = size(To);
if(~isequal(sF,sT))
    error('From and To arrays should have the same size');
end
if(sF(2) ~= 3)
    error('The size of the arrays should be N x 3');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating matrices
N    = sF(1);
Rot  = zeros(3,3,N);
% unitary cross product between initial and final orientation
% normal to the plane including both directions
CP   = cross(From,To,2);  
rcp  = sqrt(sum(CP.*CP,2));
CPu  = CP./[rcp rcp rcp];
% angle between initial and final orientation
PP   = sum(From.*To,2);
rf   = sqrt(sum(From.*From,2));
rt   = sqrt(sum(To.*To,2));
rads = acos(PP./(rf.*rt));
% useful expressions
Cos  = cos(rads);          Sen  = sin(rads);          One  = ones(size(Cos));
Ux   = CPu(:,1);           Uy   = CPu(:,2);           Uz   = CPu(:,3);
Ux2  = CPu(:,1).*CPu(:,1); Uy2  = CPu(:,2).*CPu(:,2); Uz2  = CPu(:,3).*CPu(:,3);
Uxy  = CPu(:,1).*CPu(:,2); Uxz  = CPu(:,1).*CPu(:,3); Uyz  = CPu(:,2).*CPu(:,3);

Rot(1,1,:) = Cos + Ux2.*(One - Cos);
Rot(1,2,:) = Uxy.*(One - Cos) - Uz.*Sen;
Rot(1,3,:) = Uxz.*(One - Cos) + Uy.*Sen;
   
Rot(2,1,:) = Uxy.*(One - Cos) + Uz.*Sen;
Rot(2,2,:) = Cos + Uy2.*(One - Cos);
Rot(2,3,:) = Uyz.*(One - Cos) - Ux.*Sen;
   
Rot(3,1,:) = Uxz.*(One - Cos) - Uy.*Sen; 
Rot(3,2,:) = Uyz.*(One - Cos) + Ux.*Sen; 
Rot(3,3,:) = Cos + Uz2.*(One - Cos);
end
% figure; hold on;
% plot3([0 From(1)],[0 From(2)],[0 From(3)],'-k')
% plot3([0 To(1)],[0 To(2)],[0 To(3)],'-b')
% plot3([0 CPu(1)],[0 CPu(2)],[0 CPu(3)],'-r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function substrate = SetSubstrate(ias,myelin,cells)
substrate = ones(size(ias),'uint8');
substrate( myelin ) = 0;
clear('myelin')
if(~isempty(cells))
substrate( cells  ) = 3;
clear('cells')
end % if cell
substrate( ias    ) = 2; 
clear('ias')
end % SetSubstrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function substrate = SetCleanSubstrate(ias,myelin,cells,voxelsize,dilrad,diam,nvoxels)
% volume to set values
substrate = ones(size(ias),'uint8');
% ----------------------------------------------------------------------- %
% parameters to dilate ias 
radp = round( dilrad/voxelsize );
str1 = strel('sphere',radp);
% parameters to dilate ias and close myelin
str2 = strel('sphere',1);
% parameters to close holes in fibers
str3 = strel('sphere',2);
% parameters to remove small objects
vol  = (4/3)*pi*((diam/2)^3);
vp   = round(vol/(voxelsize^3));
% ----------------------------------------------------------------------- %
% remove tiny isolated objects 
fprintf('Remove Tiny Isolated Objects  ... ')
rems          = bwareaopen(ias    ,nvoxels); 
ias(~rems)    = 0;
clear('rems')
toc
remm          = bwareaopen(myelin ,nvoxels); 
myelin(~remm) = 0;
clear('remm')
toc
% ----------------------------------------------------------------------- %
% dilate ias a lot and combine with myelin to decide myelin that remains 
fprintf('Myelin + Highly Dilated IAS  ... ')
iasd = imdilate(ias,str1);
myd  = ((myelin) & iasd);
clear('myelin','iasd')
toc
% slightly dilated ias and close pixel holes in myelin
fprintf('Dilate IAS Slightly ... ')
iass = imdilate(ias,str2);
clear('ias')
toc
fprintf('Close Holes ... ')
mys  = imclose(myd,str2);
clear('myd')
% combine ias with myelin (fiber)
fs   = iass | mys;
clear('mys')
toc
% remove small objects
fprintf('Remove Small Objects ... ')
fso  = bwareaopen(fs,vp);
clear('fs')
toc
% close holes in fibers
fprintf('Close Holes Again ... ')
fsc  = imclose(fso,str3);
substrate( fsc  ) = 0;
clear('fsc')
toc
if(~isempty(cells))
fprintf('Cells ... ')
remc         = bwareaopen(cells  ,nvoxels); 
cells(~remc) = 0;
cellsf       = imfill(cells,'holes');
clear('cells')
substrate( cellsf ) = 3;  
clear('remc','cellsf')
toc
end % if cell
substrate( iass ) = 2; 
clear('iass')
end % SetCleanSubstrate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % static methods
end % classdef   
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx %
% AnalyzeCrossSections()
% % smooth fiber along z (I DO NOT SEE ANY ADVANTAGE)
% fibercls   = zeros(nx,ny,nz,'uint8');  
% w          = 5;
% for z = (w+1):(nz-w) % dont account for extremes
% fibercls(:,:,z) = sum(fibercl(:,:,(z-w):(z+w)),3);
% end % 
% fiberclsl          = fibercls > 1;%round(0.6*(2*w+1));
% fiberclsl(:,:,1:w) = fibercl(:,:,1:w);
% fiberclsl(:,:,(nz-w+1):nz) = fibercl(:,:,(nz-w+1):nz);
% clear('fibercl');
% clear('fibercls');
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx %
