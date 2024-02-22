function script11_plot_fods()  
% plot skeletons and FOD
%
% This script needs the results from script:
% script06_axon_diam_skel_props.m
% script09_fods_dispersion.m
% and the information in the .mat files in the folder
% cc_cg_labels/
% The results of this script will be saved in folder
% White_matter_EM_06_Plot_FODs/*/
% 26-May-2023 by Ricardo Coronado-Leija
close all
clear
clc
rng(10) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Begin\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ir         = 2;
res        = {'HM','LM'};
rootProps  = ['White_matter_EM_Diam_Skel_Props_' res{ir} '/'];
rootDisp   = ['White_Matter_EM_FOD_' res{ir} '/'];
rootQuant  = 'cc_cg_labels/'; 
rootOut    = ['White_matter_EM_06_Plot_FODs_' res{ir} '/'];
if(~isfolder(rootOut)); mkdir(rootOut); end
samples    = {'Sham_25_contra','Sham_25_ipsi','Sham_49_contra','Sham_49_ipsi',...
              'TBI_24_contra','TBI_24_ipsi','TBI_28_contra','TBI_28_ipsi',...
              'TBI_2_contra','TBI_2_ipsi'}; 
idsample   = [25 25 49 49 24 24 28 28 2 2 ];
sidesample = {'contra','ipsi','contra','ipsi','contra','ipsi',...
              'contra','ipsi','contra','ipsi'};
nsamples   = length(samples);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voxsizes   = [15 50]*1e-3;
vs         = voxsizes(ir);    % voxel size in micro-meters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Analysis dispersion 
flagDisp  = rand > 0; 
% Number of axons
naxsel    = [36561 37335 27529 13840 39477 26320 33612 31103 31732 23914];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select the samples to loop
fidx = 1:nsamples;  
% fidx = [2 6];
% fidx = [2 6 1 3 4 5 7:nsamples]; 
frcc = 1;
frcg = 1;
as   = analyzeseg_modified(); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(flagDisp)   
% ----------------------------------------------------------------------- % 
maxPmany = 0.25;             % proporion of cross sections with multiple skeletons
lmin     = 1.0;              % 1 um at extremes is removed
nmin     = round(lmin/vs);   % " in voxels 
lrem     = 5.0;              % min length of axons = 2*lrem 
nrem     = round(lrem/vs)-1; % " in voxels
sigmas   = 0;                % smoothing kernel width (pixels) 
ts       = 11.5;             % for FOD SH
maxangle = 45;               % remove outliers
maxpp    = cos(pi*maxangle/180);
if(ir == 2)
perc = 0.01;
else
perc = 0.1;
end
% ----------------------------------------------------------------------- % 
for f = fidx   
sample   = [num2str(idsample(f)) '_' sidesample{f}];  
% read labels cc cg (big low res volumes)
if(ir == 2)
lblsname = [rootQuant '/LM_' sample '_cc_cg_lbls.mat'];
lblscccg = load(lblsname);
flagcc   = ~isempty(lblscccg.cc_lbl) && frcc;
flagcg   = ~isempty(lblscccg.cg_lbl) && frcg;
else
% small high res volumes only contain cc    
flagcc   = 1 && frcc;
flagcg   = 0 && frcg;
end % if
% read axon skeletons 
skelname = [rootProps '/' samples{f} '/' res{ir} '_' sample '_axons_skeletons.mat'];
skel     = load(skelname).skeletons;
nax      = size(skel,1);
% read other axon properties (because of tranformation matrix)
propname = [rootProps '/' samples{f} '/' res{ir} '_' sample '_axons_props.mat'];
axprop   = load(propname).axon_props; % 1st colum of (axprop(1,45).MatRestore{:}*RotXtoZ) is axon orientation
% get axon nz length (tricky: some axons are empty, this modifies tab datatype from array to cell)
Nzo      = zeros(nax,1);
Pmany    = zeros(nax,1);
for j = 1:nax
if(isnumeric(axprop(j,3).euclideanlen) && isnumeric(axprop(j,22).nummanyskel))  
% euclidean length 
Nzo(j)   = axprop(j,3).euclideanlen;       
% proportion of cross sections with multiple skeletons
Pmany(j) = axprop(j,22).nummanyskel/Nzo(j);    
elseif(~isnumeric(axprop(j,3).euclideanlen) && ~isnumeric(axprop(j,22).nummanyskel))  
if(~isempty(axprop(j,3).euclideanlen{:}))    
% euclidean length 
Nzo(j)   = axprop(j,3).euclideanlen{:};       
% proportion of cross sections with multiple skeletons
Pmany(j) = axprop(j,22).nummanyskel{:}/Nzo(j);
end % non empty
else % other
Nzo(j)   = 0;
Pmany(j) = maxPmany;
end % numeric
end % j
% maxlength = max(Nzo);
% ----------------------------------------------------------------------- % 
if(flagcc)
% read fod
ccfod = load([rootDisp '/' samples{f} '_cc_fod_t' num2str(ts) '.mat']).fod;
% selecting axons only in cc region   
if(ir == 2)
cclbls     = lblscccg.cc_lbl(lblscccg.cc_lbl<=naxsel(f));
else
cclbls     = 1:nax;    
end
ccaxprop1  = axprop(cclbls,:); 
ccskel     = skel(cclbls,:);
ccNzo      = Nzo(cclbls);
ccPmany    = Pmany(cclbls);
clear('cclbls')
% remove very small axons and axons with many multiple skeletons per section 
% ccidx      = ccNzo > (2*nrem+2)/2; %ccNzo > 0; %ones(size(ccNzo)) > 0; %( ccNzo > (2*nrem+2) ) & ( ccPmany < maxPmany ); 
ccidx      = ( ccNzo > (2*nrem+2) ) & ( ccPmany < maxPmany ); 
clear('ccNzo','ccPmany');
ccaxidx    = [f 1]; 
ccaxprop   = ccaxprop1(ccidx,:); clear('ccaxprop1');
ccskeleton = ccskel(ccidx,:);    clear('ccskel'); clear('ccidx');
% compute fod and plot
[ccskSmooth,~,~] = SmoothedSkeletonsTangents(ccskeleton,ccaxprop,nmin,sigmas,ccaxidx);
clear('ccskeleton','ccaxprop')
% Plot FOD SH 
figure(1);    
[~,ccomax] = max(ccfod.fodi);
ccdir      = ccfod.points(ccomax,:);
as.fodsht(ccfod.fodi,ccfod.lmax2,'smooth','false','density','low','range',[0 10],'colorbar','off','glyph','on');
cameratoolbar('togglescenelight')
axis equal
if(f == 1)
view([-100 -15])
else
view([0 0])
end
namefod = [rootOut '/' res{ir} '_' samples{f} '_cc_FOD_SH_t' num2str(ts)];
pause(1)
SaveFigurePNG([namefod '.png'],50,50);
pause(1)
SaveFigureSVG([namefod '.svg'],50,50);
pause(1)
close(1)
% Plot skeletons
figure(2); hold on
ccnsk = length(ccskSmooth.fiber);
ccidx = randsample(ccnsk,round(perc*ccnsk))';
ccminlims = [];
ccmaxlims = [];
for i = ccidx
ccsk  = ccskSmooth.fiber(i).sk;
ccori = normr(ccsk(end,:) - ccsk(1,:));
ccmu0 = mean(ccsk,1);
ccsk0 = ccsk - ccmu0;   
if( maxpp < abs(sum(ccdir.*ccori)) ) % 0.9 < norm(ccori) &&
val = 0.3 + 0.7*rand; 
scatter3(ccsk0(:,1),ccsk0(:,2),ccsk0(:,3),1,val*abs(ccori),'filled'); %axis('equal');
ccminlims = cat(1,ccminlims,min(ccsk0,1));
ccmaxlims = cat(1,ccmaxlims,max(ccsk0,1));
end
end % i
cameratoolbar('togglescenelight')
if(f == 1)
view([-100 -15])
else
view([0 0])
end
axis equal off
set(gca,'Box','off')
xlim([min(ccminlims(:,1)) max(ccmaxlims(:,1))])
ylim([min(ccminlims(:,2)) max(ccmaxlims(:,2))])
zlim([min(ccminlims(:,3)) max(ccmaxlims(:,3))])
namesk = [rootOut '/' res{ir} '_' samples{f} '_cc_FOD_sk_t' num2str(ts)];
pause(1)
SaveFigurePNG([namesk '.png'],50,50);
pause(1)
SaveFigureSVG([namesk '.svg'],50,50);
pause(1)
close(2)
end % cc
% ----------------------------------------------------------------------- % 
if(flagcg)
% read fod
cgfod = load([rootDisp '/' samples{f} '_cg_fod_t' num2str(ts) '.mat']).fod;    
% selecting axons only in cg region   
if(ir == 2)
cglbls     = lblscccg.cg_lbl(lblscccg.cg_lbl<=naxsel(f));
else
cglbls     = [];    
end
cgaxprop1  = axprop(cglbls,:); 
cgskel     = skel(cglbls,:);
cgNzo      = Nzo(cglbls);
cgPmany    = Pmany(cglbls);
clear('cglbls')
% remove very small axons and axons with many multiple skeletons per section 
cgidx      = ( cgNzo > (2*nrem+2) ) & ( cgPmany < maxPmany ); 
clear('cgNzo','cgPmany');
cgaxidx    = [f 2]; 
cgaxprop   = cgaxprop1(cgidx,:); clear('cgaxprop1');
cgskeleton = cgskel(cgidx,:);    clear('cgskel'); clear('cgidx');
% compute fod and plot
[cgskSmooth,~,~] = SmoothedSkeletonsTangents(cgskeleton,cgaxprop,nmin,sigmas,cgaxidx);
clear('cgskeleton','cgaxprop')
% Plot FOD SH 
figure(1);    
[~,cgomax] = max(cgfod.fodi);
cgdir      = cgfod.points(cgomax,:);
as.fodsht(cgfod.fodi,cgfod.lmax2,'smooth','false','density','low','range',[0 10],'colorbar','off','glyph','on');
cameratoolbar('togglescenelight')
axis equal
if(f == 1)
view([-100 -15])
else
view([0 0])
end
namefod = [rootOut '/' res{ir} '_' samples{f} '_cg_FOD_SH_t' num2str(ts)];
pause(1)
SaveFigurePNG([namefod '.png'],50,50);
pause(1)
SaveFigureSVG([namefod '.svg'],50,50);
pause(1)
close(1)
% Plot skeletons
figure(2); hold on
cgnsk = length(cgskSmooth.fiber);
cgidx = randsample(cgnsk,round(perc*cgnsk))';
cgminlims = [];
cgmaxlims = [];
for i = cgidx
cgsk  = cgskSmooth.fiber(i).sk;
cgori = normr(cgsk(end,:) - cgsk(1,:));
cgmu0 = mean(cgsk,1);
cgsk0 = cgsk - cgmu0;   
if( maxpp < abs(sum(cgdir.*cgori)) ) % 0.9 < norm(cgori) && 
val = 0.3 + 0.7*rand; 
scatter3(cgsk0(:,1),cgsk0(:,2),cgsk0(:,3),1,val*abs(cgori),'filled'); axis('equal');
cgminlims = cat(1,cgminlims,min(cgsk0,1));
cgmaxlims = cat(1,cgmaxlims,max(cgsk0,1));
end
end % i
cameratoolbar('togglescenelight')
if(f == 1)
view([-100 -15])
else
view([0 0])
end
axis equal off
set(gca,'Box','off')
xlim([min(cgminlims(:,1)) max(cgmaxlims(:,1))])
ylim([min(cgminlims(:,2)) max(cgmaxlims(:,2))])
zlim([min(cgminlims(:,3)) max(cgmaxlims(:,3))])
namesk = [rootOut '/' res{ir} '_' samples{f} '_cg_FOD_sk_t' num2str(ts)];
pause(1)
SaveFigurePNG([namesk '.png'],50,50);
pause(1)
SaveFigureSVG([namesk '.svg'],50,50);
pause(1)
close(2)
end % cg
% ----------------------------------------------------------------------- % 
if(flagcc && flagcg)
% read fod
cccgfod = load([rootDisp '/' samples{f} '_cccg_fod_t' num2str(ts) '.mat']).fod;  
% Plot FOD SH 
figure(1);    
as.fodsht(cccgfod.fodi,cccgfod.lmax2,'smooth','false','density','low','range',[0 10],'colorbar','off','glyph','on');
cameratoolbar('togglescenelight')
axis equal
if(f == 1)
view([-100 -15])
else
view([0 0])
end
namefod = [rootOut '/' res{ir} '_' samples{f} '_cccg_FOD_SH_t' num2str(ts)];
pause(1)
SaveFigurePNG([namefod '.png'],50,50);
pause(1)
SaveFigureSVG([namefod '.svg'],50,50);
pause(1)
close(1)
% Plot skeletons
figure(2); hold on
cccgminlims = [];
cccgmaxlims = [];
for i = ccidx
ccsk  = ccskSmooth.fiber(i).sk;
ccori = normr(ccsk(end,:) - ccsk(1,:));
ccmu0 = mean(ccsk,1);
ccsk0 = ccsk - ccmu0;   
if( maxpp < abs(sum(ccdir.*ccori)) ) % 0.9 < norm(ccori) && 
val = 0.3 + 0.7*rand; 
scatter3(ccsk0(:,1),ccsk0(:,2),ccsk0(:,3),1,val*abs(ccori),'filled'); axis('equal');
cccgminlims = cat(1,cccgminlims,min(ccsk0,1));
cccgmaxlims = cat(1,cccgmaxlims,max(ccsk0,1));
end
end % i
for i = cgidx
cgsk  = cgskSmooth.fiber(i).sk;
cgori = normr(cgsk(end,:) - cgsk(1,:));
cgmu0 = mean(cgsk,1);
cgsk0 = cgsk - cgmu0;   
if( maxpp < abs(sum(cgdir.*cgori)) ) % 0.9 < norm(cgori) && 
val = 0.3 + 0.7*rand; 
scatter3(cgsk0(:,1),cgsk0(:,2),cgsk0(:,3),1,val*abs(cgori),'filled'); axis('equal');
cccgminlims = cat(1,cccgminlims,min(cgsk0,1));
cccgmaxlims = cat(1,cccgmaxlims,max(cgsk0,1));
end
end % i
cameratoolbar('togglescenelight')
if(f == 1)
view([-100 -15])
else
view([0 0])
end
axis equal off
set(gca,'Box','off')
xlim([min(cccgminlims(:,1)) max(cccgmaxlims(:,1))])
ylim([min(cccgminlims(:,2)) max(cccgmaxlims(:,2))])
zlim([min(cccgminlims(:,3)) max(cccgmaxlims(:,3))])
namesk = [rootOut '/' res{ir} '_' samples{f} '_cccg_FOD_sk_t' num2str(ts)];
pause(1)
SaveFigurePNG([namesk '.png'],50,50);
pause(1)
SaveFigureSVG([namesk '.svg'],50,50);
pause(1)
close(2)
end % cccg
clear('ccskSmooth')
clear('cgskSmooth')
% ----------------------------------------------------------------------- % 
clear('skel','Nzo','Pmany','axprop');
end % f
% ----------------------------------------------------------------------- % 
end % flagDisp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('End\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % main 

function [skSmooth,tgSmooth,ttSmooth] = SmoothedSkeletonsTangents(skeleton,axprop,nmin,sigmas,axidx)
% % rotation matrix x->z
% RotXtoZ = [0 0 -1; 0 1 0; 1 0 0];
thrnrm  = 0.01; % to remove zero norm tangents 
nax     = size(skeleton,1); % number of axons     
tic
for j = nax:-1:1 %2113% %336% %112%56 %5 % 1:nax
% (*) Alignment of skeleton to z-axis is not needed here (axon should be aligned already)
skori  = skeleton(j,:).skeleton{:};
% Removing extremes (cross sections are more difficult to estimate there)
skrem  = skori(nmin:end-nmin-1,:); 
clear('skori')
% Remove regions with no skeleton (this may create outliers, hope they are minimal)
fskel  = sum(skrem,2) == 0;   
sk     = skrem(~fskel,:);
clear('skrem','fskel')
% % Save the orientation of the axon (equivalent to tangent with maximum smooouthing)
% EigVec = axprop(1,45).MatRestore{:}*RotXtoZ; 
% tg     = repmat( EigVec(:,1)' , [size(sk,1)  1] ); % from restore rotation matrix
% tgInf.fiber(j).tg = tg;
% clear('tg')
% === Smoothing the skeletons and compute tangents === %
for s = 1:length(sigmas)
if(sigmas(s) <= 1e-6) 
% No smooth: 
cms = sk;
else
% Gaussian smoothing 
xq  = imgaussfilt(sk(:,1),sigmas(s)); % smooth x
yq  = imgaussfilt(sk(:,2),sigmas(s)); % smooth y
cms = [xq yq sk(:,3)];
clear('xq','yq');   
end % if-else
% Interpolate
cms     = InterpolateSkeleton(cms,0.25);
% Compute sinuosity
[~,~,t] = RefineSeg.FiberLengthsNew(cms);
% Restore smoothed skeleton (remove rotation that aligned it to z and translation)
if(isnumeric(axprop(j,43).OriginNew))
skr = ( axprop(j,45).MatRestore{:} * (cms' + axprop(j,43).OriginNew')  ); %     
else
skr = ( axprop(j,45).MatRestore{:} * (cms' + axprop(j,43).OriginNew{:}')  ); % 
end
clear('cms','sk')
% remove tangents with zero norm
% tgdiff = diff(skr,1,2)';%
tgdiff = gradient(skr)';
tgnrm  = sqrt(sum(tgdiff.*tgdiff,2));
tgnrmm = mean(tgnrm);  
tgidx  = abs(tgnrm) > thrnrm*tgnrmm; 
tgs    = normr(tgdiff(tgidx,:)); 
% skrt = skr'; 
% figure; scatter3(skrt(:,1),skrt(:,2),skrt(:,3),10,'k','filled'); axis('equal'); xlabel('x'); ylabel('y'); zlabel('z');
% figure; hold on; scatter3(tgs(:,1),tgs(:,2),tgs(:,3),10,'k','filled'); axis('equal'); xlabel('x'); ylabel('y'); zlabel('z');
% set to outputs
ttSmooth(s).tort(j)     = t;
tgSmooth(s).fiber(j).tg = tgs;
skSmooth(s).fiber(j).sk = skr';
clear('tgdiff','tgnrm','tgnrmm','tgidx','tgs','skr');
end % s
fprintf('%d %d %d\n',axidx(1),axidx(2),j);
% % check few smoothings
% i1 = 1; i2 = 3; i3 = 7;
% figure; scatter3(skSmooth(i1).fiber(j).sk(:,1),skSmooth(i1).fiber(j).sk(:,2),skSmooth(i1).fiber(j).sk(:,3),10,'k','filled'); axis('equal');
% figure; scatter3(skSmooth(i2).fiber(j).sk(:,1),skSmooth(i2).fiber(j).sk(:,2),skSmooth(i2).fiber(j).sk(:,3),10,'k','filled'); axis('equal');
% figure; scatter3(skSmooth(i3).fiber(j).sk(:,1),skSmooth(i3).fiber(j).sk(:,2),skSmooth(i3).fiber(j).sk(:,3),10,'k','filled'); axis('equal');
% figure; scatter3(tgSmooth(i1).fiber(j).tg(:,1),tgSmooth(i1).fiber(j).tg(:,2),tgSmooth(i1).fiber(j).tg(:,3),10,'k','filled'); axis('equal');
% figure; scatter3(tgSmooth(i2).fiber(j).tg(:,1),tgSmooth(i2).fiber(j).tg(:,2),tgSmooth(i2).fiber(j).tg(:,3),10,'k','filled'); axis('equal');
% figure; scatter3(tgSmooth(i3).fiber(j).tg(:,1),tgSmooth(i3).fiber(j).tg(:,2),tgSmooth(i3).fiber(j).tg(:,3),10,'k','filled'); axis('equal');
% close all
clear('sk');
end % j each axon
toc
clear('skeleton','axprop','axidx')
end % SmoothedSkeletonsTangents

function ski = InterpolateSkeleton(sk,dz)
x   = sk(:,1);
y   = sk(:,2);
z   = sk(:,3);
nzq = round( (z(end)-z(1))/dz );
zq  = linspace(z(1),z(end),nzq);
zq  = zq(:);
xq  = interp1(z,x,zq,'spline');
yq  = interp1(z,y,zq,'spline');
ski = [xq yq zq];
end
