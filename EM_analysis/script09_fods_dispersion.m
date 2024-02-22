function script09_fods_dispersion(samplei,regioni,timei)   
% Compute FOD and some of its properties for all axons in the EM sample
%
% samplei = the index of the EM sample to analyze
% regioni = 1 == cc only , 2 == cg only , 3 == cc+cg (combined) , 4 == cc and cg (separated at same run)
% timei   = diffusion time (ms)
%
% This script needs the results from script:
% script06_axon_diam_skel_props.m
% It also needs the functions in the class: 
% analyzeseg_modified.m (original file in github.com/NYU-DiffusionMRI/RaW-seg)
% and the information in the .mat files in the folder
% cc_cg_labels/
% The results of this script will be saved in folder
% White_Matter_EM_FOD/*/
% 25-May-2023 by Ricardo Coronado-Leija
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
rootDisp   = ['White_Matter_EM_FOD_' res{ir} '/']; 
if(~isfolder(rootDisp)); mkdir(rootDisp); end
rootQuant  = [root 'White_Matter_EM_Quant/']; 
rootProps  = [root 'White_matter_EM_Diam_Skel_Props_' res{ir} '/'];
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
if(nargin > 0)
fidx = samplei;
tidx = timei;
if(regioni == 1)
fcmb = 0;    
frcc = 1;
frcg = 0;
elseif(regioni == 2)
fcmb = 0;    
frcc = 0;
frcg = 1;
elseif(regioni == 3)
fcmb = 1;    
frcc = 1;
frcg = 1;
elseif(regioni == 4)
fcmb = 0;    
frcc = 1;
frcg = 1;
else
fprintf('regioni not valid');
end
else
fidx = 1:nsamples;
% tidx = [0 1 12 25 50 75 100 250 500 1000]; 
tidx = [0 11.5];          
fcmb = 0;
frcc = 1;
frcg = 1;
end % if
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(flagDisp)   
% ----------------------------------------------------------------------- % 
lmax1    = 6; 
lmax2    = 16;
maxPmany = 0.25;                  % proporion of cross sections with multiple skeletons
lmin     = 1.0;                   % 1 um at extremes is removed
nmin     = round(lmin/vs);        % " in voxels 
lrem     = 5.0; % 1.0; %          % min length of axons = 2*lrem 
nrem     = round(lrem/vs)-1;      % " in voxels
ts       = tidx;                  % diffusion times (ms)
D        = 2;                     % intrinsic diffusivity (micron^2/ms)
Sigmas   = sqrt(2*D*ts)/2;        % smoothing kernel width (micron), sigma = L/2 -> L=sqrt(2*D*t)
sigmas   = Sigmas/vs;             % smoothing kernel width (pixels) 
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
% ----------------------------------------------------------------------- % 
if(flagcc)
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
%ccskSmooth
[~,cctgSmooth,ccttSmooth] = SmoothedSkeletonsTangents(ccskeleton,ccaxprop,nmin,sigmas,ccaxidx);
clear('ccskeleton','ccaxprop')
% concatenate
cctangent = ConcatenateTangents(cctgSmooth);
% ----------------------------------------------------------------------- %
% % saving
% for iii = 1:length(ccskSmooth)
% t_D_vs_s = [ts(iii) D vs Sigmas(iii)];    
% sksmooth = ccskSmooth(iii);
% tgsmooth = cctgSmooth(iii);
% save([rootDisp '/' samples{f} '_cc_smooth_skeletons_t' num2str(ts(iii),'%.1f') '.mat'],'sksmooth','t_D_vs_s');
% save([rootDisp '/' samples{f} '_cc_smooth_tangents_t'  num2str(ts(iii),'%.1f') '.mat'],'tgsmooth','t_D_vs_s');
% end % iii  
% ----------------------------------------------------------------------- %
clear('ccskSmooth','cctgSmooth')
if(~fcmb)  
% fods
ccfods = ComputeFODs(cctangent,lmax1,lmax2);
clear('cctangent')
% save images 
k     = 1;
lmaxs = [lmax1 lmax2];
for iii = 1:length(ccfods)
for lll = 1:2    
figure(k)
pause(1)
title(['CC ' samples{f} ' lmax = ' num2str(lmaxs(lll))  ' t = ' num2str(ts(iii)) ' ms'], ...
    'fontsize',30,'interpreter','none')
pause(1)
namefig = [rootDisp '/' samples{f} '_cc_FOD_SH_lmax' num2str(lmaxs(lll)) '_t' num2str(ts(iii))];
SaveFigurePNG([namefig '.png'],16*3,9*3);
savefig(      [namefig '.fig'])
pause(1)
close(k); k = k + 1;
end % lll
end % iii
% plots
k = 1;
for iii = 1:length(ccfods)
% lmax1
figure(k); 
set(gca,'yscale','log');
hold on
plot(0:2:lmax1,ccfods(iii).pl1,'.r','markersize',30);
plot(0:2:lmax1,ccfods(iii).C1*ccfods(iii).lambda1.^(0:2:lmax1),'--k','linewidth',3)
xlabel('l'  ,'fontweight','bold','fontsize',30)
ylabel('p_l','fontweight','bold','fontsize',30)
title(['CC ' samples{f} ' t = ' num2str(ts(iii)) ' ms: C=' num2str(ccfods(iii).C1,'%.2f') ...
    ', lambda=' num2str(ccfods(iii).lambda1,'%.2f') ' disp angle= ' num2str(ccfods(iii).dispangp1,'%.2f') '/' num2str(ccfods(iii).dispang,'%.2f')], ...
    'fontsize',30,'interpreter','none')
ylim([0 1])
xlim([0 lmax1])
% pause(1)
% namefig = [rootDisp '/' samples{f} '_cc_FOD_pl_lmax' num2str(lmax1) '_t' num2str(ts(iii))];
% SaveFigurePNG([namefig '.png'],16*3,9*3);
% savefig(      [namefig '.fig'])
% pause(1)
close(k); k = k + 1;
% lmax2
figure(k); 
set(gca,'yscale','log');
hold on
plot(0:2:lmax2,ccfods(iii).pl2,'.r','markersize',30);
plot(0:2:lmax2,ccfods(iii).C2*ccfods(iii).lambda2.^(0:2:lmax2),'--k','linewidth',3)
xlabel('l'  ,'fontweight','bold','fontsize',30)
ylabel('p_l','fontweight','bold','fontsize',30)
title(['CC ' samples{f} ' t = ' num2str(ts(iii)) ' ms: C=' num2str(ccfods(iii).C2,'%.2f') ...
    ', lambda=' num2str(ccfods(iii).lambda2,'%.2f') ' disp angle= ' num2str(ccfods(iii).dispangp2,'%.2f') '/' num2str(ccfods(iii).dispang,'%.2f')],...
'fontsize',30,'interpreter','none')
ylim([0 1])
xlim([0 lmax2])
% pause(1)
% namefig = [rootDisp '/' samples{f} '_cc_FOD_pl_lmax' num2str(lmax2) '_t' num2str(ts(iii))];
% SaveFigurePNG([namefig '.png'],16*3,9*3);
% savefig(      [namefig '.fig'])
% pause(1)
close(k); k = k + 1;
% save
t         = ts(iii);
fod       = ccfods(iii);
fod.sinu  = ccttSmooth(iii).tort';
fod.t     = t;
fod.lmax1 = lmax1;
fod.lmax2 = lmax2;
save([rootDisp '/' samples{f} '_cc_fod_t' num2str(ts(iii)) '.mat'],'fod');
% save([rootDisp '/' samples{f} '_cc_fod_t' num2str(ts(iii)) '.mat'],'fod','t','lmax1','lmax2');
clear('fod','ccfods')
end % iii 
end % if
end % cc
% ----------------------------------------------------------------------- % 
if(flagcg)
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
% cgskSmooth
[~,cgtgSmooth,cgttSmooth] = SmoothedSkeletonsTangents(cgskeleton,cgaxprop,nmin,sigmas,cgaxidx);
clear('cgskeleton','cgaxprop')
% concatenate
cgtangent = ConcatenateTangents(cgtgSmooth);
% ----------------------------------------------------------------------- %
% % saving
% for iii = 1:length(cgskSmooth)
% t_D_vs_s = [ts(iii) D vs Sigmas(iii)];
% sksmooth = cgskSmooth(iii);
% tgsmooth = cgtgSmooth(iii);
% save([rootDisp '/' samples{f} '_cg_smooth_skeletons_t' num2str(ts(iii),'%.1f') '.mat'],'sksmooth','t_D_vs_s');
% save([rootDisp '/' samples{f} '_cg_smooth_tangents_t'  num2str(ts(iii),'%.1f') '.mat'],'tgsmooth','t_D_vs_s');
% end % iii    
% ----------------------------------------------------------------------- %
clear('cgskSmooth','cgtgSmooth')
if(~fcmb)
% fods    
cgfods = ComputeFODs(cgtangent,lmax1,lmax2);
clear('cgtangent')
% save images 
k     = 1;
lmaxs = [lmax1 lmax2];
for iii = 1:length(cgfods)
for lll = 1:2    
figure(k)
pause(1)
title(['CG ' samples{f} ' lmax = ' num2str(lmaxs(lll)) ' t = ' num2str(ts(iii)) ' ms'],...
    'fontsize',30,'interpreter','none')
pause(1)
namefig = [rootDisp '/' samples{f} '_cg_FOD_SH_lmax' num2str(lmaxs(lll)) '_t' num2str(ts(iii))];
SaveFigurePNG([namefig '.png'],16*3,9*3);
savefig(      [namefig '.fig'])
pause(1)
close(k); k = k + 1;
end % lll
end % iii
% plots
k = 1;
for iii = 1:length(cgfods)
% lmax1
figure(k); 
set(gca,'yscale','log');
hold on
plot(0:2:lmax1,cgfods(iii).pl1,'.r','markersize',30);
plot(0:2:lmax1,cgfods(iii).C1*cgfods(iii).lambda1.^(0:2:lmax1),'--k','linewidth',3)
xlabel('l'  ,'fontweight','bold','fontsize',30)
ylabel('p_l','fontweight','bold','fontsize',30)
title(['CG ' samples{f} ' t = ' num2str(ts(iii)) ' ms: C=' num2str(cgfods(iii).C1,'%.2f') ...
    ', lambda=' num2str(cgfods(iii).lambda1,'%.2f') ' disp angle= ' num2str(cgfods(iii).dispangp1,'%.2f') '/' num2str(cgfods(iii).dispang,'%.2f')],...
    'fontsize',30,'interpreter','none')
ylim([0 1])
xlim([0 lmax1])
% pause(1)
% namefig = [rootDisp '/' samples{f} '_cg_FOD_pl_lmax' num2str(lmax1) '_t' num2str(ts(iii))];
% SaveFigurePNG([namefig '.png'],16*3,9*3);
% savefig(      [namefig '.fig'])
% pause(1)
close(k); k = k + 1;
% lmax2
figure(k); 
set(gca,'yscale','log');
hold on
plot(0:2:lmax2,cgfods(iii).pl2,'.r','markersize',30);
plot(0:2:lmax2,cgfods(iii).C2*cgfods(iii).lambda2.^(0:2:lmax2),'--k','linewidth',3)
xlabel('l'  ,'fontweight','bold','fontsize',30)
ylabel('p_l','fontweight','bold','fontsize',30)
title(['CG ' samples{f} ' t = ' num2str(ts(iii)) ' ms: C=' num2str(cgfods(iii).C2,'%.2f') ...
    ', lambda=' num2str(cgfods(iii).lambda2,'%.2f') ' disp angle= ' num2str(cgfods(iii).dispangp2,'%.2f') '/' num2str(cgfods(iii).dispang,'%.2f')], ...
    'fontsize',30,'interpreter','none')
ylim([0 1])
xlim([0 lmax2])
% pause(1)
% namefig = [rootDisp '/' samples{f} '_cg_FOD_pl_lmax' num2str(lmax2) '_t' num2str(ts(iii))];
% SaveFigurePNG([namefig '.png'],16*3,9*3);
% savefig(      [namefig '.fig'])
% pause(1)
close(k); k = k + 1;
% save
t         = ts(iii);
fod       = cgfods(iii);
fod.sinu  = cgttSmooth(iii).tort';
fod.t     = t;
fod.lmax1 = lmax1;
fod.lmax2 = lmax2;
save([rootDisp '/' samples{f} '_cg_fod_t' num2str(ts(iii)) '.mat'],'fod');
% save([rootDisp '/' samples{f} '_cg_fod_t' num2str(ts(iii)) '.mat'],'fod','t','lmax1','lmax2');
clear('fod','cgfods')
end % iii 
end % if
end % cg

if( fcmb && flagcg && flagcc )
% concatenate
for iii = length(sigmas):-1:1
cccgtangent(iii).tg = [];
cccgtangent(iii).tg = cat(1,cccgtangent(iii).tg,cctangent(iii).tg);
cccgtangent(iii).tg = cat(1,cccgtangent(iii).tg,cgtangent(iii).tg);
end
% fods
cccgfods = ComputeFODs(cccgtangent,lmax1,lmax2);
clear('cccgtangent')
% save images 
k     = 1;
lmaxs = [lmax1 lmax2];
for iii = 1:length(cccgfods)
for lll = 1:2    
figure(k)
pause(1)
title(['CC-CG ' samples{f} ' lmax = ' num2str(lmaxs(lll)) ' t = ' num2str(ts(iii)) ' ms'],...
    'fontsize',30,'interpreter','none')
pause(1)
namefig = [rootDisp '/' samples{f} '_cccg_FOD_SH_lmax' num2str(lmaxs(lll)) '_t' num2str(ts(iii))];
SaveFigurePNG([namefig '.png'],16*3,9*3);
savefig(      [namefig '.fig'])
pause(1)
close(k); k = k + 1;
end % lll
end % iii
% plots
k = 1;
for iii = 1:length(cccgfods)
% lmax1
figure(k); 
set(gca,'yscale','log');
hold on
plot(0:2:lmax1,cccgfods(iii).pl1,'.r','markersize',30);
plot(0:2:lmax1,cccgfods(iii).C1*cccgfods(iii).lambda1.^(0:2:lmax1),'--k','linewidth',3)
xlabel('l'  ,'fontweight','bold','fontsize',30)
ylabel('p_l','fontweight','bold','fontsize',30)
title(['CC-CG ' samples{f} ' t = ' num2str(ts(iii)) ' ms: C=' num2str(cccgfods(iii).C1,'%.2f') ...
    ', lambda=' num2str(cccgfods(iii).lambda1,'%.2f') ' disp angle= ' num2str(cccgfods(iii).dispangp1,'%.2f') '/' num2str(cccgfods(iii).dispang,'%.2f')], ...
    'fontsize',30,'interpreter','none')
ylim([0 1])
xlim([0 lmax1])
% pause(1)
% namefig = [rootDisp '/' samples{f} '_cccg_FOD_pl_lmax' num2str(lmax1) '_t' num2str(ts(iii))];
% SaveFigurePNG([namefig '.png'],16*3,9*3);
% savefig(      [namefig '.fig'])
% pause(1)
close(k); k = k + 1;
% lmax2
figure(k); 
set(gca,'yscale','log');
hold on
plot(0:2:lmax2,cccgfods(iii).pl2,'.r','markersize',30);
plot(0:2:lmax2,cccgfods(iii).C2*cccgfods(iii).lambda2.^(0:2:lmax2),'--k','linewidth',3)
xlabel('l'  ,'fontweight','bold','fontsize',30)
ylabel('p_l','fontweight','bold','fontsize',30)
title(['CC-CG ' samples{f} ' t = ' num2str(ts(iii)) ' ms: C=' num2str(cccgfods(iii).C2,'%.2f') ...
    ', lambda=' num2str(cccgfods(iii).lambda2,'%.2f') ' disp angle= ' num2str(cccgfods(iii).dispangp2,'%.2f') '/' num2str(cccgfods(iii).dispang,'%.2f')], ...
    'fontsize',30,'interpreter','none')
ylim([0 1])
xlim([0 lmax2])
% pause(1)
% namefig = [rootDisp '/' samples{f} '_cccg_FOD_pl_lmax' num2str(lmax2) '_t' num2str(ts(iii))];
% SaveFigurePNG([namefig '.png'],16*3,9*3);
% savefig(      [namefig '.fig'])
% pause(1)
close(k); k = k + 1;
% save
t         = ts(iii);
fod       = cccgfods(iii);
fod.t     = t;
fod.lmax1 = lmax1;
fod.lmax2 = lmax2;
save([rootDisp '/' samples{f} '_cccg_fod_t' num2str(ts(iii)) '.mat'],'fod');
clear('fod','cccgfods')
end % iii 

end % combined

% ----------------------------------------------------------------------- % 
clear('skel','ccskel','ccskeleton','cgskel','cgskeleton');
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

function tangent = ConcatenateTangents(tgSmooth)
tic
for s = length(tgSmooth):-1:1
nax = length(tgSmooth(s).fiber);
tg = [];
% jsample = randsample(nax,1000)';
for j = 1:nax
% for j = jsample
tg = cat(1,tg,tgSmooth(s).fiber(j).tg);    
end
tangent(s).tg = tg;
end % s
toc
end % ConcatenateTangents

%
function fods = ComputeFODs(tangent,lmax1,lmax2)
as     = analyzeseg_modified();  
flag   = 1; 
denstr = {'low','medium','high'};
k = 1;
for s = length(tangent):-1:1
[fods(s).dispang,fods(s).maindir] = as.dispersionangle3d(tangent(s).tg);
tic
[fods(s).fodi,fods(s).points,fods(s).fodori,fods(s).solang] = ...
as.fodsphere(double([tangent(s).tg;-tangent(s).tg]),'density',denstr{flag},'range',[0 10]);
toc
tic
[fods(s).pl1    ,fods(s).plm1] = as.rotinv(fods(s).fodi,lmax1,flag);
[fods(s).pl2    ,fods(s).plm2] = as.rotinv(fods(s).fodi,lmax2,flag);
[fods(s).lambda1,fods(s).C1  ] = as.rotinvpoisson(fods(s).pl1(2:end),2:2:lmax1);
[fods(s).lambda2,fods(s).C2  ] = as.rotinvpoisson(fods(s).pl2(2:end),2:2:lmax2);
fods(s).dispangp1              = as.dispersionanglep2(fods(s).pl1(2));
fods(s).dispangp2              = as.dispersionanglep2(fods(s).pl2(2));
toc
tic
figure(k); k = k + 1; 
as.fodsht(fods(s).fodi,lmax1,'smooth','false','density',denstr{flag},'range',[0 10],'colorbar','off','glyph','on');
figure(k); k = k + 1; 
as.fodsht(fods(s).fodi,lmax2,'smooth','false','density',denstr{flag},'range',[0 10],'colorbar','off','glyph','on');
toc
end % s
end % AnalysisDispersion

