function script07_compute_axon_morph_props() 
% This script computes axons morphology properties 
%
% This script needs the results from scripts:
% script06_axon_diam_skel_props.m
% It also needs the functions in: 
% RefineSeg.m
% and the information in the .mat files in the folder
% cc_cg_labels/
% The results of this script will be saved in folder
% White_matter_EM_02_CrossSections_Props/*/
% 24-April-2023
close all
clc
rng(0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Begin\n')
ir          = 2; % 1 = HM, 2 = LM
res         = {'HM','LM'};
rootIn      = ['White_matter_EM_Diam_Skel_Props_' res{ir} '/']; 
fOut        = 'White_matter_EM_02_CrossSections_Props/';
folderqn    = 'cc_cg_labels/';
if(~isfolder(fOut)); mkdir(fOut); end
% name of files needed
areasaxons  = 'axons_areas_cross_sections.mat';
skelaxons   = 'axons_skeletons.mat';
propsaxons  = 'axons_props.mat';
% sufixqn     = '000_animal_exp_axons_cells_nx_ny_nz_cclims_cglims.txt';
% samples
samples     = {'Sham_25_contra','Sham_25_ipsi','Sham_49_contra','Sham_49_ipsi',...
               'TBI_24_contra','TBI_24_ipsi','TBI_28_contra','TBI_28_ipsi',...
               'TBI_2_contra','TBI_2_ipsi'};
idsample    = [25 25 49 49 24 24 28 28 2 2 ]; 
sidesample  = {'contra','ipsi','contra','ipsi','contra','ipsi',...
               'contra','ipsi','contra','ipsi'};
nsamples    = length(samples);
voxsizes    = [15 50]; % now volumes are isotropic
vs          = voxsizes(ir)*1e-3;
Dw          = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ir == 2)
lmin     = 1.0;              % 1 um at extremes is removed
lrem     = 5.0; % 1.0; %     % min length of axons = 2*lrem 
else
lmin     = 0.5;              % 0.5 um at extremes is removed
lrem     = 2.5; % 1.0; %     % min length of axons = 2*lrem 
end
maxPmany = 0.25;             % proporion of cross sections with multiple skeletons
nmin     = round(lmin/vs);   % " in voxels 
nrem     = round(lrem/vs)-1; % " in voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f = nsamples:-1:1
sample = [num2str(idsample(f)) '_' sidesample{f}];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data
tic
csname = [rootIn samples{f} '/' res{ir} '_' sample '_' areasaxons];
csasea = load(csname).areas_cross_sections;  
skname = [rootIn samples{f} '/' res{ir} '_' sample '_' skelaxons];
skelet = load(skname).skeletons;  
prname = [rootIn samples{f} '/' res{ir} '_' sample '_' propsaxons];
axprop = load(prname).axon_props;  
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute dispersion angles caused by undulation (aligned axons)
ts       = 0;              % diffusion times (ms)
D        = 2;              % intrinsic diffusivity (micron^2/ms)
Sigmas   = sqrt(2*D*ts)/2; % smoothing kernel width (micron), sigma = L/2 -> L=sqrt(2*D*t)
sigmas   = Sigmas/vs;      % smoothing kernel width (pixels) 
axidx    = [f 0];
[~,~,~,adSmooth] = SmoothedSkeletonsTangents(skelet,axprop,nmin,sigmas,axidx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read some props and set flags
if(ir == 2)
lblsname = [folderqn '/LM_' sample '_cc_cg_lbls.mat'];
lblscccg = load(lblsname);
fcc      = ~isempty(lblscccg.cc_lbl);
fcg      = ~isempty(lblscccg.cg_lbl);
% is       = f;
% props    = load([folderqn sufixqn]); % samples properties and cc,cg limits
% propsi   = props(:,is);
else
fcc    = 1;
fcg    = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- % 
% Number of axons
naxsel = [36561 37335 27529 13840 39477 26320 33612 31103 31732 23914];
nax    = size(csasea,1);
% ----------------------------------------------------------------------- % 
% get axon nz length (tricky: some axons are empty, this modifies tab datatype from array to cell)
% anewsc09_220905_dispersion_em_da_aligned() uses a different method but it
% should provide equal results than this
Nzo    = axprop(:,'euclideanlen').euclideanlen; 
Pmany  = axprop(:,'nummanyskel').nummanyskel; 
if(~isnumeric(Nzo) && ~isnumeric(Pmany))
for j = 1:nax    
if(isempty(Nzo{j}))
Nzo{j} = 0;
end
if(isempty(Pmany{j}))
Pmany{j} = maxPmany;
end
end % j
Nzo   = cell2mat(Nzo);        % euclidean length 
Pmany = cell2mat(Pmany)./Nzo; % proportion of cross sections with multiple skeletons
end % 
% ----------------------------------------------------------------------- % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcc && fcg)
% selecting axons only in cc/cg region   
if(ir == 2)
cccglbls      = 1:naxsel(f);
else
cccglbls      = [];    
end
cccgaxprop1   = axprop(cccglbls,:); 
cccgcsasea1   = csasea(cccglbls,:);
cccgskelet1   = skelet(cccglbls,:);
cccgNzo       = Nzo(cccglbls);
cccgPmany     = Pmany(cccglbls);
cccgdang1     = adSmooth(1).adisp(cccglbls); 
clear('cccglbls')
% remove very small axons and axons with many multiple skeletons per section 
cccgidx       = ( cccgNzo > (2*nrem+2) ) & ( cccgPmany < maxPmany ); 
clear('cccgNzo','cccgPmany');
cccgaxprop    = cccgaxprop1(cccgidx,:); clear('cccgaxprop1');
cccgskelet    = cccgskelet1(cccgidx,:); clear('cccgskelet1');
cccgcsasea    = cccgcsasea1(cccgidx,:); clear('cccgcsasea1'); 
cccgdang      = cccgdang1(cccgidx);     clear('cccgdang1'); 
cccgdisp      = cosd(cccgdang').^2;
clear('cccgidx');
% remove extremes and compute parameters: diam, cvr, A/A(z)
cccgnax       = size(cccgaxprop,1);
axvola        = zeros(cccgnax,1);
euclength     = zeros(cccgnax,1);
geolength     = zeros(cccgnax,1);
sinuosity     = zeros(cccgnax,1);
diamsalignm   = zeros(cccgnax,1);
cv2ralign     = zeros(cccgnax,1);
torteucmalign = zeros(cccgnax,1);
diamsalign    = [];
% ----------------------------------------------------------------------- % 
for i = 1:cccgnax
% 01 lengths    
skeleton      = cccgskelet(i,1).skeleton{:}(nmin:end-nmin-1,:);
[a,b,c]       = RefineSeg.FiberLengthsNew(skeleton);
euclength(i)  = a; % eulidean
geolength(i)  = b; % geodesic
sinuosity(i)  = c; % sinuosity
clear('skeleton');
% 02 tortuosity (aligned)
csareasalign     = cccgcsasea(i,2).csareasalign{:}(nmin:end-nmin-1);  
csareasalign(csareasalign == 0) = [];
csareasvola      = sum(csareasalign);          % volume (sum of pixels on cross sections)
axvola(i)        = csareasvola;
csareasmeuca     = csareasvola/euclength(i);   % mean(A(z)) = volume / length
torteucm2        = csareasmeuca./csareasalign; % mean(A(z))/A(z)
torteucmalign(i) = mean(torteucm2);            % mean( mean(A)/A(z) ) 
% 03 diameter and cvr (aligned)
diams2         = 2*sqrt(csareasalign/pi);
diamsalignm(i) = mean(diams2);
diamsaligns    = std(diams2);
diamsalign     = cat(1,diamsalign,diams2);
cv2ralign(i)   = diamsaligns/diamsalignm(i);
clear('csareasalign','torteucm2','diams2');
end % i
% ----------------------------------------------------------------------- % 
clear('cccgaxprop','cccgskelet','cccgcsasea');
% Average lengths
cccgstats(f).voxelsize             = vs;
cccgstats(f).euclideanlengthmean   = mean(euclength); 
cccgstats(f).euclideanlengthmedian = median(euclength);
cccgstats(f).euclideanlengthstd    = std(euclength);
cccgstats(f).geodesiclengthmean    = mean(geolength); 
cccgstats(f).geodesiclengthmedian  = median(geolength);
cccgstats(f).geodesiclengthstd     = std(geolength);
cccgstats(f).sinuositymean         = mean(sinuosity); 
cccgstats(f).sinuositymedian       = median(sinuosity);
cccgstats(f).sinuositystd          = std(sinuosity);
% Predicted Da 
cccgstats(f).Daxonw_aligned        = sum(axvola.*(Dw*cccgdisp./torteucmalign))./sum(axvola); % <<---
% Average Tortuosity (Aligned)
cccgstats(f).tortuosityaligneucmean   = mean(torteucmalign);
cccgstats(f).tortuosityaligneucmedian = median(torteucmalign);
cccgstats(f).tortuosityaligneucstd    = std(torteucmalign);
% Average Diameters (Aligned)
cccgstats(f).diameteralignmeanW       = sum(axvola.*diamsalignm)./sum(axvola); % <<--- 
cccgstats(f).diameteralignmean        = mean(diamsalignm); 
cccgstats(f).diameteralignmedian      = median(diamsalignm); 
cccgstats(f).diameteralignstd         = std(diamsalignm); 
cccgstats(f).cv2ralignmean            = mean(cv2ralign); 
cccgstats(f).cv2ralignmedian          = median(cv2ralign); 
cccgstats(f).cv2ralignstd             = std(cv2ralign); 
clear('torteucmalign');
% Big Lengths
cccglengths(f).euclideanlength = euclength; 
cccglengths(f).geodesiclength  = geolength;
cccglengths(f).sinuosity       = sinuosity;
clear('euclength','geolength','sinuosity');
% Big Diameters
cccgdiameters(f).diameteraligned      = diamsalign;
cccgdiameters(f).diameteralignedmeans = diamsalignm;
cccgdiameters(f).cv2ralign            = cv2ralign;
clear('diamsalign','diamsalignm','cv2ralign');
end % cccg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcc) 
% selecting axons only in cc region   
if(ir == 2)
cclbls        = lblscccg.cc_lbl(lblscccg.cc_lbl<=naxsel(f));
else
cclbls        = 1:nax;    
end
ccaxprop1     = axprop(cclbls,:); 
cccsasea1     = csasea(cclbls,:);
ccskelet1     = skelet(cclbls,:);
ccNzo         = Nzo(cclbls);
ccPmany       = Pmany(cclbls);
ccdang1       = adSmooth(1).adisp(cclbls); 
clear('cclbls')
% remove very small axons and axons with many multiple skeletons per section 
ccidx         = ( ccNzo > (2*nrem+2) ) & ( ccPmany < maxPmany ); 
clear('ccNzo','ccPmany');
ccaxprop      = ccaxprop1(ccidx,:); clear('ccaxprop1');
ccskelet      = ccskelet1(ccidx,:); clear('ccskelet1');
cccsasea      = cccsasea1(ccidx,:); clear('cccsasea1'); 
ccdang        = ccdang1(ccidx);     clear('ccdang1'); 
ccdisp        = cosd(ccdang').^2;
clear('ccidx');
% remove extremes and compute parameters: diam, cvr, A/A(z)
ccnax         = size(ccaxprop,1);
axvola        = zeros(ccnax,1);
euclength     = zeros(ccnax,1);
geolength     = zeros(ccnax,1);
sinuosity     = zeros(ccnax,1);
diamsalignm   = zeros(ccnax,1);
cv2ralign     = zeros(ccnax,1);
torteucmalign = zeros(ccnax,1);
diamsalign    = [];
% ----------------------------------------------------------------------- % 
for i = 1:ccnax
% 01 lengths    
skeleton      = ccskelet(i,1).skeleton{:}(nmin:end-nmin-1,:);
[a,b,c]       = RefineSeg.FiberLengthsNew(skeleton);
euclength(i)  = a; % eulidean
geolength(i)  = b; % geodesic
sinuosity(i)  = c; % sinuosity
clear('skeleton');
% 02 tortuosity (aligned)
csareasalign     = cccsasea(i,2).csareasalign{:}(nmin:end-nmin-1);  
csareasalign(csareasalign == 0) = [];
csareasvola      = sum(csareasalign);          % volume (sum of pixels on cross sections)
axvola(i)        = csareasvola;
csareasmeuca     = csareasvola/euclength(i);   % mean(A(z)) = volume / length
torteucm2        = csareasmeuca./csareasalign; % mean(A(z))/A(z) 
torteucmalign(i) = mean(torteucm2);            % mean( mean(A)/A(z) )  
% 03 diameter and cvr (aligned)
diams2         = 2*sqrt(csareasalign/pi);
diamsalignm(i) = mean(diams2);
diamsaligns    = std(diams2);
diamsalign     = cat(1,diamsalign,diams2);
cv2ralign(i)   = diamsaligns/diamsalignm(i);
clear('csareasalign','torteucm2','diams2');
end % i
clear('ccaxprop','ccskelet','cccsasea');
% ----------------------------------------------------------------------- % 
% Average lengths
ccstats(f).voxelsize             = vs;
ccstats(f).euclideanlengthmean   = mean(euclength); 
ccstats(f).euclideanlengthmedian = median(euclength);
ccstats(f).euclideanlengthstd    = std(euclength);
ccstats(f).geodesiclengthmean    = mean(geolength); 
ccstats(f).geodesiclengthmedian  = median(geolength);
ccstats(f).geodesiclengthstd     = std(geolength);
ccstats(f).sinuositymean         = mean(sinuosity); 
ccstats(f).sinuositymedian       = median(sinuosity);
ccstats(f).sinuositystd          = std(sinuosity);
% Predicted Da 
ccstats(f).Daxonw_aligned        = sum(axvola.*(Dw*ccdisp./torteucmalign))./sum(axvola);   % <<---
% Average Tortuosity (Aligned)
ccstats(f).tortuosityaligneucmean   = mean(torteucmalign);
ccstats(f).tortuosityaligneucmedian = median(torteucmalign);
ccstats(f).tortuosityaligneucstd    = std(torteucmalign);
% Average Diameters (Aligned)
ccstats(f).diameteralignmeanW       = sum(axvola.*diamsalignm)./sum(axvola); % <<---
ccstats(f).diameteralignmean        = mean(diamsalignm); 
ccstats(f).diameteralignmedian      = median(diamsalignm); 
ccstats(f).diameteralignstd         = std(diamsalignm); 
ccstats(f).cv2ralignmean            = mean(cv2ralign); 
ccstats(f).cv2ralignmedian          = median(cv2ralign); 
ccstats(f).cv2ralignstd             = std(cv2ralign); 
clear('torteucmalign');
% Big Lengths
cclengths(f).euclideanlength = euclength; 
cclengths(f).geodesiclength  = geolength;
cclengths(f).sinuosity       = sinuosity;
clear('euclength','geolength','sinuosity');
% Big Diameters
ccdiameters(f).diameteraligned      = diamsalign;
ccdiameters(f).diameteralignedmeans = diamsalignm;
ccdiameters(f).cv2ralign            = cv2ralign;
clear('diamsalign','diamsalignm','cv2ralign');
end % cc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcg)
% selecting axons only in cg region   
if(ir == 2)
cglbls        = lblscccg.cg_lbl(lblscccg.cg_lbl<=naxsel(f));
else
cglbls        = [];    
end
cgaxprop1     = axprop(cglbls,:); 
cgcsasea1     = csasea(cglbls,:);
cgskelet1     = skelet(cglbls,:);
cgNzo         = Nzo(cglbls);
cgPmany       = Pmany(cglbls);
cgdang1       = adSmooth(1).adisp(cglbls); 
clear('cglbls')
% remove very small axons and axons with many multiple skeletons per section 
cgidx         = ( cgNzo > (2*nrem+2) ) & ( cgPmany < maxPmany ); 
clear('cgNzo','cgPmany');
cgaxprop      = cgaxprop1(cgidx,:); clear('cgaxprop1');
cgskelet      = cgskelet1(cgidx,:); clear('cgskelet1');
cgcsasea      = cgcsasea1(cgidx,:); clear('cgcsasea1'); 
cgdang        = cgdang1(cgidx);     clear('cgdang1'); 
cgdisp        = cosd(cgdang').^2;
clear('cgidx');
% remove extremes and compute parameters: diam, cvr, A/A(z)
cgnax         = size(cgaxprop,1);
axvola        = zeros(cgnax,1);
euclength     = zeros(cgnax,1);
geolength     = zeros(cgnax,1);
sinuosity     = zeros(cgnax,1);
diamsalignm   = zeros(cgnax,1);
cv2ralign     = zeros(cgnax,1);
torteucmalign = zeros(cgnax,1);
diamsalign    = [];
% ----------------------------------------------------------------------- % 
for i = 1:cgnax
% 01 lengths    
skeleton      = cgskelet(i,1).skeleton{:}(nmin:end-nmin-1,:);
[a,b,c]       = RefineSeg.FiberLengthsNew(skeleton);
euclength(i)  = a; % eulidean
geolength(i)  = b; % geodesic
sinuosity(i)  = c; % sinuosity
clear('skeleton');
% 02 tortuosity (aligned)
csareasalign     = cgcsasea(i,2).csareasalign{:}(nmin:end-nmin-1);  
csareasalign(csareasalign == 0) = [];
csareasvola      = sum(csareasalign);           % volume (sum of pixels on cross sections)
axvola(i)        = csareasvola;
csareasmeuca     = csareasvola/euclength(i);    % mean(A(z)) = volume / length
torteucm2        = csareasmeuca./csareasalign;  % mean(A(z))/A(z) 
torteucmalign(i) = mean(torteucm2);             % mean( mean(A)/A(z) )  
% 03 diameter and cvr (aligned)
diams2         = 2*sqrt(csareasalign/pi);
diamsalignm(i) = mean(diams2);
diamsaligns    = std(diams2);
diamsalign     = cat(1,diamsalign,diams2);
cv2ralign(i)   = diamsaligns/diamsalignm(i);
clear('csareasalign','torteucm2','diams2');
end % i
% ----------------------------------------------------------------------- % 
clear('cgaxprop','cgskelet','cgcsasea');
% Average lengths
cgstats(f).voxelsize             = vs;
cgstats(f).euclideanlengthmean   = mean(euclength); 
cgstats(f).euclideanlengthmedian = median(euclength);
cgstats(f).euclideanlengthstd    = std(euclength);
cgstats(f).geodesiclengthmean    = mean(geolength); 
cgstats(f).geodesiclengthmedian  = median(geolength);
cgstats(f).geodesiclengthstd     = std(geolength);
cgstats(f).sinuositymean         = mean(sinuosity); 
cgstats(f).sinuositymedian       = median(sinuosity);
cgstats(f).sinuositystd          = std(sinuosity);
% Predicted Da 
cgstats(f).Daxonw_aligned        = sum(axvola.*(Dw*cgdisp./torteucmalign))./sum(axvola);   % <<---
% Average Tortuosity (Aligned)
cgstats(f).tortuosityaligneucmean   = mean(torteucmalign);
cgstats(f).tortuosityaligneucmedian = median(torteucmalign);
cgstats(f).tortuosityaligneucstd    = std(torteucmalign);
% Average Diameters (Aligned)
cgstats(f).diameteralignmeanW       = sum(axvola.*diamsalignm)./sum(axvola); % <<---
cgstats(f).diameteralignmean        = mean(diamsalignm); 
cgstats(f).diameteralignmedian      = median(diamsalignm); 
cgstats(f).diameteralignstd         = std(diamsalignm); 
cgstats(f).cv2ralignmean            = mean(cv2ralign); 
cgstats(f).cv2ralignmedian          = median(cv2ralign); 
cgstats(f).cv2ralignstd             = std(cv2ralign); 
clear('torteucmalign');
% Big Lengths
cglengths(f).euclideanlength = euclength; 
cglengths(f).geodesiclength  = geolength;
cglengths(f).sinuosity       = sinuosity;
clear('euclength','geolength','sinuosity');
% Big Diameters
cgdiameters(f).diameteraligned      = diamsalign;
cgdiameters(f).diameteralignedmeans = diamsalignm;
cgdiameters(f).cv2ralign            = cv2ralign;
clear('diamsalign','diamsalignm','cv2ralign');
end % cg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',sample)
end % f
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
cc_stats      = struct2table(ccstats);
cc_lengths    = struct2table(cclengths);
cc_diameter   = struct2table(ccdiameters);
clear('ccstats','cclengths','ccdiameters')
%
ccnames = [fOut '/cc_stats_' res{ir} '.mat'];
save(ccnames,'cc_stats'); 
ccnamel = [fOut '/cc_lengths_' res{ir} '.mat'];
save(ccnamel,'cc_lengths');  
ccnamed = [fOut '/cc_diameter_' res{ir} '.mat'];
save(ccnamed,'cc_diameter','-v7.3'); 
clear('cc_stats','cc_lengths','cc_diameter')
%
if(ir == 2)
cg_stats      = struct2table(cgstats);
cg_lengths    = struct2table(cglengths);
cg_diameter   = struct2table(cgdiameters);
%
cgnames = [fOut '/cg_stats_' res{ir} '.mat'];
save(cgnames,'cg_stats'); 
cgnamel = [fOut '/cg_lengths_' res{ir} '.mat'];
save(cgnamel,'cg_lengths'); 
cgnamed = [fOut '/cg_diameter_' res{ir} '.mat'];
save(cgnamed,'cg_diameter','-v7.3'); 
clear('cg_stats','cg_lengths','cg_diameter')

cccg_stats      = struct2table(cccgstats);
cccg_lengths    = struct2table(cccglengths);
cccg_diameter   = struct2table(cccgdiameters);
%
cccgnames = [fOut '/cccg_stats_' res{ir} '.mat'];
save(cccgnames,'cccg_stats'); 
cccgnamel = [fOut '/cccg_lengths_' res{ir} '.mat'];
save(cccgnamel,'cccg_lengths');  
cccgnamed = [fOut '/cccg_diameter_' res{ir} '.mat'];
save(cccgnamed,'cccg_diameter','-v7.3'); 
clear('cccg_stats','cccg_lengths','cccg_diameter')
end
toc
fprintf('End\n')
end % main

function [skSmooth,tgSmooth,ttSmooth,dispSmooth] = SmoothedSkeletonsTangents(skeleton,axprop,nmin,sigmas,axidx)
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
if(size(sk,1)>(2*nmin)) %%%%%%%%%%%%%%%% (NOT IN DISP CODE) !!!!!!!!!!!!!!!!
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
% % Restore smoothed skeleton (remove rotation that aligned it to z and translation)
% if(isnumeric(axprop(j,43).OriginNew))
% skr = ( axprop(j,45).MatRestore{:} * (cms' + axprop(j,43).OriginNew')  ); %     
% else
% skr = ( axprop(j,45).MatRestore{:} * (cms' + axprop(j,43).OriginNew{:}')  ); % 
% end
skr = cms'; % aligned skeleton (only dispersion due to undulations)
clear('cms') %%%%%%%%%%%%%%%% (DIFFERENT THAN DISP CODE) !!!!!!!!!!!!!!!!
% remove tangents with zero norm
% tgdiff = diff(skr,1,2)';%
tgdiff = gradient(skr)';
tgnrm  = sqrt(sum(tgdiff.*tgdiff,2));
tgnrmm = mean(tgnrm);  
tgidx  = abs(tgnrm) > thrnrm*tgnrmm; 
tgs    = normr(tgdiff(tgidx,:)); 
% Compute dispersion angle  %%%%%%%%%%%%%%%% (NOT IN DISP CODE) !!!!!!!!!!!!!!!!
[dispang,~] = dispersionangle3d(tgs);
% skrt = skr'; 
% figure; scatter3(skrt(:,1),skrt(:,2),skrt(:,3),10,'k','filled'); axis('equal'); xlabel('x'); ylabel('y'); zlabel('z');
% figure; hold on; scatter3(tgs(:,1),tgs(:,2),tgs(:,3),10,'k','filled'); axis('equal'); xlabel('x'); ylabel('y'); zlabel('z');
% set to outputs
dispSmooth(s).adisp(j)  = dispang;
ttSmooth(s).tort(j)     = t;
tgSmooth(s).fiber(j).tg = tgs;
skSmooth(s).fiber(j).sk = skr';
clear('tgdiff','tgnrm','tgnrmm','tgidx','tgs','skr');
end % s
else %%%%%%%%%%%%%%%% (NOT IN DISP CODE) !!!!!!!!!!!!!!!! 
for s = 1:length(sigmas)    
dispSmooth(s).adisp(j)  = 90;
ttSmooth(s).tort(j)     = 0;
tgSmooth(s).fiber(j).tg = [0 0 0];
skSmooth(s).fiber(j).sk = [0 0 0]; 
end % s
end %%%%%%%%%%%%%%%% (NOT IN DISP CODE) !!!!!!!!!!!!!!!!
if(~isempty(axidx) && mod(j-nax,round(nax/10))==0) %%%%%%%%%%%%%%%% (NOT IN DISP CODE) !!!!!!!!!!!!!!!!
fprintf('%d %d %d\n',axidx(1),axidx(2),j);
end
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

function [dispang,maindir] = dispersionangle3d(tangent)
%DISPERSIONANGLE3D    Dispersion angle based on rms of cosine.
%   DISPERSIONANGLE3D(tangent) returns dispersion angle calculated by the
%   root-mean-square of the cosine of each unit vector tangent wrt the main
%   direction.
%
% author: Hong-Hsi Lee, 2018

% Flip tangent to upper-half sphere
tangent = tangent.*repmat(sign(tangent(:,3)),[1,3]);

% Mean direction
[U,~,~] = svd(tangent.'*tangent/size(tangent,1));

% === Align tangent to the main direction === %
% Lie Algebra: Lx, Ly, Lz are the most common basis for so(3) 
Lx = [0 0 0; 0 0 -1; 0 1 0]; Ly = [0 0 1; 0 0 0; -1 0 0]; Lz = [0 -1 0; 1 0 0; 0 0 0];
uL = U(2,1)*Lx - U(1,1)*Ly; uL = uL/sqrt(U(1,1)^2+U(2,1)^2);
thetar = acos(U(3,1));
tangentr = expm(thetar*uL)*tangent.'; tangentr = tangentr.';

% Calculate dispersion angle of all tangents
w = tangentr(:,3);
dispang = acosd(rms(w));

if(nargout > 1)
    maindir = U(:,1);
end

end
