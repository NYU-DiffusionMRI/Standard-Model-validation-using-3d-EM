function script08_3D_aligned_axons()
% Plot a few representative 3D axons for each sample
%
% This script needs the results from scripts:
% script06_axon_diam_skel_props.m
% script07_compute_axon_morph_props.m
% It also needs the functions: 
% SaveFigurePNG.m
% SaveFigureSVG.m
% and the information in the .mat files in the folder
% cc_cg_labels/
% The results of this script will be saved in folder
% White_matter_EM_05_Plot_Aligned_Axons/*/
% 13-Jun-2023 by Ricardo Coronado-Leija
close all
clc
rng(0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Begin\n')
ir          = 2; % 1 = HM, 2 = LM (only for 2)
res         = {'HM','LM'};
rootIn      = ['White_matter_EM_Diam_Skel_Props_' res{ir} '/']; 
rootStats   = 'White_matter_EM_02_CrossSections_Props/';
fOut        = 'White_matter_EM_05_Plot_Aligned_Axons/';
folderqn    = 'cc_cg_labels/';
if(~isfolder(fOut)); mkdir(fOut); end
% name of files needed
straxons    = 'axons_skel_areacs_props.mat';
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
Naxons      = 12; 
if(ir == 2)
qfact       = 0.25; 
sfact       = 0.125; % (0.05)
else 
qfact       = 0.250;
sfact       = 0.125;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ir == 2) % more conservative values
% lmin     = 1.0;            % 1 um at extremes is removed
lrem     = 5.0;              % min length of axons = 2*lrem 
% lplot    = 35.0; % (20) 
lplot    = 50.0; % (20) 
else
% lmin     = 0.5;              % 0.5 um at extremes is removed
lrem     = 2.5;              % min length of axons = 2*lrem 
lplot    = 10.0; 
end
maxPmany = 0.25;             % proporion of cross sections with multiple skeletons
minPprop = 0.95;             % proportion of longest object in axon
% nmin     = round(lmin/vs);   % " in voxels 
nrem     = round(lrem/vs)-1; % " in voxels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data
tic
ccsname  = [rootStats 'cc_stats_' res{ir} '.mat']; 
ccstats  = load(ccsname).cc_stats;
ccdname  = [rootStats 'cc_diameter_' res{ir} '.mat']; 
ccdiams  = load(ccdname).cc_diameter;
cclname  = [rootStats 'cc_lengths_' res{ir} '.mat']; 
cclength = load(cclname).cc_lengths;
cccv2rm  = ccstats(:,'cv2ralignmean').Variables;
cccv2rs  = ccstats(:,'cv2ralignstd').Variables;
ccsinum  = ccstats(:,'sinuositymean').Variables;
ccdiamm  = ccstats(:,'diameteralignmeanm').Variables;
ccdawm   = ccstats(:,'Daxonw_align_euc').Variables;
if(~isnumeric(cccv2rm)); cccv2rm = cell2mat(cccv2rm); end
if(~isnumeric(cccv2rs)); cccv2rs = cell2mat(cccv2rs); end
if(~isnumeric(ccsinum)); ccsinum = cell2mat(ccsinum); end
if(~isnumeric(ccdiamm)); ccdiamm = cell2mat(ccdiamm); end
if(~isnumeric(ccdawm));  ccdawm  = cell2mat(ccdawm);  end
if(ir == 2)
cccv2rm  = [cccv2rm(1:3); 0; cccv2rm(4:end)];
cccv2rs  = [cccv2rs(1:3); 0; cccv2rs(4:end)];
ccsinum  = [ccsinum(1:3); 0; ccsinum(4:end)];
ccdiamm  = [ccdiamm(1:3); 0; ccdiamm(4:end)]*0.05;
ccdawm   = [ccdawm(1:3);  0; ccdawm(4:end)];
cgsname  = [rootStats 'cg_stats_' res{ir} '.mat']; 
cgstats  = load(cgsname).cg_stats;
cgdname  = [rootStats 'cg_diameter_' res{ir} '.mat']; 
cgdiams  = load(cgdname).cg_diameter;
cglname  = [rootStats 'cg_lengths_' res{ir} '.mat']; 
cglength = load(cglname).cg_lengths;
cgcv2rm  = cgstats(:,'cv2ralignmean').Variables;
cgcv2rs  = cgstats(:,'cv2ralignstd').Variables;
cgsinum  = cgstats(:,'sinuositymean').Variables;
cgdiamm  = cgstats(:,'diameteralignmeanm').Variables*0.05;
cgdawm   = cgstats(:,'Daxonw_align_euc').Variables;
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first run everything without plotting/saving and check to set this values
fPlot = rand > 0;
% fcvr  = rand > 0; % match mean for cvr only (for Els Figure)
fcvr = rand < 0; % match means for cvr, sinu & diam 
% lplot=50, sfact=0.125, qfact=0.25 - match cvr & sinu (12 axons)
iwsize=240; ihsize=100; fslope = false;
ccszmax    = [113 237 2300];
cgszmax    = [153 207 2220];
cclplotmax = 2000*vs; %  100
cglplotmax = 1000*vs; %   50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f = nsamples:-1:1
sample = [num2str(idsample(f)) '_' sidesample{f}];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
axname   = [rootIn samples{f} '/' res{ir} '_' sample '_' straxons];
axprop   = load(axname).axon_skel_areacs_props;  
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read some props and set flags
if(ir == 2)
lblsname = [folderqn '/LM_' sample '_cc_cg_lbls.mat'];
lblscccg = load(lblsname);
fcc      = ~isempty(lblscccg.cc_lbl);
fcg      = ~isempty(lblscccg.cg_lbl);
else
fcc    = 1;
fcg    = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- % 
% Number of axons
naxsel = [36561 37335 27529 13840 39477 26320 33612 31103 31732 23914];
nax    = size(axprop,1);
% get axon nz length (tricky: some axons are empty, this modifies tab datatype from array to cell)
Nzo    = axprop(:,'euclideanlen').euclideanlen; 
Pmany  = axprop(:,'nummanyskel').nummanyskel; 
Pprop  = axprop(:,'proplargest').proplargest;
if(~isnumeric(Nzo) && ~isnumeric(Pmany) && ~isnumeric(Pprop))
for j = 1:nax
if(isempty(Nzo{j}))
Nzo{j} = 0;
end
if(isempty(Pmany{j}))
Pmany{j} = maxPmany;
end
if(isempty(Pprop{j}))
Pprop{j} = minPprop;
end
end % j
Nzo   = cell2mat(Nzo);        % euclidean length 
Pmany = cell2mat(Pmany)./Nzo; % proportion of cross sections with multiple skeletons
Pprop = cell2mat(Pprop);      % proportion largest object 
end
% ----------------------------------------------------------------------- % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcc) 
% selecting axons only in cc region   
if(ir == 2)
cclbls     = lblscccg.cc_lbl(lblscccg.cc_lbl<=naxsel(f));
else
cclbls     = 1:nax;    
end
ccaxprop1  = axprop(cclbls,1:2); 
ccNzo      = Nzo(cclbls);
ccPmany    = Pmany(cclbls);
ccPprop    = Pprop(cclbls);
clear('cclbls')
% 1st selection: remove very small axons and axons with many multiple skeletons per section
ccidx      = ( ccNzo > (2*nrem+2) ) & ( ccPmany < maxPmany ); % & (ccPprop > minPprop); 
clear('ccNzo','ccPmany');
ccaxprop   = ccaxprop1(ccidx,:);
ccPprop    = ccPprop(ccidx);
clear('ccaxprop1','ccidx');
% 2nd selection: select base in their length, diameter, cv2r and sinuosity
ccdiamsm   = cell2mat(ccdiams(f,'diameteralignedmeans').Variables)*vs;
cccv2r     = cell2mat(ccdiams(f,'cv2ralign').Variables);
ccsinu     = cell2mat(cclength(f,'sinuosity').Variables);
cclengths  = cell2mat(cclength(f,'euclideanlength').Variables)*vs;
ccidx2     = (cclengths > lplot) & (ccsinu < quantile(ccsinu,1-qfact)) & ...
             (cccv2rm(f)-sfact*cccv2rs(f) < cccv2r & cccv2r < cccv2rm(f)+sfact*cccv2rs(f));
% 3rd selection: closests to means
ccidx2     = ccidx2 & ccPprop > minPprop;
ccind      = find(ccidx2);
cclen      = cclengths(ccidx2);
if(fcvr)
cccvdif    = abs(cccv2rm(f) - cccv2r(ccidx2)); % new
else
cccvdif    = abs(cccv2rm(f)-cccv2r(ccidx2))+2*abs(ccsinum(f)-ccsinu(ccidx2))+abs(ccdiamm(f)-ccdiamsm(ccidx2)); % new
end
fprintf('cc: %d - (%f) %f %f %f - %f %f, %f %f %f\n',length(ccind),ccdawm(f),cccv2rm(f),ccsinum(f),ccdiamm(f),...
min(cclen),max(cclen),mean(ccdiamsm(ccind)),mean(cccv2r(ccind)),mean(ccsinu(ccind)))
clear('cclen');
% [~,ccsidx] = sort(cclen,'descend');
[~,ccsidx] = sort(cccvdif,'ascend'); % new
Nnaxons    = min( [Naxons length(ccsidx)] );
ccsample   = ccind(ccsidx(1:Nnaxons));
ccaxsel    = ccaxprop(ccsample,1); 
if(~isnumeric(ccaxprop(ccsample,2).size))
ccszsel    = cell2mat(ccaxprop(ccsample,2).size); 
else
ccszsel    = ccaxprop(ccsample,2).size;     
end
clear('ccaxprop','cclengths','ccPprop');
% merge axons
disp([mean(cccv2r(ccsample))   cccv2r(ccsample)'  ])
disp([mean(ccsinu(ccsample))   ccsinu(ccsample)'  ])
disp([mean(ccdiamsm(ccsample)) ccdiamsm(ccsample)'])
if(~fPlot)
ccmsizes   = max(ccszsel,[],1); fprintf('(%d,%d,%d)\n',ccmsizes(1),ccmsizes(2),ccmsizes(3));
else
ccmsizes   = ccszmax;
[nm,jm]    = max(ccmsizes);
% for ji = 1:3; if(ji~=jm); ccmsizes(ji) = 2*ccmsizes(ji); end; end
ccbigvol   = [];
for i = 1:Nnaxons
ccvol = zeros(ccszsel(i,:),'logical'); 
ccax  = zeros(ccmsizes    ,'logical');
ccszdf= floor(( ccmsizes - ccszsel(i,:) ) / 2 ); 
ccvol(ccaxsel(i,1).AxonVoxelIdxList{:}) = true;
ccax(ccszdf(1)+(1:ccszsel(i,1)),ccszdf(2)+(1:ccszsel(i,2)),ccszdf(3)+(1:ccszsel(i,3))) = imclose(ccvol,ones(7,7,7)); 
ccbigvol = cat(2,ccbigvol,ccax);
clear('ccvol','ccax')
end % i
ccidxcut = round(nm/2-cclplotmax/(2*vs)):round(nm/2+cclplotmax/(2*vs));
if(jm == 1); ccbigvol = ccbigvol(ccidxcut,:,:); end
if(jm == 2); ccbigvol = ccbigvol(:,ccidxcut,:); end
if(jm == 3); ccbigvol = ccbigvol(:,:,ccidxcut); end
clear('ccaxsel','ccszsel')
% plot and save
figure(1)
RefineSeg.PlotAxons(imresize3(ccbigvol,2,'Method','linear'),samples{f},fslope);
clear('ccbigvol');
pause(1)
SaveFigurePNG([fOut res{ir} '_' sample '_cc_axons.png'],iwsize,ihsize)
pause(1)
SaveFigureSVG([fOut res{ir} '_' sample '_cc_axons.svg'],iwsize,ihsize)
% pause(1)
close(1)
end % if
end % cc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcg)
% selecting axons only in cc region   
if(ir == 2)
cglbls     = lblscccg.cg_lbl(lblscccg.cg_lbl<=naxsel(f));
else
cglbls     = [];    
end
cgaxprop1  = axprop(cglbls,1:2); 
cgNzo      = Nzo(cglbls);
cgPmany    = Pmany(cglbls);
cgPprop    = Pprop(cglbls);
clear('cglbls')
% 1st selection: remove very small axons and axons with many multiple skeletons per section
cgidx      = ( cgNzo > (2*nrem+2) ) & ( cgPmany < maxPmany ); % & (cgPprop > minPprop); 
clear('cgNzo','cgPmany');
cgaxprop   = cgaxprop1(cgidx,:); 
cgPprop    = cgPprop(cgidx);
clear('cgaxprop1','cgidx');
% 2nd selection: select base in their length, diameter, cv2r and sinuosity
cgdiamsm   = cell2mat(cgdiams(f,'diameteralignedmeans').Variables)*vs;
cgcv2r     = cell2mat(cgdiams(f,'cv2ralign').Variables);
cgsinu     = cell2mat(cglength(f,'sinuosity').Variables);
cglengths  = cell2mat(cglength(f,'euclideanlength').Variables)*vs;
cgidx2     = (cglengths > lplot) & (cgsinu < quantile(cgsinu,1-qfact)) & ...
             (cgcv2rm(f)-sfact*cgcv2rs(f) < cgcv2r & cgcv2r < cgcv2rm(f)+sfact*cgcv2rs(f));
% 3rd selection: closests to means
cgidx2     = cgidx2 & cgPprop > minPprop;
cgind      = find(cgidx2); 
cglen      = cglengths(cgidx2);
if(fcvr)
cgcvdif    = abs(cgcv2rm(f) - cgcv2r(cgidx2)); % new
else
cgcvdif    = abs(cgcv2rm(f)-cgcv2r(cgidx2))+2*abs(cgsinum(f)-cgsinu(cgidx2))+abs(cgdiamm(f)-cgdiamsm(cgidx2)); % new
end
fprintf('cg: %d - (%f) %f %f %f - %f %f, %f %f %f\n',length(cgind),cgdawm(f),cgcv2rm(f),cgsinum(f),cgdiamm(f),...
min(cglen),max(cglen),mean(cgdiamsm(cgind)),mean(cgcv2r(cgind)),mean(cgsinu(cgind)))
clear('cglen');
% [~,cgsidx] = sort(cglen,'descend');
[~,cgsidx] = sort(cgcvdif,'ascend'); % new
Nnaxons    = min( [Naxons length(cgsidx)] );
cgsample   = cgind(cgsidx(1:Nnaxons)); 
cgaxsel    = cgaxprop(cgsample,1); 
if(~isnumeric(cgaxprop(cgsample,2).size))
cgszsel    = cell2mat(cgaxprop(cgsample,2).size); 
else
cgszsel    = cgaxprop(cgsample,2).size;     
end
clear('cgaxprop','cglengths','cgPprop');
% merge axons
disp([mean(cgcv2r(cgsample))   cgcv2r(cgsample)'  ])    
disp([mean(cgsinu(cgsample))   cgsinu(cgsample)'  ]) 
disp([mean(cgdiamsm(cgsample)) cgdiamsm(cgsample)']) 
if(~fPlot)
cgmsizes   = max(cgszsel,[],1); fprintf('(%d,%d,%d)\n',cgmsizes(1),cgmsizes(2),cgmsizes(3));
else
cgmsizes   = cgszmax;
[nm,jm]    = max(cgmsizes);
% for ji = 1:3; if(ji~=jm); cgmsizes(ji) = 2*cgmsizes(ji); end; end
cgbigvol   = [];
for i = 1:Nnaxons
cgvol = zeros(cgszsel(i,:),'logical'); 
cgax  = zeros(cgmsizes    ,'logical'); 
cgszdf= floor(( cgmsizes - cgszsel(i,:) ) / 2 ); 
cgvol(cgaxsel(i,1).AxonVoxelIdxList{:}) = true;
cgax(cgszdf(1)+(1:cgszsel(i,1)),cgszdf(2)+(1:cgszsel(i,2)),cgszdf(3)+(1:cgszsel(i,3))) = imclose(cgvol,ones(7,7,7)); 
cgbigvol = cat(2,cgbigvol,cgax);
clear('cgvol','cgax')
end % i
cgidxcut = round(nm/2-cglplotmax/(2*vs)):round(nm/2+cglplotmax/(2*vs));
if(jm == 1); cgbigvol = cgbigvol(cgidxcut,:,:); end
if(jm == 2); cgbigvol = cgbigvol(:,cgidxcut,:); end
if(jm == 3); cgbigvol = cgbigvol(:,:,cgidxcut); end
clear('cgaxsel','cgszsel')
% plot and save
figure(1)
RefineSeg.PlotAxons(imresize3(cgbigvol,2,'Method','linear'),samples{f},fslope);
clear('cgbigvol');
pause(1)
SaveFigurePNG([fOut res{ir} '_' sample '_cg_axons.png'],iwsize,ihsize)
pause(1)
SaveFigureSVG([fOut res{ir} '_' sample '_cg_axons.svg'],iwsize,ihsize)
% pause(1)
close(1)
end % if
end % cg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',sample)
end % f
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('End\n')
end % main
