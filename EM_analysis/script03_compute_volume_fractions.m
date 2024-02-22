function script03_compute_volume_fractions()
% Compute volume fractions 
% This script needs results from script:
% script02_proofread_vols.m
% It also needs the information in the text file
% 000_animal_exp_axons_cells_nx_ny_nz_cclims_cglims.txt
% The results of this script will be saved in folder
% White_matter_EM_01_Volume_Props/*/
%
% 02-Sep-2022 by Ricardo Coronado-Leija
close all
clc
rng(0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% folders needed
fprintf('Begin\n')
foldersegs  = 'White_matter_EM/';
foldervols  = 'White_matter_EM_FullSubstrates/';
folderqn    = '';
fOut        = 'White_matter_EM_01_Volume_Props/';
if(~isfolder(fOut)); mkdir(fOut); end
% name of files needed
volaxons    = 'myelinated_axons.mat';
volcells    = 'nucleus.mat';
volcleaned  = 'substrate_cleaned.mat';
sufixqn     = '000_animal_exp_axons_cells_nx_ny_nz_cclims_cglims.txt';
% samples
samples     = {'Sham_25_contra','Sham_25_ipsi','Sham_49_contra','Sham_49_ipsi',...
               'TBI_24_contra','TBI_24_ipsi','TBI_28_contra','TBI_28_ipsi',...
               'TBI_2_contra','TBI_2_ipsi'};
idsample    = [25 25 49 49 24 24 28 28 2 2 ]; 
sidesample  = {'contra','ipsi','contra','ipsi','contra','ipsi',...
               'contra','ipsi','contra','ipsi'};
nsamples    = length(samples);
ir          = 2; % 1 = HM, 2 = LM
res         = {'HM','LM'};
voxsizes    = [15 50]; % now volumes are isotropic
vs          = voxsizes(ir)*1e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f = nsamples:-1:1
sample = [num2str(idsample(f)) '_' sidesample{f}];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data
tic
soname = [foldervols samples{f} '/' res{ir} '_' sample '_' volcleaned];
SubstrateCleaned  = load(soname).SubstrateCleaned;  
lbname = [foldersegs samples{f} '/' res{ir} '_' sample '_' volaxons]; 
ncname = [foldersegs samples{f} '/' res{ir} '_' sample '_' volcells];
if(ir == 2)
% big volumes: low resolution    
Axons = load(lbname).final_lbl; 
Cells = load(ncname).final_lbl; 
else
% small volumes: high resolution    
Axons = load(lbname).myelinated_axons; 
Cells = [];
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ir == 2)
is     = f;
props  = load([folderqn sufixqn]); % samples properties and cc,cg limits
propsi = props(:,is);
vcc    = length(propsi(8):propsi(9))*length(propsi(10):propsi(11))*length(propsi(12):propsi(13));
vcg    = length(propsi(14):propsi(15))*length(propsi(16):propsi(17))*length(propsi(18):propsi(19));
fcc    = vcc > 10; % nx*ny*nz > 10 big volume to exclue 49cgci
fcg    = vcg > 10; % nx*ny*nz > 10 big volume to exclue 49cgci
else
fcc    = 1;
fcg    = 0;
propsi = zeros(20,1);
[nx,ny,nz] = size(SubstrateCleaned);
propsi(8:13) = [1 nx 1 ny 1 nz];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% divide the volumes into regions
if(fcc && fcg)
CCCG_vol  = SubstrateCleaned;
CCCGaxons = Axons;
if(ir == 2)
CCCGcells = Cells; 
end
end
if(fcc) 
CC_col = SubstrateCleaned(propsi(8):propsi(9),propsi(10):propsi(11),propsi(12):propsi(13));
if(ir == 2)
CCcells = Cells           (propsi(8):propsi(9),propsi(10):propsi(11),propsi(12):propsi(13));
CCaxons = Axons           (propsi(8):propsi(9),propsi(10):propsi(11),propsi(12):propsi(13));
else
CCaxons = Axons           (propsi(8):propsi(9),propsi(10):propsi(11),propsi(12):285); % all have nz=285
end
end
if(fcg)
CG_vol  = SubstrateCleaned(propsi(14):propsi(15),propsi(16):propsi(17),propsi(18):propsi(19));
CGaxons = Axons           (propsi(14):propsi(15),propsi(16):propsi(17),propsi(18):propsi(19));
if(ir == 2)
CGcells = Cells           (propsi(14):propsi(15),propsi(16):propsi(17),propsi(18):propsi(19));
end
end
clear('SubstrateCleaned');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcc && fcg)
% ==== Volumes ===== %    
cccg.vs       = vs;
cccg.sz       = size(CCCG_vol);
cccg.vol      = numel(CCCG_vol);
% axons volume
cccg.ias_vol = sum(CCCG_vol == 2,'all');
% myelin volume 
cccg.my_vol  = sum(CCCG_vol == 0,'all');
% nucleus volume
if(ir == 2)
cccg.nuc_vol = sum(CCCG_vol == 3,'all');
else
cccg.nuc_vol = 0;
end
clear('CCCG_vol')
% extra-axonal volume
cccg.eas_vol   = cccg.vol - cccg.ias_vol - cccg.my_vol;
% ==== Fractions ==== %
% axon fraction
cccg.iavf_vol  = cccg.ias_vol/cccg.vol;
% myelin fraction
cccg.mvf_vol   = cccg.my_vol/cccg.vol;
% nucleus fraction
cccg.nvf_vol   = cccg.nuc_vol/cccg.vol;
% extra-axonal fraction
cccg.eavf_vol  = cccg.eas_vol/cccg.vol;
% ==== dMRI ==== %
% intra axonal water fraction
cccg.f_vol     = cccg.iavf_vol/(cccg.iavf_vol+cccg.eavf_vol);
% ==== Counts ==== %
% number axons cells
cccg.naxon = length( unique( CCCGaxons(:) ) );   
if(ir == 2)
cccg.ncell = length( unique( CCCGcells(:) ) ); 
else
cccg.ncell = 0;    
end
clear('CCCGaxons','CCCGcells') 
% concatenate into a structure
CCCGvolprops(f) = cccg;
end % cgcg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcc) 
% ==== Volumes ===== %    
cc.vs       = vs;
cc.sz       = size(CC_col);
cc.vol      = numel(CC_col);
% axons volume
cc.ias_vol = sum(CC_col == 2,'all');
% myelin volume 
cc.my_vol  = sum(CC_col == 0,'all');
% nucleus volume
if(ir == 2)
cc.nuc_vol = sum(CC_col == 3,'all');
else
cc.nuc_vol = 0;
end
clear('CC_col')
% extra-axonal volume
cc.eas_vol   = cc.vol - cc.ias_vol - cc.my_vol;
% ==== Fractions ==== %
% axon fraction
cc.iavf_vol  = cc.ias_vol/cc.vol;
% myelin fraction
cc.mvf_vol   = cc.my_vol/cc.vol;
% nucleus fraction
cc.nvf_vol   = cc.nuc_vol/cc.vol;
% extra-axonal fraction
cc.eavf_vol  = cc.eas_vol/cc.vol;
% ==== dMRI ==== %
% intra axonal water fraction
cc.f_vol     = cc.iavf_vol/(cc.iavf_vol+cc.eavf_vol);
% ==== Counts ==== %
% number axons cells
cc.naxon = length( unique( CCaxons(:) ) );   
if(ir == 2)
cc.ncell = length( unique( CCcells(:) ) ); 
else
cc.ncell = 0;    
end
clear('CCaxons','CCcells')
% concatenate into a structure
CCvolprops(f) = cc;
end % cc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcg)
% ==== Volumes ===== %    
cg.vs       = vs;
cg.sz       = size(CG_vol);
cg.vol      = numel(CG_vol);
% axons volume
cg.ias_vol = sum(CG_vol == 2,'all');
% myelin volume 
cg.my_vol  = sum(CG_vol == 0,'all');
% nucleus volume
if(ir == 2)
cg.nuc_vol = sum(CG_vol == 3,'all');
else
cg.nuc_vol = 0;
end
clear('CG_vol')
% extra-axonal volume
cg.eas_vol   = cg.vol - cg.ias_vol - cg.my_vol;
% ==== Fractions ==== %
% axon fraction
cg.iavf_vol  = cg.ias_vol/cg.vol;
% myelin fraction
cg.mvf_vol   = cg.my_vol/cg.vol;
% nucleus fraction
cg.nvf_vol   = cg.nuc_vol/cg.vol;
% extra-axonal fraction
cg.eavf_vol  = cg.eas_vol/cg.vol;
% ==== dMRI ==== %
% intra axonal water fraction
cg.f_vol     = cg.iavf_vol/(cg.iavf_vol+cg.eavf_vol);
% ==== Counts ==== %
% number axons cells
cg.naxon = length( unique( CGaxons(:) ) );   
if(ir == 2)
cg.ncell = length( unique( CGcells(:) ) ); 
else
cg.ncell = 0;    
end
clear('CGaxons','CGcells')
% concatenate into a structure
CGvolprops(f) = cg;
end % cg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',sample)
end % f
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc_vol_properties = struct2table(CCvolprops);
clear('CCvolprops')
tic
ccname = [fOut '/cc_vol_properties_' res{ir} '.mat'];
save(ccname,'cc_vol_properties'); 
clear('cc_vol_properties')
%
if(ir == 2)
cg_vol_properties = struct2table(CGvolprops);
clear('CGvolprops')
cgname = [fOut '/cg_vol_properties_' res{ir} '.mat'];
save(cgname,'cg_vol_properties'); 
clear('cg_vol_properties')
%
cccg_vol_properties = struct2table(CCCGvolprops);
clear('CCCGvolprops')
cccgname = [fOut '/cccg_vol_properties_' res{ir} '.mat'];
save(cccgname,'cccg_vol_properties');  
clear('cccg_vol_properties')
toc
end
fprintf('End\n')
end % main