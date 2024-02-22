function script10_save_computed_disp()
% Concatenate FOD and dispersion properties for all samples
%
% This script needs the results from script:
% script09_fods_dispersion.m
% and the information in the .mat files in the folder
% cc_cg_labels/
% The results of this script will be saved in folder
% White_matter_EM_03_Dispersion/*/
% 26-May-2023 by Ricardo Coronado-Leija
close all
clc
rng(0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% folders needed
fprintf('Begin\n')
ir          = 2; % 1 = HM, 2 = LM
res         = {'HM','LM'};
fIn         = ['White_Matter_EM_FOD_' res{ir} '/']; 
folderqn    = 'cc_cg_labels/';
fOut        = 'White_matter_EM_03_Dispersion/';
if(~isfolder(fOut)); mkdir(fOut); end
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
% dtime       = 0;
dtime       = 11.5;
% dtime       = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f = nsamples:-1:1
sample = [num2str(idsample(f)) '_' sidesample{f}];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read labels and set flags
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
if(fcc && fcg)
% read results   
ccggname = [fIn '/' samples{f} '_cccg_fod_t' num2str(dtime) '.mat'];
cccgfod  = load(ccggname).fod;

% set values
cccgdispersion(f).time         = cccgfod.t; 
cccgdispersion(f).vs           = vs; 
% fod 
cccgdispersion(f).fod          = cccgfod.fodi;
cccgdispersion(f).points       = cccgfod.points;
cccgdispersion(f).dir          = cccgfod.maindir';
cccgdispersion(f).adisp        = cccgfod.dispang;
% spherical harmonics (largest order)
cccgdispersion(f).lmax         = cccgfod.lmax2;
cccgdispersion(f).plm          = cccgfod.plm2;
cccgdispersion(f).pl           = cccgfod.pl2;
cccgdispersion(f).C            = cccgfod.C2;
cccgdispersion(f).lambda       = cccgfod.lambda2;
cccgdispersion(f).adisp_p2     = cccgfod.dispangp2;

end % cccg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcc) 
% read results    
ccname = [fIn '/' samples{f} '_cc_fod_t' num2str(dtime) '.mat'];
ccfod  = load(ccname).fod;    
 
% set values
ccdispersion(f).time         = ccfod.t; 
ccdispersion(f).vs           = vs; 
% fod 
ccdispersion(f).fod          = ccfod.fodi;
ccdispersion(f).points       = ccfod.points;
ccdispersion(f).dir          = ccfod.maindir';
ccdispersion(f).adisp        = ccfod.dispang;
% spherical harmonics (largest order)
ccdispersion(f).lmax         = ccfod.lmax2;
ccdispersion(f).plm          = ccfod.plm2;
ccdispersion(f).pl           = ccfod.pl2;
ccdispersion(f).C            = ccfod.C2;
ccdispersion(f).lambda       = ccfod.lambda2;
ccdispersion(f).adisp_p2     = ccfod.dispangp2;
 
end % cc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcg)
% read results    
cgname = [fIn '/' samples{f} '_cg_fod_t' num2str(dtime) '.mat'];
cgfod  = load(cgname).fod;  
 
% set values
cgdispersion(f).time         = cgfod.t; 
cgdispersion(f).vs           = vs; 
% fod  
cgdispersion(f).fod          = cgfod.fodi;
cgdispersion(f).points       = cgfod.points;
cgdispersion(f).dir          = cgfod.maindir';
cgdispersion(f).adisp        = cgfod.dispang;
% spherical harmonics
cgdispersion(f).lmax         = cgfod.lmax2;
cgdispersion(f).plm          = cgfod.plm2;
cgdispersion(f).pl           = cgfod.pl2;
cgdispersion(f).C            = cgfod.C2;
cgdispersion(f).lambda       = cgfod.lambda2;
cgdispersion(f).adisp_p2     = cgfod.dispangp2;
 
end % cg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',sample)
end % f
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc_dispersion = struct2table(ccdispersion);
clear('ccdispersion')
tic
ccname = [fOut '/cc_dispersion_' res{ir} '.mat'];
save(ccname,'cc_dispersion'); 
clear('cc_dispersion')
%
if(ir == 2)
cg_dispersion = struct2table(cgdispersion);
clear('cgdispersion')
cgname = [fOut '/cg_dispersion_' res{ir} '.mat'];
save(cgname,'cg_dispersion'); 
clear('cg_dispersion')
%
cccg_dispersion = struct2table(cccgdispersion);
clear('cccgdispersion')
cccgname = [fOut '/cccg_dispersion_' res{ir} '.mat'];
save(cccgname,'cccg_dispersion'); 
clear('cccg_dispersion')
end
toc
fprintf('End\n')
end % main