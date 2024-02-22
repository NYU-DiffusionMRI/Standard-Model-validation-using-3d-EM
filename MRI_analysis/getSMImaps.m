function getSMImaps(pathIn,pathScheme,pathSigma,pathMask,folderOut,Basename,smioptions)
% Fit and get Standard Model parameter maps using the SMI toolbox (fixed TE )
%
% pathIn: path to dwi (.nii) file     
% pathScheme: path to protocol (.txt) file (MRtrix style)
% pathSigma: sigma noise map (.nii), such as obtained from MP-PCA during pre-processing
% pathMask: path to mask (.nii) file
% folderOut: output directory
% Basename: output basename
% smioptions: flags/values for SMI fitting 
%         lmax:   lmax for spherical harmonics
%         Dmax:   max diffusivity for Da, Depar and Deperp (override by ul)
%         Diso:   isotropic diffusivity (e.g. free water)
%         fiso:   flag (0/1), estimate isotropic fraction? (e.g. free water fraction)
%         rc:     flag (0/1), apply rician bias correction?
%         fod:    flag (0/1), compute fod using spherical deconvolution?
%         fodseg: flag (0/1), segment fod lobes ? very slow
%         ul:     two row array indicating the lower (1,:) and upper (2,:) 
%                 bounds for the parameters (uniform distribution prior)
% This script calls the SMI toolbox:
% https://github.com/NYU-DiffusionMRI/SMI
% It needs tools for reading/writting nifti files
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% By Ricardo Coronado-Leija
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max spherical harmonics
if(isfield(smioptions,'lmax') )
lmax = smioptions.lmax;   
else
lmax = 4; % dki protocol    
end % lmax
% max diffusivity (Da, Depar and Deperp)  
if(isfield(smioptions,'Dmax') )
Dmax = smioptions.Dmax;   
else
Dmax = 3;   
end % Dmax
% isotropic diffusivity (e.g. free water)  
if(isfield(smioptions,'Diso') )
Diso = smioptions.Diso;   
else
Diso = 3;   
end % Diso
% flag(1/0): estimate isotropic compartment fraction? (e.g. free water fraction)  
if(isfield(smioptions,'fiso') )
fiso = smioptions.fiso;   
else
fiso = 0;   
end % fiso
% flag(1/0): use rician bias correction?  
if(isfield(smioptions,'rc') )
rc = smioptions.rc;   
else
rc = 0;   
end % rc
% flag(1/0): compute full FOD?  
if(isfield(smioptions,'fod') )
fod = smioptions.fod;   
else
fod = 0;   
end % fod
% flag(1/0): segment FOD? very very slow 
if(isfield(smioptions,'fodseg') )
fodseg = smioptions.fodseg;   
else
fodseg = 0;   
end % fodseg
% upper and lower bounds for all parameters
if(isfield(smioptions,'ul') )
ul     = smioptions.ul;   
ulflag = 1;
else
ulflag = 0;    
end % ul
% create output folder
if(~isfolder(folderOut)); mkdir(folderOut); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMask = 1;
% read 4D volume (nii)
tic
diw_nii      = load_untouch_nii(pathIn); 
dwi          = single(diw_nii.img);
[nr,nc,ns,~] = size(dwi);
clear('diw_nii');
% read scheme
scheme       = readmatrix(pathScheme);
bvec         = double(scheme(:,1:3));
bval         = double(scheme(:,4));
% use or not the mask
if(fMask && ~isempty(pathMask))
% read mask
niimask      = load_untouch_nii(pathMask);    
mask         = niimask.img > 0.1;
clear('niimask')
else
% mask includes all voxels    
mask         = true(nr,nc,ns); 
pathMask     = pathIn;
end
% read sigma
niisigma = load_untouch_nii(pathSigma); 
sigma    = niisigma.img;
clear('niisigma')
% Remove negative/zero values
fprintf('Correcting\n')
mindw    = min(dwi(dwi > 0));
dwi(dwi <= 0) = mindw;
dwi      = removebad(dwi,mindw);
idx0     = bval < 50; 
S0m      = mean(dwi(:,:,:,idx0) ,4);
fprintf('voxels %d\n',sum(mask(:)));
t = toc;
fprintf('Time read data %f s\n',t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bval = round(bval/1000); %1000*
bval = bval/1000; 

% % normalization ??
% dwi   = dwi./S0m;
% sigma = sigma./S0m;

% Specify mask and noise map
options.mask  = mask > 0;
options.sigma = abs(sigma);

% Specify protocol information (no need to specify beta and TE if all
% measurements are LTE and have the same TE)
options.b    = bval;
options.beta = [];
options.dirs = bvec;
options.TE   = [];
options.MergeDistance = 0.1; % If []: default is 0.05 [ms/um^2], this 
% is the threshold for considering different b-values as the same shell

% Specify options for the fit
if(fiso == 1)
options.compartments = {'IAS','EAS','FW'}; % The order does not matter
ffw_max = 0.50;
else
options.compartments = {'IAS','EAS'};      % The order does not matter
ffw_max = 0.01;
end
%
% Removing Rician Bias
options.NoiseBias = 'None';
if(rc == 1)
% options.NoiseBias    = 'Rician'; % this data has Rician noise bias
val = dwi.^2 - repmat(sigma.^2,[1 1 1 length(bval)]);
val(val < 0) = 0;
dwi = sqrt( val );   
% else
% options.NoiseBias    = 'None';   % this data does not have Rician noise bias    
end

if(ulflag)
options.MLTraining.bounds = ul;
else
options.MLTraining.bounds = [0.01, 0.01, 0.01, 0.01,       0,  50,  30, 0.00; ...
                             0.99, Dmax, Dmax, Dmax, ffw_max, 150, 100, 1.00];
end
% options.MLTraining.bounds = [0.05, 1, 1, 0.1, 0  ,  50,  50, 0.05; ...
%                              0.95, 3, 3, 1.2, 0.5, 150, 120, 0.99];

options.MLTraining.Ntraining = 1e5;

plflag = lmax;
options.RotInv_Lmax=plflag;

options.D_FW = Diso;

options.Lmax=lmax; 

% options.Nlevels = 64;
options.Nlevels = 31;

% FOD
options.flag_fit_fODF = fod; 

% Condon-Shortley phase factor should be zero to visualize fod in mrtrix
options.CS_phase = 0; 

% Run SM fitting
tic
[out] = SMI.fit(dwi,options);
t=toc;
fprintf('Time SM fit %f s\n',t)

% 
out.RotInvs.S0 = removebad(out.RotInvs.S0);
out.RotInvs.S2 = removebad(out.RotInvs.S2);
if(plflag > 2)
out.RotInvs.S4 = removebad(out.RotInvs.S4);
end
if(plflag > 4)
out.RotInvs.S6 = removebad(out.RotInvs.S6);
end
% 
out.kernel = removebad(out.kernel);
nkernel    = size(out.kernel,4);

% ======================================================================= %
% ======================================================================= %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================================================================= %

% Save SM parameters
tag_parameter={'$f$','$D_a$','$D_e^{||}$','$D_e^{\perp}$','$f_{FW}$','$p_2$','$p_4$','$p_6$'};
% kernel mask pdf 
Ylims=[0 1;0 3;0 3;0 2;0 0.5;0 1;0 1;0 1];
fig=figure('Name','1','Position',[241 343 1894 759]);
namepar = {'f','D_a','D_e_par','D_e_perp','ffw','p2','p4','p6'};
fprintf('Dfw=%.2f:\n',Diso)
for ii=1:nkernel
current_param=out.kernel(:,:,:,ii);
current_param(~mask) = 0;
% pdf
[pdf_fit,x_fit]=ksdensity(current_param(reshape(mask(:)==1,[],1)));
currLims=[0 1.1*max([pdf_fit(:)])];
mu_fits=mean(current_param(reshape(mask(:)==1,[],1)))*[1 1];
mu_prior=mu_fits;%mean(prior(:,ii))*[1 1];
subplot(2,4,ii), plot(x_fit,pdf_fit,mu_fits,currLims,'r-.',mu_prior,currLims,'b-.','LineWidth',2),set(gca,'FontSize',25)
title(tag_parameter{ii},'interpreter','latex'), xlim(Ylims(ii,:))
SaveNII(current_param ,[folderOut '/' Basename '_' namepar{ii} '.nii.gz'] , pathMask,[]);
fprintf('%s=%.2f ',namepar{ii},mu_fits(1))
end
fprintf('\n')
legend('fits','\mu fits','\mu prior')
sgt = sgtitle('Kernel estimates - mask'); sgt.FontSize = 30;
saveas(fig,[folderOut '/' Basename '_sm_fits.png'])
close('1')

% Save RotInvs
SaveNII(S0m,[folderOut '/' Basename '_b0.nii.gz'] , pathMask,[]);
SaveNII(out.RotInvs.S0(:,:,:,1),[folderOut '/' Basename '_S0m.nii.gz'] , pathMask,[]);
SaveNII(out.RotInvs.S0./out.RotInvs.S0(:,:,:,1),[folderOut '/' Basename '_S0.nii.gz'] , pathMask,[]);
SaveNII(out.RotInvs.S2./out.RotInvs.S0(:,:,:,1),[folderOut '/' Basename '_S2.nii.gz'] , pathMask,[]);
if(lmax > 2)
SaveNII(out.RotInvs.S4./out.RotInvs.S0(:,:,:,1),[folderOut '/' Basename '_S4.nii.gz'] , pathMask,[]);
end
if(plflag > 4)
SaveNII(out.RotInvs.S6./out.RotInvs.S0(:,:,:,1),[folderOut '/' Basename '_S6.nii.gz'] , pathMask,[]);
end

% Power law (smi)
if(plflag > 2)
% polynomial regression    
plpr   = out.kernel(:,:,:,6:end);
plpr2D = vectorize_mask(plpr,mask);
lplpr  = log(plpr2D);
if(plflag <= 4)
lplpr = lplpr(1:2,:);
end
% matrices for power law
if(plflag > 4)
Apl  = [1 2; 1 4; 1 6];
else
Apl  = [1 2; 1 4];    
end
% polynomial regression
xpr      = Apl\lplpr;
plpre    = vectorize_mask( exp(Apl*xpr)  ,mask);
Cpr      = vectorize_mask( exp(xpr(1,:)) ,mask);
lambdapr = vectorize_mask( exp(xpr(2,:)) ,mask);
SaveNII(Cpr     ,[folderOut '/' Basename '_C_pr.nii.gz']      , pathMask,[]);
SaveNII(lambdapr,[folderOut '/' Basename '_lambda_pr.nii.gz'] , pathMask,[]);
SaveNII(plpre   ,[folderOut '/' Basename '_ple_pr.nii.gz']    , pathMask,[]);
end 
% ======================================================================= %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================================================================= %
if(fod)
% fod and its rotinv    
plm  = cat(4,double(mask>0),out.plm);
pl   = cat(4,double(mask>0),out.pl);
% compute angles of dispersion   
p2   = out.kernel(:,:,:,6); % extract p2 from kernel
adp2 = acosd(sqrt(            2*p2/3+1/3));
adpl = acosd(sqrt(   2*pl(:,:,:,2)/3+1/3));
% Save FOD SH
SaveNII(plm   ,[folderOut '/' Basename '_plm.nii.gz']    , pathMask,[]);
% Save FOD RotInv 
SaveNII(pl    ,[folderOut '/' Basename '_pl.nii.gz']     , pathMask,[]);
% Save PLR Dispersion Angle  
SaveNII(adp2  ,[folderOut '/' Basename '_p2_adisp.nii.gz']   , pathMask,[]);
% Save FOD Dispersion Angle  
SaveNII(adpl  ,[folderOut '/' Basename '_pl_adisp.nii.gz']   , pathMask,[]);
% ----------------------------------------
if(fodseg)
% points
points        = load('PuntosElectroN6000.txt');
% inner products
pps           = single(points*points');
% duplicate
points2       = [points; -points];
% np2           = 2*np;
pps2          = [pps -pps; -pps pps]; 
% pps2          = points2*points2';
% imshow(pps2-[pps -pps; -pps pps])
[spps2,ipps2] = sort(abs(pps2),2,'descend');
% ============== Create SH Basis ============== %    
% Spherical Harmonics 
CSph     = 0;
SH_D     = SMI.get_even_SH(points,lmax,CSph);
SH_D2    = [SH_D; SH_D];
Nl       = SHnorms(lmax);    
plmallD  = vectorize_mask(plm,mask > 0); 
nvall    = sum(mask(:)>0);
% remove normalization plms, so normalization ( 4*pi*fod/length(fod) ) = 1
plmallD  = plmallD.*(1./Nl);    
%
nmax      = 3;     % max number of peaks
minAng2   = 25;    % angle for neighboor 
tolAng    = 0.1;   % tolerance angle
tic
[peaksseg,peakampseg,volfrac1seg,~,plmsegseg,plssegseg] = ... 
SegmentFODsVoxels(plmallD,points2,spps2,ipps2,minAng2,tolAng,SH_D2,SH_D2,nmax);
toc
% select just the ones in the mask that contains all
peaksall    = reshape( permute(peaksseg(:,:,1:nvall),[2 1 3]) , [ nmax*3 nvall] ); % (x1,y1,z1,x2,y2,z2,x3,y3,z3)
peakampall  = peakampseg(:,1:nvall);
volfrac1all = volfrac1seg(:,1:nvall);
plmsegall   = plmsegseg(:,1:nvall,:).*(Nl);
plssegall   = plssegseg(:,1:nvall,:);
adispall    = squeeze( acosd(sqrt(2*plssegall(2,:,:)/3+1/3)) )';
% get volumes
adispvol    = vectorize_mask(adispall   ,mask>0);
peaksvol    = vectorize_mask(peaksall   ,mask>0);
peakampvol  = vectorize_mask(peakampall ,mask>0);
volfrac1vol = vectorize_mask(volfrac1all,mask>0);
plmsegpk1   = vectorize_mask(plmsegall(:,:,1),mask>0);
plmsegpk2   = vectorize_mask(plmsegall(:,:,2),mask>0);
plmsegpk3   = vectorize_mask(plmsegall(:,:,3),mask>0);
plssegpk1   = vectorize_mask(plssegall(:,:,1),mask>0);
plssegpk2   = vectorize_mask(plssegall(:,:,2),mask>0);
plssegpk3   = vectorize_mask(plssegall(:,:,3),mask>0);
%
SaveNII(adispvol   ,[folderOut '/' Basename '_peaks_p2_adisp.nii.gz'], pathMask,[]);
SaveNII(peaksvol   ,[folderOut '/' Basename '_peaks.nii.gz']         , pathMask,[]);
SaveNII(peakampvol ,[folderOut '/' Basename '_peaks_amp.nii.gz']     , pathMask,[]);
SaveNII(volfrac1vol,[folderOut '/' Basename '_peaks_vol.nii.gz']     , pathMask,[]);
SaveNII(plmsegpk1  ,[folderOut '/' Basename '_peak1_plm.nii.gz']     , pathMask,[]);
SaveNII(plmsegpk2  ,[folderOut '/' Basename '_peak2_plm.nii.gz']     , pathMask,[]);
SaveNII(plmsegpk3  ,[folderOut '/' Basename '_peak3_plm.nii.gz']     , pathMask,[]);
SaveNII(plssegpk1  ,[folderOut '/' Basename '_peak1_pls.nii.gz']     , pathMask,[]);
SaveNII(plssegpk2  ,[folderOut '/' Basename '_peak2_pls.nii.gz']     , pathMask,[]);
SaveNII(plssegpk3  ,[folderOut '/' Basename '_peak3_pls.nii.gz']     , pathMask,[]);
end % fodseg
% ----------------------------------------
% Power law
if(lmax > 2)
% spherical deconvolution
plsd2D    = vectorize_mask(pl   (:,:,:,2:end),mask);
lplsd     = log(plsd2D);
if(lmax <= 4)
lplsd    = lplsd(1:2,:);
end
% matrices for power law
ls   = (2:2:lmax)';
Apl  = [ones(size(ls)) ls];
% spherical deconvolution
xsd         = Apl\lplsd;
plsde       = vectorize_mask( exp(Apl*xsd)     ,mask);
Csd         = vectorize_mask( exp(xsd(1,:))    ,mask);
lambdasd    = vectorize_mask( exp(xsd(2,:))    ,mask);
SaveNII(Csd        ,[folderOut '/' Basename '_C_sd.nii.gz']         , pathMask,[]);
SaveNII(lambdasd   ,[folderOut '/' Basename '_lambda_sd.nii.gz']    , pathMask,[]);
SaveNII(plsde      ,[folderOut '/' Basename '_ple_sd.nii.gz']       , pathMask,[]);
end 
end
end % main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function volout = removebad(volin,val)
if(nargin == 1)
val = 0;    
end
volout = volin;
volout(isnan(volin) | isinf(volin)) = val;
end % removebad

function lmax = nsh2lmax(nsh)
lmax = 2*(floor((sqrt(1+8*nsh)-3)/4));
end

function nSH = lmax2nsh(lmax)
nSH = ((lmax+1).*(lmax+2)/2);
end

function nzh = lmax2nzh(lmax)
nzh = floor(lmax/2) + 1;
end

function N_l = SHnorms(lmax,feven)
if(nargin < 2)
feven = 1;
end
if(feven)    
L_unique = 0:2:lmax;
else
L_unique = 0:lmax;    
end
L_all    = repelem(L_unique,2*L_unique+1)';
N_l      = sqrt((2*L_all+1)*(4*pi));
end % SHnorms

function pl = plm2pl(plm)
% compute rotational invariants of plm 
% plm should be normalized (plm(0)~=1)

[nsh,nv] = size(plm);
lmax     = nsh2lmax(nsh);
nZH      = lmax2nzh(lmax);

pl = zeros(nZH,nv);
k      = 0;
for l = 0:2:lmax
m    = -l:l;
nm   = length(m);
list = (k+1):(k+nm);
pl(floor(l/2)+1,:) = sqrt(sum(abs(plm(list,:)).^2,1));
k = k + nm;
end % l

end % plm2pl

function [peaks,peakamp,volfrac1,volfrac2,plmseg,plsseg] = SegmentFODsVoxels(plms,points,SortPPs,iPPs,minang,tolang,SHori,SHnew,nmax)
% shange to only use SHori

[nshori,nvoxels] = size(plms);
nshnew           = size(SHnew,2);
lmaxori          = nsh2lmax(nshori);
nzhori           = lmax2nzh(lmaxori);
lmaxnew          = nsh2lmax(nshnew);
fods             = SHori*plms;%HarmonicsFourierTransform(plms,SHori,0);

peaks    = zeros(nmax,3,nvoxels);
peakamp  = zeros(nmax,nvoxels);
volfrac1 = zeros(nmax,nvoxels);
volfrac2 = zeros(nmax,nvoxels);
plmseg   = zeros(nshnew,nvoxels,nmax);
plsseg   = zeros(nzhori,nvoxels,nmax);

for v = 1:nvoxels
[pk,fv,hi]         = SegmentFODMaxFastOff(fods(:,v),points,SortPPs,iPPs,minang,tolang);
nseg               = min([nmax size(pk,1)]);
[~,plm,pl,vf]      = getFODLobesOff(fods(:,v),hi,nseg,lmaxnew,SHnew);
peaks(1:nseg,:,v)  = pk(1:nseg,:);
peakamp(1:nseg,v)  = fv(1:nseg,2);
volfrac1(1:nseg,v) = fv(1:nseg,1);
volfrac2(1:nseg,v) = vf(1:nseg);
plmseg(:,v,1:nseg) = plm(:,1:nseg);
plsseg(:,v,1:nseg) = pl(1:nzhori,1:nseg);
% fprintf('%d/%d\n',v,nvoxels)
end

end % SegmentFODsVoxels

function [peaks,fodvals,H] = SegmentFODMaxFastOff(fod,points,SortPPs,iPPs,minang,tolang)
% [peaks,fodvals,H] = SegmentFODMaxFast(fod,points,minang,tolang)
% This function takes a fod and the points where its sampled and tries to
% divide it into lobes by clustering its points to the ascociated peaks
% ======================================================================= %
% Input:
% fod     -> fiber orientation distribution
% points  -> points in which fod is sampled
% SortPPs -> if PPs = abs(points*points'), SortPPs have the values sortted
%            by row in descending order (closest to farest)
% iPPs    -> Indicated the original indices of the sorted points
% minang  -> angle of the neighborhoods considered for each point to look
%            for the fod maxima
% tolang  -> tolerance angle between the current and previous maxima, if
%            the angle becomes lower the searching stops
%
% Output:
% H(i,1:3) -> maxima point towards the point points(i,:) is atracted
% H(i,4)   -> maxima fod value
% H(i,5)   -> maxima points index
% H(i,6)   -> cluster
% fodvals  -> sum of fod lobe values , peak fod value
% peaks    -> directions of the 'unique' lobe peaks
% ======================================================================= %
% By Ricardo Coronado Leija: 2023-May-04
% ======================================================================= %

% Difference between this code and SegmentFODMaxFast: considering offset!

% (***) Remove offset
fodoff = min(fod);
fod    = fod - fodoff;

% Set some values
fPlot  = rand < 0; % if true, PlotSpherical() function should be on path
fPrint = rand < 0;
maxit  = 100;
PPtol  = cos(tolang*pi/180);
PPmin  = cos(minang*pi/180);
N      = size(points,1);
H      = zeros(N,6);      % point -> fod -> idx

% crop SortPPs, iPPs to only include neighborhoods with pps > PPmin (ang < minang )
ii      = SortPPs > PPmin; 
jj      = sum(ii,1) > 0;
npp     = sum(jj);
SortPPs = SortPPs(:,1:npp); 
iPPs    = iPPs(:,1:npp); 
iMax    = zeros(N,1);
mMax    = zeros(N,1);

% ======================================================================= %
% ======================================================================= %
% ======================================================================= %

% Find maximums for each neighborhood
for i = 1:N
% find neighborhood with ang < minang    
idx    = iPPs(i,:);
% max fod in neighborhood
[~,mi]  = max(fod(idx));
iMax(i) = mi; 
mMax(i) = idx(mi); 
end

% ======================================================================= %
% ======================================================================= %
% ======================================================================= %

for i = 1:N
% === first step === %    
% original point
p      = points(i,:);
v      = fod(i);
% find neighborhood with ang < minang
PPs    = SortPPs(i,:);
% max new point
mi     = iMax(i);
mm     = mMax(i);
pn     = points(mm,:);
vn     = fod(mm);
% pp of new vs old point
pp     = PPs(mi);
% ================== %
%
if(fPrint)
fprintf('%d-(%.2f %.2f %.2f) = %.3f -> (%.2f %.2f %.2f) = %.3f\n',i,p(1),p(2),p(3),v,pn(1),pn(2),pn(3),vn);
end
if(fPlot)
figure(1); hold on
PlotSpherical(fod,points,0.1);   
scatter3(v*p(1)  ,v*p(2)  ,v*p(3)  ,100,'k');
scatter3(vn*pn(1),vn*pn(2),vn*pn(3),100,'k','filled');
end
%
% === start iterative process === %
niter = 0;
while(pp < PPtol && niter < maxit)
% current point    
p     = pn;
v     = vn;
% find neighborhood with ang < minang
PPs   = SortPPs(mm,:);
% max new point
mi    = iMax(mm);
mm    = mMax(mm);
pn    = points(mm,:);
vn    = fod(mm);
% change
pp    = PPs(mi);
niter = niter + 1;
%
if(fPrint)
fprintf('%d-(%.2f %.2f %.2f) = %.3f -> (%.2f %.2f %.2f) = %.3f\n',i,p(1),p(2),p(3),v,pn(1),pn(2),pn(3),vn);
end
if(fPlot)
figure(niter); hold on
PlotSpherical(fod,points,0.1);   
scatter3(v*p(1)  ,v*p(2)  ,v*p(3)  ,100,'k');
scatter3(vn*pn(1),vn*pn(2),vn*pn(3),100,'k','filled');
end % plot
%
end % while
% === end iterative process === %
% ========= save maxima ========= %   
H(i,1:3) = pn;      % point 
H(i,4)   = fod(mm); % fod value at point
H(i,5)   = mm;      % point index 
% ========= save maxima ========= %
end % i

% ======================================================================= %
% ======================================================================= %
% ======================================================================= %

% (***) Restore offset to peaks 
fod    = fod + fodoff;
H(:,4) = H(:,4) + fodoff;
sfod   = sum(fod);   

% clusterize
[~,~,c] = unique(H(:,5));
nlobes  = max(c); % number of clusters
% sort by volume (sum of fod values for the respective points)
nc      = zeros(nlobes,1);
for i = 1:nlobes
nc(i) = sum( fod( c == i ) )./sfod;        
end
[ncs,ic] = sort(nc,'descend');
cs       = c;
for i = 1:length(ic)
cidx = c == ic(i);
cs(cidx) = i;
end
% save sorted indices
H(:,6)  = cs;
% get direction of peaks
[~,b,~] = unique(cs);
peaks   = zeros(nlobes,3);
fodvals = zeros(nlobes,2);
for i = 1:nlobes
peaks(i,:)   = H(b(i),1:3);   
fodvals(i,1) = ncs(i); 
fodvals(i,2) = H(b(i),4); 
end % i
% % check sorting
% ncs       = histcounts(cs,edges);
% figure; 
% subplot(1,2,1); bar(1:nlobes,nc);
% subplot(1,2,2); bar(1:nlobes,ncs);
end % SegmentFODFast

function [fodLs,plmLs,plLs,vfLs] = getFODLobesOff(fod,infolobes,nmax,lmax,SHBasis)
% get info from lobe segmentation
ilobes   = infolobes(:,6);
[~,id,~] = unique(ilobes);
nlobes   = length(id);
nmax     = min([nmax nlobes]); % max number of lobes
N_l      = SHnorms(lmax);      % plm normalization factors
lF       = length(fod);        % number of points in which fod is sampled
nF       = 4*pi/lF;            % Norm factor: 4*pi*sum(fod)/length(fod)=1 => 4*pi*plm/Nl=>p(0,0)=1;
thr      = 0.65;               % threshold for spurious peaks (good range for l=6 (noisy) 0.55-0.75), 
nsh      = lmax2nsh(lmax);     % num of SH basis
sF       = sum(fod);           % to get volume fraction for each lobe
                               
% ===================
% (***) Remove offset
fodo   = fod;
% fodoff = min(fod);
% fod    = fod - fodoff;
% ===================

vfLs  = zeros(nmax,1);   % volume fraction each lobe
plmLs = zeros(nsh,nmax); % plms each lobe
fodLs = zeros(lF,nmax);  % amplitudes each lobe 

for i = 1:nmax
% make non lobe values equal to zero 
idxL               = ilobes == i;
% FODL               = fod;     % (1) without offset, always >= 0
% FODL(~idxL)        = 0;       % (1) 
FODLo              = fodo;    % (2) with offset, sometimes < 0 
FODLo(~idxL)       = 0;       % (2) 
% FODLnn             = fodo;    % (3) with offset, always >= 0
% FODLnn(~idxL)      = 0;       % (3) 
% FODLnn(FODLnn<0)   = 0;       % (3) 
% FODLth             = fodo;    % (4) with offset, always >= 0.1 
% FODLth(~idxL)      = 0;       % (4) 
% ----------------------------------------------------------------------- %
thr1 = thr*max(FODLo);         
FODLo(FODLo<thr1) = 0;      
FODLnm = FODLo;
% ---
FODLm      = RescaleFOD(FODLnm,SHBasis); % rescale fod, so max(fod) keeps the same after SH
sL         = sum(FODLm);    % lobe volume
vfLs(i)    = sL/sF;         % lobe volume fraction (wrt original fod)
fodLs(:,i) = (FODLm/sL)/nF; % Normalization should be: 4*pi*fod/length(fod) = 1;  
fodfit     = fodLs(:,i) - (1/sqrt(4*pi))*SHBasis(:,1);
plms       = pinv(SHBasis(:,2:end))*fodfit; % HarmonicsFourierTransform(fodfit,SHBasis(:,2:end),1);
plmLs(:,i) = [1/sqrt(4*pi); plms]; % to normalize as: 4*pi*plm/Nl => p(0,0) = 1;
end % i
% Compute pls (normalized)
plLs = plm2pl(4*pi*plmLs./N_l);
end % getFODLobes

function fods = RescaleFOD(fodori,SH)
% interpolate fod using SH
plmint = pinv(SH)*fodori; % HarmonicsFourierTransform(fodori,SH,1); 
fodint = SH*plmint; % HarmonicsFourierTransform(plmint,SH,0);
% rescale interpolated lobe so plm peaks have same amplitude as original
mori   = max(fodori);
mint   = max(fodint);
fods   = mori*fodori/mint;
% % check
% plmchk = HarmonicsFourierTransform(fods,SH,1); 
% fodchk = HarmonicsFourierTransform(plmchk,SH,0);
% mchk   = max(fodchk);
% fprintf('pk1=%f, pk2=%f, pk3=%f\n',mori,mint,mchk);
% lmax   = nsh2lmax(length(plmint));
% N_l    = SHnorms(lmax);
% pln    = plm2pl(4*pi*plmint./N_l);
% pls    = plm2pl(4*pi*plmchk./N_l);
end % RescaleFOD
