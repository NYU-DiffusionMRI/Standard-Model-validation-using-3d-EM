function getRobustDKImaps(pathIn,pathScheme,pathMask,folderOut,Basename,bthresh,flagB0,fWMTI)
% Fit and get DKI parameter maps using the RobustDKI method
%
% pathIn: path to dwi (.nii) file     
% pathScheme: path to protocol (.txt) file (MRtrix style)
% pathMask: path to mask (.nii) file
% folderOut: output directory
% Basename: output basename
% bthresh: maximum b-value threshold
% flagB0: (flag 1/0) exclude b0s from fitting 
% fWMTI: (flag 1/0) compute WMTI parameters 
%
% This script calls the RobustDKI function:
% https://github.com/jelleveraart/RobustDKIFitting
% and compute DKI parameter maps using the utils functions from DESIGNER 
% https://github.com/NYU-DiffusionMRI/DESIGNER/tree/master/utils
% It needs tools for reading/writting nifti files
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% By Ricardo Coronado-Leija
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMask = 1; 
if(~isfolder(folderOut));  mkdir(folderOut); end

% read 4D volume (nii)
diw_nii      = load_untouch_nii(pathIn); 
dwi          = diw_nii.img;
scheme       = readmatrix(pathScheme);
[nr,nc,ns,~] = size(dwi);
% use or not the mask
if(fMask && isfile(pathMask))
% read mask
niimask      = load_untouch_nii(pathMask);    
mask         = logical(niimask.img);
else
pathMask     = pathIn;       
mask         = true(nr,nc,ns); 
end
% change scheme
btab         = scheme;
% btab, bthresh should be 500,1000,1500,2000 etc ..
btab(:,4)    = btab(:,4)/1000;
bthresh      = bthresh/1000;   
% remove high shells
if(bthresh > 0)
fprintf('Remove big bs\n')          
idxb = btab(:,4) < bthresh; 
else
idxb = true(size(btab,1),1);    
end

% do not use b0 values in the fitting
if(flagB0)
fprintf('Remove b0s\n')     
idxb0 = btab(:,4) > 0.05; 
idxb  = idxb & idxb0; 
end

% Remove negative/zero values
fprintf('Correcting\n')
mindw = min(dwi(dwi > 0));
dwi(dwi <= 0) = mindw;
dwi   = removebad(dwi,mindw);
fprintf('voxels %d\n',sum(mask(:)));

% Fit DKI
fprintf('DKI fit\n')
tic
[dt, b0, ~, ~] = RobustDKIFitting(dwi(:,:,:,idxb), btab(idxb,:), mask);
toc
dt = removebad(dt);
% Get DKI Parameters
tic
[fa,md,rd,ad,fe,mk,rk,ak] = dki_parameters(dt,mask);
toc

% Remove bad voxels
b0  = removebad(b0);
ad  = removebad(ad);
rd  = removebad(rd);
md  = removebad(md);
fa  = removebad(fa);
ak  = removebad(ak);
rk  = removebad(rk);
mk  = removebad(mk);
fe  = removebad(fe);
b0  = removebad(b0);
dti = nyu2mrtrix(dt);
dti = removebad(dti);

% SaveNII(rmse,[folderOut '/' Basename '_rmse.nii.gz'], pathMask,[]);
SaveNII(dt  ,[folderOut '/' Basename '_dt.nii.gz']  , pathMask,[]);
SaveNII(b0  ,[folderOut '/' Basename '_b0.nii.gz']  , pathMask,[]);
SaveNII(ad  ,[folderOut '/' Basename '_ad.nii.gz']  , pathMask,[]);
SaveNII(rd  ,[folderOut '/' Basename '_rd.nii.gz']  , pathMask,[]);
SaveNII(md  ,[folderOut '/' Basename '_md.nii.gz']  , pathMask,[]);
SaveNII(fa  ,[folderOut '/' Basename '_fa.nii.gz']  , pathMask,[]);
SaveNII(ak  ,[folderOut '/' Basename '_ak.nii.gz']  , pathMask,[]);
SaveNII(rk  ,[folderOut '/' Basename '_rk.nii.gz']  , pathMask,[]);
SaveNII(mk  ,[folderOut '/' Basename '_mk.nii.gz']  , pathMask,[]);
SaveNII(fe  ,[folderOut '/' Basename '_pdd.nii.gz'] , pathMask,[]);
SaveNII(dti ,[folderOut '/' Basename '_tensor.nii.gz'] , pathMask,[]);

% Get WMTI Parameters
if(fWMTI)
tic
[awf1,eas1,ias1] = wmti_parameters(dt,mask); % branch Da < De_par
[nr,nc,ns]       = size(awf1);
p2_1             = zeros(nr,nc,ns);
% compute p2
for r = 1:nr
for c = 1:nc  
for s = 1:ns      
    % branch Da < De_par
    Daij          = diag([ias1.da1(r,c,s) ias1.da2(r,c,s) ias1.da3(r,c,s)]);
    Daij_stf      = Daij-diag(ias1.Da(r,c,s)/3*ones(1,3)); % only subtract the diagonal line
    Daij_stf_prod = transpose(Daij_stf(:))*Daij_stf(:);
    p2_1(r,c,s)   = sqrt((3/2)*Daij_stf_prod/(ias1.Da(r,c,s)^2));
end % r
end % c
end % s
% branch Da < De_par
awf1         = removebad(awf1);
p2_1         = removebad(p2_1);
ias1.Da      = removebad(ias1.Da);
eas1.de1     = removebad(eas1.de1);
eas1.de_perp = removebad(eas1.de_perp);
adp2b1       = acosd(sqrt(2*p2_1/3+1/3)); % dispersion angle
toc
% branch Da < De_par
SaveNII(awf1        ,[folderOut '/' Basename '_wmti_f.nii.gz']       , pathMask,[]);
SaveNII(p2_1        ,[folderOut '/' Basename '_wmti_p2.nii.gz']      , pathMask,[]);
SaveNII(adp2b1      ,[folderOut '/' Basename '_wmti_p2_adisp.nii.gz'], pathMask,[]);
SaveNII(ias1.Da     ,[folderOut '/' Basename '_wmti_Da.nii.gz']      , pathMask,[]);
SaveNII(eas1.de1    ,[folderOut '/' Basename '_wmti_Depar.nii.gz']   , pathMask,[]);
SaveNII(eas1.de_perp,[folderOut '/' Basename '_wmti_Deperp.nii.gz']  , pathMask,[]);
end

end % main

function volout = removebad(volin,val)
if(nargin == 1)
val = 0;    
end
volout = volin;
volout(isnan(volin) | isinf(volin)) = val;
end