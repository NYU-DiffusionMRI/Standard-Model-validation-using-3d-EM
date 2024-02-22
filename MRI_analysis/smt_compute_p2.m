function smt_compute_p2(pathSMT,pathS2,pathScheme,folderOut)
% compute FOD rotational invariants from kernel obtained with SMT ( fit with l=0 )
% Sl(b) = pl * Kl(b,smt)
% pl    = "Sl(b)./Kl(b)"
%
% pathSMT     = path to SMT results file (.nii)
% pathS2      = path to Signals Second Order Rotational Invariant (.nii)
%               this is obtained using SMI
% pathScheme  = path to protocol (mrtrix style)
% folderOut   = folder to save output files (usually same directory as pathSMT)
% By Ricardo Coronado-Leija
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~isfolder(folderOut)); mkdir(folderOut); end

% read SMT results
SMTnii    = load_untouch_nii(pathSMT);
smtf      = SMTnii.img(:,:,:,1);
smtda     = SMTnii.img(:,:,:,2)*1e3; 
smtdepar  = smtda;
smtdeperp = SMTnii.img(:,:,:,3)*1e3;
mask      = smtf > 0;

% SMI rotational invariants
S2nii     = niftiread(pathS2);

% protocol (1000,2000,3000)
scheme     = readmatrix(pathScheme);
% bvec       = double(scheme(:,1:3));
b          = double(scheme(:,4))/1e3;
beta       = ones(size(b));
TE         = [];
[bb,~,~]   = SMI.Group_dwi_in_shells_b_beta_TE(b,beta,TE,0.2);
nbs        = size(bb,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
fmask      = vectorize_mask(smtf,mask);
damask     = vectorize_mask(smtda,mask);
deparmask  = vectorize_mask(smtdepar,mask);
deperpmask = vectorize_mask(smtdeperp,mask);
Nvoxels    = length(fmask);

% kernel rotational invariants (uses SMI toolbox)
Dfw    = 2; % irrelevant (SMT does not consider iso) 
s2mask = vectorize_mask(S2nii,mask)';
K2_all = zeros(Nvoxels,nbs);
for ii = 1:Nvoxels
f1           =  fmask(ii);
Da1          =  damask(ii);
Depar1       =  deparmask(ii);
Deperp1      =  deperpmask(ii);
current_x    = [f1,Da1,Depar1,Deperp1,0,1,1];
K2_all(ii,:) = SMI.RotInv_Kell_wFW_b_beta_TE_numerical(2,bb(1,:),bb(2,:),bb(4,:),current_x,Dfw);
end % ii

% remove b0 and compute dot products
k2mask2 = K2_all(:,2:end).*K2_all(:,2:end);
K2_s2_  = abs(K2_all(:,2:end).*s2mask(:,2:end));

% compute p2
p2_all = zeros(Nvoxels,1,'single');
for ii = 1:Nvoxels
p2_all(ii,1) = sum( K2_s2_(ii,:) ) / sum( k2mask2(ii,:) );
end % ii

% save
p2_vol           = vectorize_mask(p2_all',mask);
adisp_vol        = acosd(sqrt(2*p2_vol/3+1/3)); % dispersion angle
p2name           = [folderOut '/smt_p2.nii.gz'];
adispname        = [folderOut '/smt_p2_adisp.nii.gz'];
SaveNII(p2_vol   ,p2name   ,pathSMT);
SaveNII(adisp_vol,adispname,pathSMT);

end % main

