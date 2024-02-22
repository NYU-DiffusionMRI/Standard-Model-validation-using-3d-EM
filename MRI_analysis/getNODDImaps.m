function getNODDImaps(pathIn,pathbval,pathbvec,pathMask,folderOut,Basename,NODDImodel,NODDIfix,NODDIval)
% Fit and get NODDI parameter maps 
%
% pathIn: path to dwi (.nii) file     
% pathbval: path to bvalues (.bval) file (FSL style)
% pathbval: path to bvectors (.bvec) file (FSL style)
% pathMask: path to mask (.nii) file
% folderOut: output directory
% Basename: output basename
% NODDImodel: name of the model to fit with NODDI toolbox, see GetParameterStrings.m file
% NODDIfix: array of 0/1 indicating the parameter to fix, check GetParameterStrings(NODDImodel) 
% NODDIval: array of doubles indicating the values for the fixed parameter
%
% This script calls the NODDI functions:
% http://mig.cs.ucl.ac.uk/index.php?n=Tutorial.NODDImatlab
% It needs tools for reading/writting nifti files
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% By Ricardo Coronado-Leija
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Model Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noddi = MakeModel(NODDImodel);
fprintf('NODDI model: %s\n',noddi.name);
fprintf('Parameters:\n');
for i = 1:(noddi.numParams)
if(NODDIfix(i) > 0.1)
noddi.GS.fixed(i)     = NODDIfix(i);
noddi.GD.fixed(i)     = NODDIfix(i);
noddi.GS.fixedvals(i) = NODDIval(i);
noddi.GD.fixedvals(i) = NODDIval(i);    
end % if    
fprintf('%d %.12f %d %.12f : %s\n',noddi.GS.fixed(i),noddi.GS.fixedvals(i), ...
noddi.GD.fixed(i),noddi.GD.fixedvals(i),noddi.paramsStr{i});
end % i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fMask = 1;

% Create output folder
if(~isfolder(folderOut));  mkdir(folderOut); end

% read scheme/protocol
protocol = FSL2Protocol(pathbval,pathbvec,50);
% protocol = BTAB2Protocol(pathScheme,50);

% NODDI create ROI
% read dwi & mask
% ROIout = CreateROI(dwi,mask);
pathROI = [folderOut '/' Basename '_roi.mat'];
CreateROI(pathIn,pathMask,pathROI);
ROIout = load(pathROI);
% NODDI fit 
% NODDIout = batch_fitting(ROIout,protocol,noddi);
pathFit = [folderOut '/' Basename '_fit.mat'];
batch_fitting(pathROI,protocol,noddi,pathFit);
NODDIout = load(pathFit);

% Correct and save Parameters
ficvfF        = 0;
diF           = 0;
kappaF        = 0;
NODDIout.mlps = removebad(NODDIout.mlps);
for i = 1:NODDIout.model.numParams
vol = 0*ROIout.mask;
for ii = 1:size(ROIout.idx,1)
x = ROIout.idx(ii,1);
y = ROIout.idx(ii,2);
z = ROIout.idx(ii,3);
vol(x,y,z) = NODDIout.mlps(ii,i);
end  % ii

if(strcmp(NODDIout.model.paramsStr{i},'ficvf'))
ficvf  = vol;    
ficvfF = 1;
end
if(strcmp(NODDIout.model.paramsStr{i},'di'))
di  = vol;    
diF = 1;
end
if(strcmp(NODDIout.model.paramsStr{i},'kappa'))
vol    = 10*vol;
kappa  = vol;    
kappaF = 1;
end

SaveNII(vol,[folderOut '/' Basename '_' NODDIout.model.paramsStr{i} '.nii.gz'] , pathMask,[]);

end % i


if(ficvfF && diF)
deperp = (1-ficvf).*di;
SaveNII(deperp,[folderOut '/' Basename '_deperp.nii.gz'] , pathMask,[]);
end

if(kappaF)   
odi    = (2/pi)*atan(1./kappa); 
p2     = 1/4*(3./(sqrt(kappa).*dawson(sqrt(kappa)))-2-3./kappa);
p2(isnan(p2)) = median(p2(:),'omitnan');
p4     = (1./(32*kappa.*kappa)).*( 105 + 12*kappa.*(5+kappa) + ( 5*sqrt(kappa).*(2*kappa-21) )./( dawson(sqrt(kappa)) ) );
p4(isnan(p4)) = median(p4(:),'omitnan');
adisp  = acosd(sqrt(2*p2/3+1/3));
SaveNII(odi  ,[folderOut '/' Basename '_odi.nii.gz']   , pathMask,[]);
SaveNII(p2   ,[folderOut '/' Basename '_p2.nii.gz']    , pathMask,[]);
SaveNII(p4   ,[folderOut '/' Basename '_p4.nii.gz']    , pathMask,[]);
SaveNII(adisp,[folderOut '/' Basename '_adisp.nii.gz'] , pathMask,[]);
end


end % main

function volout = removebad(volin,val)
if(nargin == 1)
val = 0;    
end
volout = volin;
volout(isnan(volin) | isinf(volin)) = val;
end