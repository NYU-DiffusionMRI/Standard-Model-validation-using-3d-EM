function SaveNII(var, fname, example, bgvalue)
% savenii(var, fname, example)
%
% Saves variable "var", sized [x, y, z, m], as a NiFTi file.
% NiFTi filename "fname" must be provided by the user. A NiFTi
% template must be given as well "example". Typically, this is
% the original diffusion-weigthed data to ensure all maps are
% aligned with the original data. variables are converted to
% singles.

if ~exist('bgvalue', 'var') || isempty(bgvalue)
   bgvalue = 0;
end

var(isnan(var)) = bgvalue;

nii = load_untouch_nii(example);
if(size(var, 4) == 1)
nii.hdr.dime.dim(1)   = 3;
else
nii.hdr.dime.dim(1)   = 4;
end
fov = (nii.hdr.dime.dim(2:4).*nii.hdr.dime.pixdim(2:4));
nii.hdr.dime.dim(2)   = size(var, 1);
nii.hdr.dime.dim(3)   = size(var, 2);
nii.hdr.dime.dim(4)   = size(var, 3);
nii.hdr.dime.dim(5)   = size(var, 4);
nii.hdr.dime.pixdim(2:4) = fov./nii.hdr.dime.dim(2:4);
nii.hdr.dime.datatype = 16;
nii.hdr.dime.bitpix   = 32;
nii.hdr.dime.scl_slope= 1; 
nii.hdr.dime.scl_inter= 0; 
nii.img = single(var);
save_untouch_nii(nii, fname);

end