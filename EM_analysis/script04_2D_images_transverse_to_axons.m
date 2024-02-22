function script04_2D_images_transverse_to_axons()
% get 2D original/segmented images transverse to the axons
%
% This script needs the original volumes and the results from script:
% script02_proofread_vols.m
% It also needs the information in the text file
% 000_animal_exp_axons_cells_nx_ny_nz_cclims_cglims.txt
% The results of this script will be saved in folder
% White_matter_EM_04_2D_Transverse_Sections/*/
% 16-Sep-2022 by Ricardo Coronado-Leija
close all
clc
rng(0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% folders needed
fprintf('Begin\n')
folderqn    = '';
ir          = 2; % 1 = HM, 2 = LM
res         = {'HM','LM'};
foldersegs  = 'White_matter_EM/';
foldervols  = 'White_matter_EM_FullSubstrates/';
fOut        = 'White_matter_EM_04_2D_Transverse_Sections/';
if(~isfolder(fOut)); mkdir(fOut); end
% name of files needed
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
voxsizes    = [15 50]; % now volumes are isotropic
vs          = voxsizes(ir)*1e-3;
% bsbig       = round(15/vs); % pixels
bsbig       = round(25/vs); % pixels
pslice      = 0.25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f = nsamples:-1:1
sample = [num2str(idsample(f)) '_' sidesample{f}];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data
tic
if(ir == 2)
ccdir     = [2 1 1 1 1 1 1 1 1 1]; % z-direction (for cg its always 3)   
dataname  = [foldersegs samples{f} '/' res{ir} '_' sample '.h5']; 
else
ccdir     = [2 1 1 3 1 1 1 2 2 1]; % z-direction       
dataname  = [foldervols samples{f} '/' res{ir} '_' sample '_data.h5']; 
end % name original
data      = squeeze(h5read(dataname,'/raw'));
subsname  = [foldervols samples{f} '/' res{ir} '_' sample '_' volcleaned];
Substrate = load(subsname).SubstrateCleaned;  
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
[nx,ny,nz] = size(Substrate);
propsi(8:13) = [1 nx 1 ny 1 nz];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% divide the volumes into regions
if(fcc) 
CCsubs = Substrate(propsi(8):propsi(9),propsi(10):propsi(11),propsi(12):propsi(13));
CCdata = data     (propsi(8):propsi(9),propsi(10):propsi(11),propsi(12):propsi(13));
end
if(fcg)
CGsubs = Substrate(propsi(14):propsi(15),propsi(16):propsi(17),propsi(18):propsi(19));
CGdata = data     (propsi(14):propsi(15),propsi(16):propsi(17),propsi(18):propsi(19));
end
clear('Substrate','data');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcc)
% ======================================================================= %    
% get middle slice (always transverse to the axons)
% ======================================================================= %
[ccnr,ccnc,ccns] = size(CCdata);
if(ccdir(f) == 1)
islice      = round(pslice*ccnr);    
ccslicedata = squeeze( CCdata(islice,:,:) );
ccslicesubs = squeeze( CCsubs(islice,:,:) );
elseif(ccdir(f) == 2)
islice      = round(pslice*ccnc);
ccslicedata = squeeze( CCdata(:,islice,:) );
ccslicesubs = squeeze( CCsubs(:,islice,:) );
elseif(ccdir(f) == 3)
islice      = round(pslice*ccns);    
ccslicedata = squeeze( CCdata(:,:,islice) );
ccslicesubs = squeeze( CCsubs(:,:,islice) );
else
fprintf('Not valid')
end % get
clear('CCdata','CCsubs')
% ======================================================================= %
% cut images
% ======================================================================= %
[ccnx,ccny] = size(ccslicedata);
ccnx2       = round(ccnx/2);
ccny2       = round(ccny/2);
if(ir == 2)
ccspan = bsbig;
else
ccspan = round( min([ccnx ccny])/2 )-1;
end % cut
ccxidx = max([1 (ccnx2-ccspan)]):min([ccnx (ccnx2+ccspan)]);
ccyidx = max([1 (ccny2-ccspan)]):min([ccny (ccny2+ccspan)]);
ccslicedatacut =       ( ccslicedata(ccxidx,ccyidx) );
ccslicesubscut = double( ccslicesubs(ccxidx,ccyidx) )*(255/2) ;
clear('ccslicedata','ccslicesubs')
% ======================================================================= %
% conf images
% ======================================================================= %
ccR = ccslicesubscut;
ccG = ccslicesubscut;
ccB = ccslicesubscut;
% === new 13-Jun-2023 (color)
ccR(ccslicesubscut == 0)     = 255;
ccG(ccslicesubscut == 0)     = 255;
ccB(ccslicesubscut == 0)     = 255;
ccR(ccslicesubscut == 127.5) = 232;
ccG(ccslicesubscut == 127.5) = 141;
ccB(ccslicesubscut == 127.5) = 169;
ccR(ccslicesubscut == 255)   = 255;
ccG(ccslicesubscut == 255)   = 192;
ccB(ccslicesubscut == 255)   = 0;
% === new (end)
ccR(ccslicesubscut > 255)    = 122; % (64)
ccG(ccslicesubscut > 255)    = 102; % (64)
ccB(ccslicesubscut > 255)    = 70;  % (64)
ccRGB  = cat(3,ccR,ccG,ccB);
ccsz   = size(ccRGB);    
ccdata = imresize (uint8(ccslicedatacut),2,'nearest');
ccRGB2 = imresize3(uint8(ccRGB),[2*ccsz(1) 2*ccsz(2) ccsz(3)],'nearest');
clear('ccR','ccG','ccB')
clear('ccslicedatacut','ccslicesubscut','ccRGB');
% ======================================================================= %
% save images
% ======================================================================= %
imwrite(ccdata,[fOut '/' res{ir} '_' sample '_cc_data.png']);
imwrite(ccRGB2,[fOut '/' res{ir} '_' sample '_cc_seg.png']);
clear('ccdata','ccRGB2');
end % cc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(fcg)
% ======================================================================= %    
% get middle slice (always transverse to the axons)
% ======================================================================= %
[~,~,cgns]  = size(CGdata);  
islice      = round(pslice*cgns); 
cgslicedata = squeeze( CGdata(:,:,islice) );
cgslicesubs = squeeze( CGsubs(:,:,islice) );
clear('CGdata','CGsubs')
% ======================================================================= %
% cut images
% ======================================================================= %
[cgnx,cgny] = size(cgslicedata);
cgnx2       = round(cgnx/2);
cgny2       = round(cgny/2);
if(ir == 2)
cgspan = bsbig;
else
cgspan = round( min([cgnx cgny])/2 )-1;
end % cut
% cgxidx = (cgnx2-cgspan):(cgnx2+cgspan);
% cgyidx = (cgny2-cgspan):(cgny2+cgspan);
cgxidx = max([1 (cgnx2-cgspan)]):min([cgnx (cgnx2+cgspan)]);
cgyidx = max([1 (cgny2-cgspan)]):min([cgny (cgny2+cgspan)]);
cgslicedatacut =       ( cgslicedata(cgxidx,cgyidx) );
cgslicesubscut = double( cgslicesubs(cgxidx,cgyidx) )*(255/2) ;
clear('cgslicedata','cgslicesubs')
% ======================================================================= %
% conf images
% ======================================================================= %
cgR = cgslicesubscut;
cgG = cgslicesubscut;
cgB = cgslicesubscut;
% === new 13-Jun-2023 (color)
cgR(cgslicesubscut == 0)     = 255;
cgG(cgslicesubscut == 0)     = 255;
cgB(cgslicesubscut == 0)     = 255;
cgR(cgslicesubscut == 127.5) = 232;
cgG(cgslicesubscut == 127.5) = 141;
cgB(cgslicesubscut == 127.5) = 169;
cgR(cgslicesubscut == 255)   = 255;
cgG(cgslicesubscut == 255)   = 192;
cgB(cgslicesubscut == 255)   = 0;
% === new (end)
cgR(cgslicesubscut > 255)    = 122; % (64)
cgG(cgslicesubscut > 255)    = 102; % (64)
cgB(cgslicesubscut > 255)    = 70;  % (64)
cgRGB  = cat(3,cgR,cgG,cgB);
cgsz   = size(cgRGB);    
cgdata = imresize (uint8(cgslicedatacut),2,'nearest');
cgRGB2 = imresize3(uint8(cgRGB),[2*cgsz(1) 2*cgsz(2) cgsz(3)],'nearest');
clear('cgR','cgG','cgB')
clear('cgslicedatacut','cgslicesubscut','cgRGB');
% ======================================================================= %
% save images
% ======================================================================= %
imwrite(cgdata,[fOut '/' res{ir} '_' sample '_cg_data.png']);
imwrite(cgRGB2,[fOut '/' res{ir} '_' sample '_cg_seg.png']);
clear('cgdata','cgRGB2');
end % cg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',sample)
end % f
%
fprintf('End\n')
end % main