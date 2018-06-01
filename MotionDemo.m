%% Demo of the motion correction algorithm
%
%
% Disclaimer: This code is just an early version.
%             It is slow, and may crash, as it has not been polished
%
% 
% This temo will generate the needed data to simulate and test the motion
% correction algorithm proposed in article
%
% -------------------------------------------------------------------------
%   A general method for motion compensation in x-ray computed tomography
%
%   Ander Biguri1, Manjit Dosanjh2, Steven Hancock2 and Manuchehr Soleimani1
%   Published 24 July 2017 • © 2017 Institute of Physics and Engineering in Medicine 
%   Physics in Medicine & Biology, Volume 62, Number 16
% -------------------------------------------------------------------------
%
% The code will generate projections, motion maps and will run
% recosntruction. It has description sof the code along the way.
%
% Warning: this demo can take a lot of time to run. Consider running it in
% pieces.
%
%
%
%% Section 0: Set up. 
% This only needs to be ran once, compiles and sets up everything needed.
clear;clc;
% make sure we are in the right folder

% Compile Motion
mex -largeArrayDims ./Motion_Source/AxM.cpp   ./Motion_Source/ray_interpolated_projection_motion.cu ./Motion_Source/Siddon_projection.cu -outdir ./Motion_Correction
mex -largeArrayDims ./Motion_Source/AtbM.cpp ./Motion_Source/voxel_backprojection_motion.cu ./Motion_Source/voxel_backprojection2.cu ./Motion_Source/voxel_backprojection_parallel.cu -outdir ./Motion_Correction
addpath('./Motion_Correction')
% % % Compile Deformation Vector Field inverse
% 
% 
% % Done!

%% Section 1: Generation of needed data 
%
%
% The demo will describe how to get the results from the first experiment
% in the paper descriving this method. 
% In the article, 100 projections are used, with 100 "motion states",
% however, for this demo only 2 will be used. That will help by haing an
% easier example to follow, and a less memory expensive one on that. 
% However, note that the exact same logic followed here can be used for any
% arbitrary amount of "motion states".
%


% This code is a snapshot of the code used to generate the Deformation
% Vector fields (DVFs). It is the same as the ones in the article. For
% simplification, the code to invert the DVF has not been added and instead
% the resul itself is added 

% ------------------------------------------------------------------------
% parameter for the motion maps. They will be 128 voxels wide.
% L=127/2;
% [x,y,z]=meshgrid(0:127,0:127,0:127);
% 
% t=4; 
% % in the paper t=linspace(0,8,100); and the following equation is evaluated
% % for each t.
% 
% 
%     B(:,:,:,1)=-t.*sin(x.*pi/L).*sin(y.*pi/L).*sin(z.*pi/L);
%     B(:,:,:,2)=-t.*sin(x.*pi/L).*sin(y.*pi/L).*sin(z.*pi/L);
%     B(:,:,:,3)=-t.*sin(x.*pi/L).*sin(y.*pi/L).*sin(z.*pi/L);
% % clear; % we do not need it now.
% F=backwards2forwards(single(B));
% F=single(F);
% B=single(B);
% dvf_F.x=squeeze(F(:,:,:,1));dvf_F.y=squeeze(F(:,:,:,2));dvf_F.z=squeeze(F(:,:,:,3));
% dvf_B.x=squeeze(B(:,:,:,1));dvf_B.y=squeeze(B(:,:,:,2));dvf_B.z=squeeze(B(:,:,:,3));
% 
% save(['./Motion_Correction/motionBF0' num2str(1) '.mat'],'dvf_F','dvf_B')
% ------------------------------------------------------------------------


% This file contains 2 DVFs, B and F. One of them is the forward
% devormation and the other the inverse map of that deformation. These are
% needed
load('.\Motion_Correction\motionBF01.mat')

%% Section 2; Generate projections

% This following code is standard TIGRE stuff. We are just creating a
% geometry for the CBCT problem. 

% --------------------------------------------------------------------------------

% define geometry
geo.DSD = 1536;                             % Distance Source Detector      (mm)
geo.DSO = 1000;                             % Distance Source Origin        (mm)
% Detector parameters
geo.nDetector=[512; 512]/2;					% number of pixels              (px)
geo.dDetector=[0.8; 0.8]*2; 					% size of each pixel            (mm)
geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
% Image parameters
geo.nVoxel=[128;128;128];                   % number of voxels              (vx)
geo.sVoxel=[256;256;256];                   % total size of the image       (mm)
geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
% Offsets
geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)
geo.offDetector=[0; 0];                     % Offset of Detector            (mm)
% Auxiliary
geo.accuracy=0.1;

% --------------------------------------------------------------------------------

% Lets generate the data. First, we are going to use 100 uniformly
% distributed angles.

nangles=100;
angles=linspace(0,2*pi-2*pi/nangles,nangles);


% load our reference image
load('./Motion_Correction/img128.mat');
thorax=img; clear img %Just renamig the variable

% Now, we are going to generate samples "moved". Half of the projections we
% generate will be in one of the motion states, the other half in the
% other. In order to link the concepts here with the concepts in the paper in Fig 2,
% lets call the state on where the thorax is "static" (a) and the state on
% where the thorax has moved (with the equation from the previous ection)
% (b). 
%
% This means that our variable "thorax" is now in state (a).
% lets simulate data for that state

projections=Ax(thorax,geo,angles);

% Now, lets warp the image to obtain it after motion, in state (b)

[x,y,z]=meshgrid(1:128,1:128,1:128);
thoraxM=interp3(thorax,    x+dvf_F.x,...
                           y+dvf_F.y,...
                           z+dvf_F.z,    'linear',0);
clear x y z % we dont need them anymore.

% We can see the images in state (a) and state (b) with the following code
plotImg([thorax thoraxM],'dim',3,'clims',[0, 0.02])

%% Small Note
% We have been working on "voxel space" until now, but the rest of the code
% expects the units to be in mm, thus we need to multiply the number of
% voxels moved by the size of the voxel.


dvf_F.x=dvf_F.x*geo.dVoxel(1);
dvf_F.y=dvf_F.y*geo.dVoxel(2);
dvf_F.z=dvf_F.z*geo.dVoxel(3);

dvf_B.x=dvf_B.x*geo.dVoxel(1);
dvf_B.y=dvf_B.y*geo.dVoxel(2);
dvf_B.z=dvf_B.z*geo.dVoxel(3);


%% Section 3: proof of principle

% now lets generate the other half of the projections.
projectionsM=Ax(thoraxM,geo,angles);

% Right. Now we can start linking what we have with what its in Figure 2 in
% the article. Thorax and thoraxM correspond to the blob in image (a) and
% (b) respectively, an static object and a moved one.
%
% In the same exaple, we now have "projections" and "projectionsM". The
% first one corresponds to the blue line in (a) and the second one to the
% red line in (b). 
% 
% lets visualize this.
plotProj([projections projectionsM],angles)

% The main scientific contribution to the method is hidden in this
% diagrams and this lines. The integral over both blue lines yield the
% exact same result, as does the integral over the red lines. One of them
% takes a straigth path and the other doesnt (depending if its in (a) or
% (b)), but both yield the same result. 
% 
% Thus, if we are able to "warp" the way we simulate the X-rays in the
% function Ax(), then we should be able to generate the projections from (b),
% using the thorax from (a), or viceversa. Lets try that with the modified
% Ax, AxM().


projections_moved_b=AxM(thorax,geo,angles,'interpolated',dvf_F);
projections_moved_a=AxM(thoraxM,geo,angles,'interpolated',dvf_B);

% In the following figure you can see
%
% (1) (2)
% (3) (4)
% 
% (1) projections in state (a) using straight lines and the image from (a).
%  Blue line in (a)
% (2) projections in state (b) using straight lines and the image from (b).
%  Red line in (b)
% (3) projections in state (a) using warped lines and the image from (b).
%  Red line in (a)
% (4) projections in state (b) using warped lines and the image from (a).
%  Blue line in (b)
plotProj([projections_moved_a projections_moved_b;projections projectionsM;],angles,'clims',[0 3.5],'colormap','magma','slice',75)

% Here it is the difference. Note how they are almost the same. Some errors
% due to interpolations are normal.
plotProj([projections-projections_moved_a ; projectionsM-projections_moved_b],angles,'colormap','magma','slice',75)

% As a last test, lets "FDK"  the 4 datasets. FDK is not the correct
% algorithm to use with the motion warped backprojection though, as the
% filtering step is only mathematically applicable with the radon
% transform, i.e. with straight lines. However, we can see how we can
% recosntruct the image in another motion state using projections that are
% generated from the image in a completely diferent motion state. 

FDK_state_a_data_a=FDK(projections,geo,angles);
FDK_state_b_data_b=FDK(projectionsM,geo,angles);

% Remember that we have created  "projections_moved_a" just using thoraxM,
% so its been created solely using the image is state (b) 
% See line 192
FDK_state_a_data_b=FDK(projections_moved_a,geo,angles);

FDK_state_b_data_a=FDK(projections_moved_b,geo,angles);


% Top row: reference images
% mid row: normal FDK, generated from projections normaly
% bottom row: FDK using projections generated with motion modelling
plotImg([FDK_state_a_data_b FDK_state_b_data_a;FDK_state_a_data_a FDK_state_b_data_b; thorax thoraxM],'dim',3,'clims',[0, 0.02],'slice',100)
%% Section 4: Using it in an iterative algorithm
% Up to here, the code focuses on showing the principles behind the motion correction
% Lest now run a proper demo with iterative recosntruction algorithms.

% we are, for simplicity, going to assume that the patient did not move for
% half the projections, and then moved and stayed still the rest of it.
% This is far from reality, but proves the fucntionallity. Do change to
% different approaches if you want to see more. In the paper, each of the
% 100 projections has a different motion state, and there is a DVF for each
% of them.

data=cat(3,projections(:,:,1:50), projectionsM(:,:,51:100));
% now data is the data with motion. 
% we are going to recosntruct the image in state (a). Lets use SART as an
% algorithm, but note that any algorithm can be use, just replace Ax and
% Atb by AxM and AtbM and give the correct vector field for each
% projection. 

% Also, note that SART is a quite slow algorithm, and that TIGRE is not
% supposed to be a comercial product now. Please, do not consider the
% current computational times as a reference for using the code. Faster
% algorithms are there, and these algorithms can be certainly accelerated.

% This will be our reference image, without motion. this is the "ideal
% case" to compare against.
% imgSART_nomotion=SART(projections,geo,angles,100);

% This is what SART looks like with no motion compensation and motion in
% the projections
% imgSART_motion=SART(data,geo,angles,100);

%% intermission
% lets clear unwanted data, to save up some memory
clearvars -except geo angles data dvf_B dvf_F imgSART_nomotion imgSART_motion

%%  Reconstruct with motion correction.
%
% Note, this is the code for SART, but modified for motion correction. As
% this is unoptimized and just a demo, there is no "motion corrected"
% algorithms as functions. However, note that the chagnes to the algorithms
% are minor. The only diferent thing is loading of the DVF and the use of
% AtbM() and AxM() instead of Atb() and Ax().

%%  SART motion corrected

% Note that this code is essentially lines 78-190 of SART, with redundant
% code removed (code for optional parameters) and DVF choosing added (each
% projection needs its corresponding DVF).


geoaux=geo;
geoaux.sVoxel(3)=geo.sDetector(2);
geoaux.nVoxel=[2,2,2]'; % accurate enough?
geoaux.dVoxel=geoaux.sVoxel./geoaux.nVoxel;
W=Ax(ones(geoaux.nVoxel','single'),geoaux,angles,'ray-voxel');  %
W(W<min(geo.dVoxel)/4)=Inf;
W=1./W;
[x,y]=meshgrid(geo.sVoxel(1)/2-geo.dVoxel(1)/2+geo.offOrigin(1):-geo.dVoxel(1):-geo.sVoxel(1)/2+geo.dVoxel(1)/2+geo.offOrigin(1),...
    -geo.sVoxel(2)/2+geo.dVoxel(2)/2+geo.offOrigin(2): geo.dVoxel(2): geo.sVoxel(2)/2-geo.dVoxel(2)/2+geo.offOrigin(2));
A = permute(angles+pi/2, [1 3 2]);
V = (geo.DSO ./ (geo.DSO + bsxfun(@times, y, sin(-A)) - bsxfun(@times, x, cos(-A)))).^2;
V=single(V);
% V=single(sum(V,3));

clear A x y dx dz;

res=zeros(geo.nVoxel','single');

% mvf.x=zeros([128,128,300],'single');
% mvf.y=zeros([128,128,300],'single');
% mvf.z=zeros([128,128,300],'single');
lambda=1;
niter=100;
errorL2=[];
% res=thorax;
for ii=1:niter
    if ii==1; tic; end
    for jj=1:length(angles)
 
        % select the correct DVF
        if ~floor(jj/51) % the first 50 projections
            dvf.x=zeros(size(dvf_F.x),'single');
            dvf.y=zeros(size(dvf_F.x),'single');
            dvf.z=zeros(size(dvf_F.x),'single');
        else
            dvf.x=dvf_F.x;
            dvf.y=dvf_F.y;
            dvf.z=dvf_F.z;
        end
        % --------------------

        proj_err=data(:,:,jj)-AxM(res,geo,angles(jj),'interpolated',dvf);
        
        weighted_err=W(:,:,jj).*proj_err;

        % select the correct DVF
        if ~floor(jj/51)
            dvf.x=zeros(size(dvf_F.x),'single');
            dvf.y=zeros(size(dvf_F.x),'single');
            dvf.z=zeros(size(dvf_F.x),'single');
        else
            dvf.x=dvf_B.x;
            dvf.y=dvf_B.y;
            dvf.z=dvf_B.z;
        end
        % --------------------
        backprj=AtbM(weighted_err,geo,angles(jj),'FDK',dvf);
        
        weigth_backprj=bsxfun(@times,1./V(:,:,jj),backprj);
        res=res+lambda*weigth_backprj;
        
        res(res<0)=0;
    end
    
     if (ii==1);
        expected_time=toc*niter;
        disp('SART Motion');
        disp(['Expected duration  :    ',secs2hms(expected_time)]);
        disp(['Exected finish time:    ',datestr(datetime('now')+seconds(expected_time))]);
        disp('');
    end
end
imgSART_corrected=res;

%% clear everything
clearvars -except  imgSART_nomotion imgSART_motion imgSART_corrected geo angles data dvf_B dvf_F

%% section 5: Show result!
%
% This is a similar result that figure 4. 
plotImg([imgSART_nomotion imgSART_motion imgSART_corrected],'Dim',3,'Slice', 100,'clims',[0, 0.02])
%
% Note that with this method
%
% 1) you use all projections, regarthless of their moment in motion, to
% recosntruct static images, as long as you know the motion for each of the
% projections.
% 
% We have recosntructed the state (a) of the image, but changin lines 315 and 331
% to " if floor(jj/51) " and then  change the lines 320-322 and 336-338 to
% use the other vector fields. That will reconstruct the image in state
% (b).
%
% 2) as the modification to the algorithm is the Ax() and Atb() fucntions,
% used for all algorithms, then any other algorithm can be changed to a
% motion compensated one. In the paper we show in figure 8 (d) the TV
% spatial regularization algorithm. (MA)-ROOSTER is essentially reconstructing
% N of these at the same time and then temporally regularizing it. With
% motion information, technically, one can recosntruct the N frames of the
% 4D image each with the full projection dataset, and then regularize
% acordingly. ROOSTER is just an example here.
%
%



