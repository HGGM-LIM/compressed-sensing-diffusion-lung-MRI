function Demo_SIDER
% function Demo_SIDER.m
%
% The proposed method incorporates the knowledge of the signal decay into
% the reconstruction (SIDER) to accelerate the acquisition of MR diffusion
% data by undersampling in both spatial and b-value dimensions.
%
% Demo for reconstructing diffusion lung MR images with spatial total
% variation (S-TV) and spatiotemporal total variation (ST-TV) using the
% Split Bregman formulation. SIDER combines total variation (TV) with a
% penalty function that promotes sparsity across the b-direction as follows
% min_u beta|grad_x,y u|_1 + gamma|M u|_1 st. ||Fu-f||^2 < delta, where the
% first term corresponds to spatial TV and M is is an operator that encodes
% the relationship between ventilation images for consecutives values of b,
% based on a stretched exponential model.
%
% The demo loads fully sampled data, undersampled data and reconstruct
% images using both TV and SIDER methods.
% Data is undersampled using a modified version of Lustig's variable
% density pdf, downloaded from (SparseMRI V0.2)
% http://web.stanford.edu/~mlustig/SparseMRI.html, see M. Lustig, D.L
% Donoho and J.M Pauly "Sparse MRI: The Application of Compressed Sensing
% for Rapid MR Imaging" Magnetic Resonance in Medicine, 2007 Dec;
% 58(6):1182-1195.
%
% Code downloaded from repository:
% https://github.com/HGGM-LIM/compressed-sensing-diffusion-lung-MRI
%
% If you use this code, please, cite the following paper:
% JFPJ Abascal, M Desco, J Parra-Robles. Incorporation of prior knowledge
% of the signal behavior into the reconstruction to accelerate the
% acquisition of MR diffusion data,(submitted for publication) 2016.
%
% Departamento de Bioingeniería e Ingeniería Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% juanabascal78@gmail.com, jmparra@hggm.es, desco@hggm.es

% -------------------------------------------------------------------------
% Load complex image
load('DataControl','b','uTarget');

% Choose one slice
zSlice      = 4;

uTarget     = uTarget(:,:,:,zSlice); 
b           = b(1:5);
N           = size(uTarget);

[nel_Lm,xcenters_Lm]    = hist(150:20:450,30);

% -------------------------------------------------------------------------
% MASK FOR THE LUNG to define the region where to estimate D, alpha and Lm
% maps (this estimation is a pixel-by-pixel fitting so its calculation is
% not affected by the mask). The mask is automatically computed by
% thresholding the ventilation image
se      = strel('disk',3);
im      = imclose(abs(uTarget(:,:,1)),se);   % imdilate
[val, ind] = sort(im(:),'descend');
im      = im/val(5);
mask    = im>=0.3;                  % thresholding
se      = strel('disk',1);          % remove small objects
mask    = (imopen(abs(mask),se))>0;

h=figure;
subplot(2,1,1);
imagesc(abs((uTarget(:,:,1))));colorbar; colormap hot; axis image;
title('b=0');
subplot(2,1,2);
imagesc(abs((mask(:,:,1))));colorbar; colormap hot; axis image;
title('Mask');
pause(2);
close(h);

% -------------------------------------------------------------------------
% ESTIMATE D, alpha and Lm MAPS
%
% Ventilation images are filtered before the estimation (as images are very
% noisy especially for large values of b)
k       = fspecial('gaussian',3,1); % For patients (noisier)
h=figure;
for ip = 1:N(3)
    uTargetSmooth(:,:,ip) = imfilter(uTarget(:,:,ip),k);
    subplot(2,1,1); imagesc(abs(uTarget(:,:,ip)));colorbar; colormap hot; axis image
    title(['Original ventilation image, b ' num2str(ip)]);
    subplot(2,1,2); imagesc(abs(uTargetSmooth(:,:,ip)));colorbar; colormap hot; axis image;
    title('Filtered ventilation image');
    pause(1);
end
close(h);

% Estimate D and alpha using a STRETCHED EXPONENTIAL MODEL
options     = optimset('display','iter','Algorithm','levenberg-marquardt','Jacobian','on');
x0_D        = 0.2*ones(N(1:2));
x0_D        = x0_D(mask);
x0_a        = 0.8*ones(N(1:2));
x0_a        = x0_a(mask);
x0          = [x0_D; x0_a];
x           = lsqnonlin(@objVecPotExpJac,x0,[],[],options,...
    abs(uTargetSmooth(:,:,1)),abs(uTargetSmooth(:,:,2:end)),b(2:end),mask);

Drec0        = zeros(N(1:2));
Drec0(mask)  = abs(x(1:nnz(mask)));
alpharec0    = zeros(N(1:2));
alpharec0(mask)= abs(x(nnz(mask)+1:2*nnz(mask)));

% Apply a threshold to the estimated values
Drec0(abs(Drec0)>0.9)=0; Drec0(abs(Drec0)<0.05)=0;
alpharec0(abs(alpharec0)>1.3)=0; alpharec0(abs(alpharec0)<0.3)=0;

% ESTIMATE Lm map
tmp         = alpharec0;
tmp(logical(mask.*alpharec0>1)) = 0.9; % It doesnt work for values of alpha larger than 1
Lm0         = ComputeLmImage(N,Drec0,tmp);

h=figure;
subplot(2,2,1); imagesc(abs(Drec0)); colormap hot; axis image; colorbar; title('D fully sampled');caxis([0 0.9]);
subplot(2,2,2); imagesc(abs(alpharec0)); colormap hot; axis image; colorbar; title('alpha');caxis([0 1.3]);
subplot(2,2,3); imagesc(abs(Lm0)); colormap hot; axis image; colorbar; title('Lm');caxis([150 450]);
pause(2);
close(h);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMPRESSED SENSING
%
% Pseudo-random undersampling pattern for all values of b
rand('state',1);

% Select parameters for Lustig's variable density pdf: radius, sparsity and
% P. Below we provide the set of paired sparsity and P values used in the
% paper:
% sparsity [0.5 0.25 0.20 0.15 0.1 0.07]
% P        [4   4    4    6    9   13]
radius      = 0.017;
DN          = [N(1),1];
sparsity    = 0.15;      % Corresponds to acceleration factor x7
P           = 6;

RAll    = zeros(N);
h=figure;
for ip = 1:N(3)
    % Loop across of values of b
    pdf     = genPDF(DN,P,sparsity,2,radius,0);
    temp    = genSampling_LIM(pdf,20,1);
    indR    = temp(:,1)==1;
    R       = zeros(N(1:2));
    R(indR,:)= 1;
    R(1,:)= 1;
    spy(R); title(['Undersampling b ' num2str(ip)]);
    RAll(:,:,ip) = R;
    data    = fft2(uTarget(:,:,ip)).*R;
    dataAll(:,:,ip) = data;
    pause(1);
end % ip
close(h);

% -------------------------------------------------------------------------
% ZERO-FILLING RECONSTRUCTION
u_fft       = ifft2(dataAll);

% Smooth ventilation images
for ip = 1:N(3)
    u_fft_s(:,:,ip) = imfilter(u_fft(:,:,ip),k);
end

% POTENTIAL MODEL
options     = optimset('display','iter','Algorithm','levenberg-marquardt','Jacobian','on');
x0_D        = 0.2*ones(N(1:2));
x0_D        = x0_D(mask);
x0_a        = 0.8*ones(N(1:2));
x0_a        = x0_a(mask);
x0          = [x0_D; x0_a];
x           = lsqnonlin(@objVecPotExpJac,x0,[],[],options,abs(u_fft_s(:,:,1)),abs(u_fft_s(:,:,2:end)),b(2:end),mask);

% D and A
Drecfft      = zeros(N(1:2));
Drecfft(mask)= abs(x(1:nnz(mask)));
alpharecfft  = zeros(N(1:2));
alpharecfft(mask)= abs(x(nnz(mask)+1:2*nnz(mask)));
Drecfft(abs(Drecfft)>0.9)=0; Drecfft(abs(Drecfft)<0.05)=0;
alpharecfft(abs(alpharecfft)>1.3)=0; alpharecfft(abs(alpharecfft)<0.3)=0;

% Lm
tmp         = alpharecfft;
tmp(logical(mask.*alpharecfft>1)) = 0.9;%!!!!1 Othewise, Lms doesnt work!!!
Lmfft       = ComputeLmImage(N,Drecfft,tmp);
Lmfft(isnan(Lmfft)) = 0;

h=figure;
subplot(2,2,1); imagesc(abs(Drecfft)); colormap hot; axis image; colorbar; title('D IFFT');caxis([0 0.9]);
subplot(2,2,2); imagesc(abs(alpharecfft)); colormap hot; axis image; colorbar; title('alpha');caxis([0 1.3]);
subplot(2,2,3); imagesc(abs(Lmfft)); colormap hot; axis image; colorbar; title('Lm');caxis([150 450]);
pause(2);
close(h);
% -------------------------------------------------------------------------
% STATIC Spatial Total Variation reconstruction using Split Bregman
% Code download from
% http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html
%
% Goldstein's spatial TV using the Split Bregman formulation
% u_tv = mrics(RAll(:,:,1),data(:,:,1), mu, lambda, gamma, nInner, nBreg);
%
% SpatialTVSB.m: same as mrics.m but it computes the solution error
% norm
mu          = 1;
lambda      = 1;
gamma       = 0.01;
nInner      = 1;
nBreg       = 800;
% Loop across all values of b to reconstruct ventilation images
for ip = 1:N(3)
    % The figure displayed shows (at several iteration numbers) the
    % reconstructed image (u), the dummy variable x that outputs the
    % shrinkage (for the gradient along x-direction) indicating the
    % sparsity as the number of nonzero coefficients, and the solution
    % error norm for all ventilation images
    [uThis,errThis] = SpatialTVSB(RAll(:,:,ip),dataAll(:,:,ip), mu, lambda, gamma, nInner, nBreg,uTargetSmooth(:,:,ip),mask);
    u_tv(:,:,ip) = uThis;
    err_tv(:,ip) = errThis;
end % ip

% Smooth ventilation images
for ip = 1:N(3)
    u_tv_s(:,:,ip) = imfilter(u_tv(:,:,ip),k);
end

% POTENTIAL MODEL
options     = optimset('display','iter','Algorithm','levenberg-marquardt','Jacobian','on');
x0_D        = 0.2*ones(N(1:2));
x0_D        = x0_D(mask);
x0_a        = 0.8*ones(N(1:2));
x0_a        = x0_a(mask);
x0          = [x0_D; x0_a];
x           = lsqnonlin(@objVecPotExpJac,x0,[],[],options,abs(u_tv_s(:,:,1)),abs(u_tv_s(:,:,2:end)),b(2:end),mask);

% D and A
Drectv      = zeros(N(1:2));
Drectv(mask)= abs(x(1:nnz(mask)));
alpharectv  = zeros(N(1:2));
alpharectv(mask)= abs(x(nnz(mask)+1:2*nnz(mask)));
Drectv(abs(Drectv)>0.9)=0; Drectv(abs(Drectv)<0.05)=0;
alpharectv(abs(alpharectv)>1.3)=0; alpharectv(abs(alpharectv)<0.3)=0;

% Lm
tmp         = alpharectv;
tmp(logical(mask.*alpharectv>1)) = 0.9;
Lmtv        = ComputeLmImage(N,Drectv,tmp);
Lmtv(isnan(Lmtv)) = 0;

h=figure;
subplot(2,2,1); imagesc(abs(Drectv)); colormap hot; axis image; colorbar; title('D TV');caxis([0 0.9]);
subplot(2,2,2); imagesc(abs(alpharectv)); colormap hot; axis image; colorbar; title('alpha');caxis([0 1.3]);
subplot(2,2,3); imagesc(abs(Lmtv)); colormap hot; axis image; colorbar; title('Lm');caxis([150 450]);
pause(2);
close(h);
% -------------------------------------------------------------------------
% SIDER
%
% Get mean estimated values of D and alpha for the model
Dest        = mean(Drectv(Drectv>0));
Alphaest    = mean(alpharectv(alpharectv>0));

mu      = 1;
lambda  = 1;
gamma   = 1e-2;
nBreg   = 800;
alpha   = 0.2;
beta    = 0.2;
% The figure displayed shows (at several iteration numbers) the
% reconstructed image (u), the dummy variable x that outputs the
% shrinkage (for the gradient along x-direction) and the dummy variable p
% (for the operator along b-direction), and the solution error norm for all
% ventilation images
[u_t,err_t] = SIDER(Dest,Alphaest,b,dataAll,RAll,N(1:3),mu,lambda,gamma,alpha,beta,nBreg,mask,uTargetSmooth);

% Smooth ventilation images
for ip = 1:N(3)
    u_t_s(:,:,ip) = imfilter(u_t(:,:,ip),k);
end

% POTENTIAL MODEL
options     = optimset('display','iter','Algorithm','levenberg-marquardt','Jacobian','on','TolFun',1e-06);
x0_D        = Dest*ones(N(1:2));
x0_D        = x0_D(mask);
x0_a        = 0.8*ones(N(1:2));
x0_a        = x0_a(mask);
x0          = [x0_D; x0_a];
x           = lsqnonlin(@objVecPotExpJac,x0,[],[],options,abs(u_t_s(:,:,1)),abs(u_t_s(:,:,2:end)),b(2:end),mask);

% D and A
Drect       = zeros(N(1:2));
Drect(mask)= abs(x(1:nnz(mask)));
alpharect   = zeros(N(1:2));
alpharect(mask)= abs(x(nnz(mask)+1:2*nnz(mask)));
Drect(abs(Drect)>0.9)=0; Drect(abs(Drect)<0.05)=0;
alpharect(abs(alpharect)>1.3)=0; alpharect(abs(alpharect)<0.3)=0;

% Lm
tmp         = alpharect;
tmp(logical(mask.*alpharect>1)) = 0.9;
Lmt         = ComputeLmImage(N,Drect,tmp);
Lmt(isnan(Lmt)) = 0;

h=figure;
subplot(2,2,1); imagesc(abs(Drect)); colormap hot; axis image; colorbar; title('D SIDER'); caxis([0 0.9]);
subplot(2,2,2); imagesc(abs(alpharect)); colormap hot; axis image; colorbar; title('alpha'); caxis([0 1.3]);
subplot(2,2,3); imagesc(abs(Lmt)); colormap hot; axis image; colorbar; title('Lm'); caxis([150 450]);
pause(2);
close(h);
% --------------------------------------------------------------------
% Plot solution error norm
figure; plot(mean(err_tv,2)); hold on; plot(mean(err_t,2),'r--');
xlabel('Iteration number'); ylabel('Solution error norm'); legend('TV','SIDER');
title('Solution error vs. iteration number (mean of ventilation images)');

% Ventilation images
figure;
subplot(3,2,1);
imagesc(abs(uTargetSmooth(:,:,1))); colormap hot; axis image; axis off; colorbar; title('u(b=0) fully sampled');
subplot(3,2,2);
imagesc(abs(uTargetSmooth(:,:,5))); colormap hot; axis image; axis off; colorbar; title('u(b=6.4) fully sampled');
subplot(3,2,3);
imagesc(abs(u_tv(:,:,1))); colormap hot; axis image; axis off; colorbar; title('u(b=0) TV');
subplot(3,2,4);
imagesc(abs(u_tv(:,:,5))); colormap hot; axis image; axis off; colorbar; title('u(b=6.4) TV');
subplot(3,2,5);
imagesc(abs(u_t(:,:,1))); colormap hot; axis image; axis off; colorbar; title('u(b=0) SIDER');
subplot(3,2,6);
imagesc(abs(u_t(:,:,5))); colormap hot; axis image; axis off; colorbar; title('u(b=6.4) SIDER');

figure;
subplot(3,3,1); imagesc(abs(Drec0)); colormap hot; axis image; axis off; colorbar; title('D fully sampled'); caxis([0 0.9]);
subplot(3,3,2); imagesc(abs(alpharec0)); colormap hot; axis image; axis off; colorbar; title('\alpha fully sampled'); caxis([0 1.3]);
subplot(3,3,3); imagesc(abs(Lm0)); colormap hot; axis image; axis off; colorbar; title('Lm fully sampled'); caxis([150 450]);
subplot(3,3,4); imagesc(abs(Drectv)); colormap hot; axis image; axis off; colorbar; title('D TV'); caxis([0 0.9]);
subplot(3,3,5); imagesc(abs(alpharectv)); colormap hot; axis image; axis off; colorbar; title('\alpha TV'); caxis([0 1.3]);
subplot(3,3,6); imagesc(abs(Lmtv)); colormap hot; axis image; axis off; colorbar; title('Lm TV');caxis([150 450]);
subplot(3,3,7); imagesc(abs(Drect)); colormap hot; axis image; axis off; colorbar; title('D SIDER');caxis([0 0.9]);
subplot(3,3,8); imagesc(abs(alpharect)); colormap hot; axis image; axis off; colorbar; title('\alpha SIDER');caxis([0 1.3]);
subplot(3,3,9); imagesc(abs(Lmt)); colormap hot; axis image; axis off; colorbar; title('Lm SIDER');caxis([150 450]);

pause(1);

end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [obj,JacThis] = objVecPotExpJac(x,A0,A,b,mask)
% Fit vntilation images to a stretch exponential model to get maps
% of diffusion (D) and heterogeneity index (alpha)
obj     = objVecPotExp(x,A0,A,b,mask);

JacThis = objVecPotJac(x,A0,A,b,mask);
%     end
end

function obj = objVecPotExp(x,A0,A,b,mask)
% Cost function for nonlinear LS
n           = size(A);
nm          = nnz(mask);
Dval        = x(1:nm);
alphaval    = x(nm+1:2*nm);
DThis       = zeros(n(1:2));
DThis(mask) = Dval;
alphaThis   = zeros(n(1:2));
alphaThis(mask) = alphaval;
obj         = [];
for iq = 1:size(A,3)
    tmp     = A0.*exp( -( (b(iq)*DThis).^alphaThis) ) - A(:,:,iq);
    obj     = [obj; (tmp(mask))];
end
end

function JacThis = objVecPotJac(x,A0,A,b,mask)
% Jacobian for nonlinear LS
nm          = nnz(mask);
Dval        = x(1:nm);
alphaval    = x(nm+1:2*nm);
JacThis     = [];
for ig  = 1:size(A,3)
    % Derivative wrt D
    A0e         = A0(mask).*exp(-((b(ig)*Dval).^alphaval));
    ba          = b(ig).^alphaval;
    Da          = Dval.^(alphaval-1);
    JacThisD    = -alphaval.*ba.*Da.*A0e;
    
    % Derivative wrt a
    logbD       = log(b(ig)*Dval);
    bDa         = (b(ig)*Dval).^alphaval;
    JacThisA    = -bDa.*logbD.*A0e;
    
    tmp         = [spdiags(JacThisD,0,nm,nm) spdiags(JacThisA,0,nm,nm)];
    JacThis     = [JacThis; tmp];
end
end

function [Lm,Lstd] = ComputeLmImage(N,Drec0,alpharec0)
% Compute maps of mean alveolar length (Lm) given maps of diffusion
% (D) and heterogeneity index (alpha)
%
Lm      = zeros(N(1:2));
Lstd    = zeros(N(1:2));
indThis = find( (Drec0 > 0).*(alpharec0 > 0) );
for ia = 1:length(indThis)
    [Lm_val,Lstd_val] = Lmdiststretchedexp(Drec0(indThis(ia)),alpharec0(indThis(ia)));
    Lm(indThis(ia))      = Lm_val;
    Lstd(indThis(ia))    = Lstd_val;
end % ia
end

function [Lm,Lmstd,Hk,ldd]=Lmdiststretchedexp(DDC,alpha,diftime, Do)
if nargin==3,
    Do=0.9;
end
if nargin==2,
    Do=0.88;
    diftime=1.6e-3;
end

D=0.01:0.001:Do;
beta=alpha;

b=0:0.5:10;

t0=1./DDC;
k=D;

betavect=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
betai=0.01:0.01:1;
B=[0.145 0.197 0.243 0.285 0.382 0.306 0.360 0.435 0.7];
C=[0.89 0.50 0.35 0.25 0 0.13 0.22 0.4015 0.33];

Bi=interp1(betavect,B,betai,'cubic');
Ci=interp1(betavect,C,betai,'cubic');

delta=beta.*(beta-0.5)./(1-beta);
fk=1+Ci(round(beta.*100)).*(k.*t0).^delta;

Hk=(t0.*Bi(round(beta.*100))./(k*t0).^((1-beta/2)/(1-beta))).*exp(-((1-beta).*beta.^(beta./(1-beta)))./(k.*t0).^(beta./(1-beta))).*fk;
Hk=Hk/sum(Hk);
ldd=1e4*sqrt(2*D*diftime);
Lm=sum(Hk.*ldd);
Lmstd=sqrt(sum(((ldd-Lm).^2).*Hk));
end

%