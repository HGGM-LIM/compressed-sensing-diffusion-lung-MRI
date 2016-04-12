% [u] = SpatialTVSB(R,f, mu, lambda, gamma, nInner,
% nBreg)
% [u,errAll] = SpatialTVSB(R,f, mu, lambda, gamma, nInner,
% nBreg,imTrue)
%
% Inputs:
%
% R         = undersampling matrix, same size as f
% f         = 2D data, which corresponds to fft2(imTrue)+noise
% dimIm     = Nx*Ny*frames
% mu        = parameter weighting the data fidelity term, use mu=1
% lambda    = parameter weighting the TV constraints, use lambda=1
% gamma     = parameter to improve the conditioning, use between mu/100 and
% mu
% nInner    = inner iterations, use n=1
% nBreg     = number of (outer) iterations
% imTrue    = target image to compute the error at each iteration
% 
% Outputs: 
%
% u         = reconstructed image, size dimIm
% errAll    = Relative solution error norm at each iteration for all frames
%
%
% This code is a modified version (to compute the error at each iteration)
% of spatial TV Goldstein'n code mrics.m downloaded from  
% (http://www.ece.rice.edu/~tag7/Tom_Goldstein/Split_Bregman.html), see Tom
% Goldstein and Stanley Osher. The Split Bregman Method for L1-Regularized
% Problems. SIAM J. Imaging Sci., 2(2), 323�343.  
%
% Juan Felipe P�rez-Juste Abascal, Paula Montesinos
% Departamento de Bioingenier�a e Ingenier�a Aeroespacial
% Universidad Carlos III de Madrid, Madrid, Spain
% paumsdv@gmail.com, juanabascal78@gmail.com, desco@hggm.es

function [u,errAll] = SpatialTVSB(R,f, mu, lambda, gamma, nInner, nBreg,varargin)

[rows,cols] = size(f);

% normalize the data so that standard parameter values work
normFactor  = getNormalizationFactor(R,f);
f           = normFactor*f;

% Reserve memory for the auxillary variables
f0          = f;
u           = zeros(rows,cols);
x           = zeros(rows,cols);
y           = zeros(rows,cols);
bx          = zeros(rows,cols);
by          = zeros(rows,cols);

switch nargin
    case 8        
        imTrue      = varargin{1};
        imTrue      = normFactor*imTrue;
        imTrueNorm  = norm(imTrue(:));
        errAll      = zeros(nBreg,1);
        errBest     = inf;
        uBest       = u;  
    case 9
        imTrue      = varargin{1};
        mask        = varargin{2};
        imTrue      = normFactor*imTrue;
        imTrueNorm  = norm(imTrue(mask));
        errAll      = zeros(nBreg,1);
        errBest     = inf;
        uBest       = u;  
end % switch

% RHS of the linear system
murf = ifft2(mu*R.*f);

% Build Kernels in the Fourier space
uker = zeros(rows,cols);
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
uker = mu*R+lambda*fft2(uker)+gamma;

h   = waitbar(0,'TV for each b value...');
h2  = figure;

%  Do the reconstruction
for outer = 1:nBreg;
    for inner = 1:nInner;
        % update u
        rhs     = murf+lambda*Dxt(x-bx)+lambda*Dyt(y-by)+gamma*u;
        u       = ifft2(fft2(rhs)./uker);

        % update x and y
        dx      = Dx(u);
        dy      = Dy(u);
        [x,y]   = shrink2( dx+bx, dy+by,1/lambda);

        % update bregman parameters
        bx      = bx+dx-x;
        by      = by+dy-y;

    end % inner

    f           = f+f0-R.*fft2(u);
    murf        = ifft2(mu*R.*f);

    % Compute the error
    if nargin == 8
        errAll(outer) =  norm(u(:)-imTrue(:))/imTrueNorm;
        if (errAll(outer) <= errBest)
            uBest       = u;
            errBest     = errAll(outer);
        end %errThis 
    elseif nargin >= 9
        % Solution error norm in a mask
        errAll(outer)    = norm(u(mask)-imTrue(mask))/imTrueNorm; 
        if (errAll(outer) <= errBest)
            uBest       = u;
            errBest     = errAll(outer);
        end %errThis 
    end % nargin
%    
    if any([outer==1 outer==25 rem(outer,200)==0])
        figure(h); waitbar(outer/nBreg,h);        
        figure(h2);
        subplot(2,2,1);
        imagesc(abs(u)); colorbar; title(['Iter ' num2str(outer)]);
        subplot(2,2,2);
        imagesc(abs(x)); colorbar; title(['x, nnz(x) % ' num2str(ceil(100*nnz(x(:))/prod([rows cols])))]); colorbar;
        subplot(2,2,3);
        plot(errAll(1:outer)); axis tight; title(['Sol. error' ]);
        colormap hot;
        drawnow;          
    end

end % outer

close(h); close(h2);

if nargin >= 8
    u = uBest;
end

% undo the normalization so that results are scaled properly
u = u/normFactor;

return;

function normFactor = getNormalizationFactor(R,f)

normFactor = 1/norm(f(:)/size(R==1,1));

return;

function d = Dx(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
return

function d = Dxt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
d(:,cols) = u(:,cols)-u(:,1);
return

function d = Dy(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
return

function d = Dyt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);
d(rows,:) = u(rows,:)-u(1,:);
return

function [xs,ys] = shrink2(x,y,lambda)

s = sqrt(x.*conj(x)+y.*conj(y));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;

return;

