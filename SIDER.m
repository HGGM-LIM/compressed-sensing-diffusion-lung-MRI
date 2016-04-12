function [u,errAll] = SIDER(Dest,Alphaest,b,f,R,N,mu,lambda,gamma,alpha,beta,nBreg,mask,varargin)
% function [u,errAll] =
% SIDER(Dest,Alphaest,b,f,R,N,mu,lambda,gamma,alpha,beta,nBreg,mask,var
% argin)
%
% The proposed method incorporates the knowledge of the signal decay into
% the reconstruction (SIDER) to accelerate the acquisition of MR diffusion
% data by undersampling in both spatial and b-value dimensions.
%
% SIDER combines total variation (TV) with a
% penalty function that promotes sparsity across the b-direction as follows
% min_u beta|grad_x,y u|_1 + gamma|M u|_1 st. ||Fu-f||^2 < delta, where the
% first term corresponds to spatial TV and M is is an operator that encodes
% the relationship between ventilation images for consecutives values of b,
% based on a stretched exponential model.
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

switch nargin
    case 14
        uTarget     = varargin{1};      
end % nargin

tolKrylov   = 1e-2;% 1e-2 % Krylov convergence criterion

rows        = N(1);
cols        = N(2);
numTime     = N(3);
Nf          = size(f);

R           = reshape(R,Nf);
f           = f.*R;

% Normalize data
normFactor  = getNormalizationFactor(R,f);
f           = normFactor*f;
uTarget     = uTarget*normFactor;

% Calculate the norm of the target on the given mask 
% for ip = 1:numTime
%     uTargetNorm(ip) = norm(reshape(uTarget(:,:,ip),[],1));
% end
for ip = 1:numTime
    tmp             = uTarget(:,:,ip);
    uTargetNorm(ip) = norm(tmp(mask));
end

% Signal decay operator
dm          = -[exp(-( (Dest*b(2:end)).^Alphaest - (Dest*b(1:end-1)).^Alphaest )) 0]';
L           = spdiags([dm ones(N(3),1)],[-1 0],N(3),N(3));

LT          = L';
LTL         = LT*L;

LMat        = @(x)reshape((L*reshape(x,[prod(N(1:2)),N(3)])')',N);
LTMat       = @(x)reshape((LT*reshape(x,[prod(N(1:2)),N(3)])')',N);
LTLMat      = @(x)reshape((LTL*reshape(x,[prod(N(1:2)),N(3)])')',N);

uBest       = zeros(rows,cols,numTime);
errBest     = inf;

% Reserve memory for the auxillary variables
f0      = f;
v       = zeros(rows,cols,numTime);
u       = zeros(rows,cols,numTime);
x       = zeros(rows,cols,numTime);
y       = zeros(rows,cols,numTime);
bx      = zeros(rows,cols,numTime);
by      = zeros(rows,cols,numTime);
dx      = zeros(rows,cols,numTime);
dy      = zeros(rows,cols,numTime);

p       = zeros(rows,cols,numTime);
bp      = zeros(rows,cols,numTime);

murf    = mu*ifft2(f);

h   = waitbar(0,'SIDER');
h2  = figure;

%  Do the reconstruction
for outer = 1:nBreg
    rhs_wt  = gamma*u;
    rhs_p   = lambda*LTMat(p-bp);
    rhs_tv  = lambda*(Dxt(x-bx)+Dyt(y-by));
    rhs     = murf+rhs_wt+rhs_p+rhs_tv;
    u       = reshape(krylov(rhs(:)),N);
    
    % Derivatives
    dx      = Dx(u);
    dy      = Dy(u);
    dp      = LMat(u);
    
    % update x and y
    [x,y]   = shrink2(dx+bx, dy+by,alpha/lambda);
    p       = shrink1(dp+bp, beta/lambda);
    
    % update bregman parameters
    bx      = bx+dx-x;
    by      = by+dy-y;
    bp      = bp+dp-p;
    
    % Bregman iteration for the data constraint
    fForw   = fft2(u).*R;
    f       = f + f0-fForw;
    murf    = mu*(ifft2(f));
        
%     for ip = 1:numTime
%         % Solution error norm    
%         tmp                 = u(:,:,ip)-uTarget(:,:,ip);
%         errAll(outer,ip)    = norm(tmp(:))/uTargetNorm(ip);
%     end
    % Solution error norm in a mask
    for ip = 1:numTime 
        tmp                 = u(:,:,ip);
        tmp2                = uTarget(:,:,ip);
        errAll(outer,ip)    = norm(tmp(mask)-tmp2(mask))/uTargetNorm(ip);
    end
    
    if nargin >= 14
        if any([outer ==1, outer == 100, outer == 500, rem(outer,250)==0])
            figure(h); waitbar(outer/nBreg,h);        
            ip = 1;
            close(h2);
            h2=figure;
            subplot(2,2,1);
            imagesc(abs(u(:,:,ip))); colorbar; title(['Iter ' num2str(outer)]);
            subplot(2,2,2);
            imagesc(abs(x(:,:,1))); colorbar; title(['x, nnz(x) % ' num2str(ceil(100*nnz(x(:))/prod(N)))]); colorbar;
            subplot(2,2,3);
            imagesc(abs(p(:,:,ip))); colorbar; title(['p, nnz(p) % ' num2str(ceil(100*nnz(p(:))/prod(N)))]); colorbar;
            subplot(2,2,4); plot(errAll(1:outer,:)); axis tight; title(['Sol. error' ]);
            colormap hot;
            drawnow;            
        end % rem
    end % plots
    
    if (mean(errAll(outer,:)) <= errBest)
        uBest       = u;
        errBest     = mean(errAll(outer,:));
    end %errThis
    
end % outer

close(h2); close(h);

if (nargin >= 15)
    u = uBest;
end % nargin

% undo the normalization so that results are scaled properly
u = u/normFactor;

% =====================================================================
% =====================================================================
    function normFactor = getNormalizationFactor(R,f)        
        normFactor = 1/norm(f(:)/size(R==1,1));        
    end

    function d = Dx(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,2:cols,:) = u(:,2:cols,:)-u(:,1:cols-1,:);
        d(:,1,:) = u(:,1,:)-u(:,cols,:);
    end

    function d = Dxt(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(:,1:cols-1,:) = u(:,1:cols-1,:)-u(:,2:cols,:);
        d(:,cols,:) = u(:,cols,:)-u(:,1,:);
    end

    function d = Dy(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(2:rows,:,:) = u(2:rows,:,:)-u(1:rows-1,:,:);
        d(1,:,:) = u(1,:,:)-u(rows,:,:);
    end

    function d = Dyt(u)
        [rows,cols,height] = size(u);
        d = zeros(rows,cols,height);
        d(1:rows-1,:,:) = u(1:rows-1,:,:)-u(2:rows,:,:);
        d(rows,:,:) = u(rows,:,:)-u(1,:,:);
    end
 
    function [xs,ys] = shrink2(x,y,lambda)        
        s = sqrt(x.*conj(x)+y.*conj(y));
        ss = s-lambda;
        ss = ss.*(ss>0);        
        s = s+(s<lambda);
        ss = ss./s;        
        xs = ss.*x;
        ys = ss.*y;
    end

    function xs = shrink1(x,lambda)        
        s = abs(x);
        xs = sign(x).*max(s-lambda,0);
    end

    % Krylov solver subroutine
    function dx = krylov(r)
        %dx = gmres (@jtjx, r, 30, tolKrylov, 100);
        [dx,flag] = bicgstab(@jtjx, r, tolKrylov, 100);
    end

    % Callback function for matrix-vector product (called by krylov)
    function b = jtjx(sol)
        solMat  = reshape(sol,N);
        
        % Laplacian part
        btv     = lambda*(Dyt(Dy(solMat)) + Dxt(Dx(solMat)));
        bP      = lambda*(LTLMat(solMat));
        
        % Jacobian u part
        bF      = mu*(ifft2(R.*fft2(solMat)));
        bwt     = gamma*sol;
        
        b       = bwt(:) + btv(:) + bP(:) + bF(:);
    end
% =====================================================================
end
%
