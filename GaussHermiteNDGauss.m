function [x,xn,xm,w,wn] = GaussHermiteNDGauss(n,mu,sigma)

% This function determines the abscisas (x) and weights (w) for the
% Gauss-Hermite quadrature of order n>1, on the interval [-INF, +INF].
    % This function is valid for any degree n>=2, as the companion matrix
    % (of the n'th degree Hermite polynomial) is constructed as a
    % symmetrical matrix, guaranteeing that all the eigenvalues (roots)
    % will be real.

% Input:
    % n = Number of quadrature points, i.e. degree of quadrature.
    % mu = Mean of Gaussian distribution. 
    % sigma = standard dev of Gaussian distribution. If diagonal, may
    % be entered as a vector. Otherwise, must be square and dimension must
    % match mu.
    
% Output:
    % x = Quadrature abscisas for standard normal distribution.
    % xn = Quadrature abscisas in ndgrid format. Each dimension is in a
    % separate cell.
    % xm = Quadrature abscisas for each dimension. The matrix will be n x length(mu).
    % w = 1-D quadrature weights.
    % wn = Normalized N-D quadrature weights.
    
% Error Handling
if  ~isvector(mu) || ~ismatrix(sigma) || (min(size(sigma))~=1 && size(sigma,1)~=size(sigma,2)) || length(mu)~=length(sigma)
    disp('Input error.');
    return;
end
if size(mu,1)~=1
    mu=mu';
end
if size(sigma,1)==size(sigma,2)
    sigma=diag(sigma)';
end
if size(sigma,1)~=1
    sigma=sigma';
end

% 1-D Gauss-Hermite Quadrature Points and Weights Solve
i   = 1:n-1;
a   = sqrt(i/2);
CM  = diag(a,1) + diag(a,-1);

[V, L]   = eig(CM);
[x, ind] = sort(diag(L));
V        = V(:,ind)';
w        = sqrt(pi) * V(:,1).^2;

% Generate N-D Quadrature Point Grid for Gaussian Function
ndim=length(mu);
xn=cell(1,ndim);
xn{1}=sqrt(2)*kron(x,sigma)+repmat(mu,n,1);
xm=xn{1};
xn=mat2cell(xn{1},n,ones(1,ndim));
[xn{:}]=ndgrid(xn{:});

% Generate N-D Weighting Grid for Gaussian Function
wn=w;
for iii=1:ndim-1
    wn=kron(wn,w);
end
if ndim>1
    wn=reshape(wn,repmat(n,1,ndim));
end
wn=wn/pi^(ndim/2);
