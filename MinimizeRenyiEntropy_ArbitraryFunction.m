% In this script we use fmincon to minimize the differential entropy for
% some optimization variables 
% clc;



% MC parameters 
N_in = [5000]; % number of theta realizations
size_N_in = size(N_in,2);
N_out = [5000]; % number of z realizations 
size_N_out = size(N_out,2); 
N_trials = 1; 

% properties of renyi entropy
alpha = 0.9;

% properties of z 
d_z = 2; 
cov_z = eye(d_z); 

% properties of theta
d_theta = d_z;               % for now, dim_z = dim_theta
mu_theta = zeros(d_theta,1); % mean of theta 
cov_theta = eye(d_theta);    % this is Cov(P) 

% rng(0); % fix random seed for eval purposes
% note: predictions only apply to linear functions f(x) = Mx + b and z ~ N(f(theta), Sigma_z) 

% optimization variables
k0 = randn(d_theta,1); % we make k which will set the diagonal transformation M = diag(k)
                       % in this particular case only, it happens to be
                       % same dim as theta. 

% mapping \mu_z = f(theta,k) = \mathcal{G}(theta,k) = Mx + b 
% Change this function to create custom \mu_z = f(theta)
f_theta = (@(theta,k,b) diag(k)*theta+b); % the mean of z is a linear function of theta here.
b = ones(d_theta,1); % randn(d_theta,1); % some offset 

% please note: in this special evaluation case we have M = diag(k0),
% i.e. f_theta = diag(k0)*theta+b. This is how we get the optimization
% variables in. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx_n_in = 1:size_N_in
    for idx_n_out = 1:size_N_out
%             tic
        N_in(idx_n_in);
        N_out(idx_n_out); 
        % call fmincon, treat b as the variables to optimize over
        %OPTIONS = optimoptions('fmincon','Algorithm','interior-point');
        tic; % clocking about 4 minutes execution in two dimensions 
        [k_fin, fval, exitflag] = fmincon(@(k) RenyiEntropy_ArbitraryFunction(alpha,mu_theta, cov_theta, d_theta, d_z, f_theta, cov_z, N_in(idx_n_in), N_out(idx_n_out), k, b, N_trials), k0) %, [],[],[],[],[],[], [], OPTIONS) % no constraints required 
        toc; 
%         diff_ent =  0.5*log((2*pi*exp(1))^d_z * det(cov_z)) % in other words the contribution of cov_mu_z is 0 
        renyi_ent_measured = fval 
    end
end
