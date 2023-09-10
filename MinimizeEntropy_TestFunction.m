% In this script we use fmincon to minimize the differential entropy for
% some optimization variables 
% clc;

dims = 2; % the dimensions of z to consider
len_dims = size(dims,2); 

for z = 1:len_dims
    d_z = dims(z); 
    d_theta = d_z; % for now, dim_z = dim_theta
    cov_z = eye(d_z); 
    N_in = [3000]; % number of theta realizations
    size_N_in = size(N_in,2);
    N_out = [3000]; % number of z realizations 
    size_N_out = size(N_out,2); 
    N_trials = 1; 

    rng(0); % fix random seed 

    % note: predictions only apply to linear functions f(x) = Mx + b and z ~ N(f(theta), Sigma_z) 
    
    % distribution properties
    mu_theta = zeros(d_theta,1); % mean of theta 
    cov_theta = eye(d_theta); % this is Cov(P) 
    %cov_mu_z = M*eye(d_z)*M';  % this is MCov(P)M^* = the random component of the mean of Z

    % function properties
    b = ones(d_theta,1); % randn(d_theta,1); % some offset 
    % Change this function to create custom \mu_z = f(theta)
    f_theta = (@(theta,M,b) M*theta+b); % the mean of z is a linear function of theta here.

    % optimization variables
    k0 = randn(d_theta,1); % we make k which will set the diagonal transformation M = diag(k)  
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
            [k_fin, fval, exitflag] = fmincon(@(k) DifferentialEntropy_TestFunction(mu_theta, cov_theta, d_theta, d_z, f_theta, cov_z, N_in(idx_n_in), N_out(idx_n_out), k, b, N_trials), k0) %, [],[],[],[],[],[], [], OPTIONS) % no constraints required 
            toc; 
            ent_expected =  0.5*log((2*pi*exp(1))^d_z * det(cov_z)) % in other words the contribution of cov_mu_z is 0 
            ent_measured = fval 
        end
    end
end