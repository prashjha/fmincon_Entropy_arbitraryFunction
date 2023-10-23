% In this script we use fmincon to minimize the differential entropy for
% some optimization variables 
clc;
clear all;

%% Tissue Parameters
T1pmean = [ 30 ]; % s
T1pstdd = [ 10 ]; % s
T1lmean = [ 25 ]; % s
T1lstdd = [ 10 ]; % s
kplmean = [ .15 ];       % s
kplstdd = [ .03 ];       % s
kvemean = [ 0.05 ];       % s
kvestdd = [ .01  ];       % s
t0mean  = [ 4    ];       % s
t0stdd  = [ 1.3  ] ;       % s
alphamean  =  [2.5];
alphasttd  =  [.3];
betamean  =  [4.5];
betasttd  =  [.3];
tisinput=[T1pmean; T1pstdd; T1lmean; T1lstdd; kplmean; kplstdd; kvemean; kvestdd;t0mean;t0stdd;alphamean; alphasttd; betamean ; betasttd ];

%% Variable Setup
lnsr2pi = 0.9189385332046727; % log(sqrt(2*pi))
Ntime = 30;
TR = 3;
TR_list = (0:(Ntime-1))*TR;
nsubstep = 3;
SubTimeList = (0:(Ntime-1)*nsubstep )*TR/nsubstep;
M0 = [0,0];
ve = 0.95;
%ve = 1.;
VIF_scale_fact = [100;0];
    
%% setup optimization variables
Nspecies = 2
FaList = optimvar('FaList',Nspecies,Ntime,'LowerBound',0, 'UpperBound',35*pi/180);
subFaList = optimexpr(      [Nspecies,(Ntime-1)*nsubstep+1 ]);
subFaList(:,1:nsubstep:(Ntime-1)*nsubstep+1) = FaList ;
TRList = optimvar('TR',Ntime-1,1,'LowerBound',0, 'UpperBound',5);%TR_list;

statevariable    = optimvar('state',Nspecies,(Ntime-1)*nsubstep+1 ,'LowerBound',0);
auxvariable      = optimexpr(      [Nspecies,(Ntime-1)*nsubstep+1 ]);
stateconstraint  = optimconstr(    [Nspecies,(Ntime-1)*nsubstep+1 ]);


%% reference solution
opts.TolFun = 1e-09;
opts.TolX = 1e-09;
opts.Display = 'off';
params = struct('t0',[t0mean(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
    'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplmean(1) ;0,0],...
    'TRList',TR_list,'PerfusionTerms',[kvemean(1),0],'volumeFractions',ve,...
    'fitOptions', opts);
model = HPKinetics.NewMultiPoolTofftsGammaVIF();
for n = 1:Ntime
    flips(2,n) = 30*pi/180;
    flips(1,n) = 20*pi/180;
end
params.FaList = flips;
[t_axis,Mxy,Mz] = model.compile(M0.',params);


modelSNR = 20 ; % TODO - FIXME
signuImage = (max(Mxy(1,:))+max(Mxy(2,:)))/2/modelSNR;
% walker paper is peak pyruvate only
signuImage = max(Mxy(1,:))/modelSNR;
% variance for Gauss RV is sum. sqrt for std
signu = sqrt(2* Ntime) * signuImage;

% properties of z 
d_z = 1; 
cov_z = signu* eye(d_z); 


disp('build state variable')
stateconstraint(:,1)  = statevariable(:,1) ==0;
auxvariable(:,1) =0;

    NumberUncertain=3;
    switch (NumberUncertain)
       case(3)
         T1Pqp   = T1pmean;
         T1Lqp   = T1lmean;
         kplqp   = optimvar('kpl');
         klpqp   =    0 ;     % @cmwalker where do I get this from ? 
         kveqp   = optimvar('kve');
         t0qp    = optimvar('t0');
         % properties of theta
         d_theta = NumberUncertain;               % for now, dim_z = dim_theta
         mu_theta = [kplmean;kvemean ;t0mean  ]; %zeros(d_theta,1); % mean of theta 
         cov_theta = diag([kplstdd ;kvestdd ;t0stdd  ]);    % this is Cov(P) 
       %case(4)
       %  T1Pqp   = xn{1}(iqp);
       %  T1Lqp   = xn{2}(iqp);
       %  kplqp   = xn{3}(iqp);
       %  klpqp   =    0 ;     % @cmwalker where do I get this from ? 
       %  kveqp   = xn{4}(iqp);
       %  t0qp    = t0mean(1); 
    end 

    % precompute
    A = [(-kplqp -kveqp/ve -1/T1Pqp), 0; kplqp, -1/T1Lqp ]
    % syms a  kpl d subTR    T1P kveqp T1L 
    % A_inv = inv([a,  0; kpl, d ])
    % 
    % A_inv =
    % 
    % [       1/a,   0]
    % [-kpl/(a*d), 1/d]
    % >> a = -1/T1P - kpl - kveqp/ve
    % >> d = -1/T1L
    % >> eval(A_inv )
    %    ans =
    %    
    %    [        -1/(kpl + kveqp/ve + 1/T1P),    0]
    %    [-(T1L*kpl)/(kpl + kveqp/ve + 1/T1P), -T1L]
    A_inv = [        -1/(kplqp + kveqp/ve + 1/T1Pqp),    0; -(T1Lqp*kplqp)/(kplqp + kveqp/ve + 1/T1Pqp), -T1Lqp];
    A_inv_sq = A_inv^2
    % >> syms a  kpl d subTR    T1P kveqp T1L 
    % >> expATR = expm([a,  0; kpl, d ] * subTR )
    % 
    % expATR =
    % 
    % [                                     exp(a*subTR),                0]
    % [(kpl*exp(a*subTR) - kpl*exp(subTR*d))/(a - d), exp(subTR*d)]
    % 
    % >> a = -1/T1P - kpl - kveqp
    % >> d = -1/T1L
    % >> eval(expATR)
    % 
    % ans =
    % 
    % [                                                              exp(-subTR*(kpl + kveqp + 1/T1P)),                   0]
    % [(kpl*exp(-subTR/T1L) - kpl*exp(-subTR*(kpl + kveqp + 1/T1P)))/(kpl + kveqp - 1/T1L + 1/T1P), exp(-subTR/T1L)]
    %    

    % loop over time
    for jjj = 1:Ntime-1
      subTR = TRList(jjj) /nsubstep;
      expAsubTR = [ exp(-subTR*(kplqp + kveqp/ve + 1/T1Pqp)),                   0; (kplqp*exp(-subTR/T1Lqp) - kplqp*exp(-subTR*(kplqp + kveqp/ve + 1/T1Pqp)))/(kplqp + kveqp/ve - 1/T1Lqp + 1/T1Pqp), exp(-subTR/T1Lqp)];
      for kkk = 1:nsubstep
        iii = (jjj-1)*nsubstep + kkk;
        aifterm   = - kveqp/ve*A_inv*(eye(2) - expAsubTR)*VIF_scale_fact.*[(SubTimeList(iii)+t0qp ).^(alphamean-1).* exp(-(SubTimeList(iii)+t0qp )/ betamean) /(betamean^alphamean* gamma(alphamean));0] ...
                  + A_inv_sq*(expAsubTR-(A*subTR)-eye(2))*kveqp/ve/subTR*VIF_scale_fact.* [(SubTimeList(iii+1)+t0qp ).^(alphamean-1).* exp(-(SubTimeList(iii+1)+t0qp )/ betamean) /(betamean^alphamean* gamma(alphamean))-(SubTimeList(iii)+t0qp ).^(alphamean-1).* exp(-(SubTimeList(iii)+t0qp )/ betamean) /(betamean^alphamean* gamma(alphamean));0];
        auxvariable(:,iii+1) = expAsubTR*(cos(subFaList(:,iii)).*(auxvariable(:,iii))) + aifterm          ;
        stateconstraint(:,iii+1) = statevariable(:,iii+1) ==  expAsubTR *(cos(subFaList(:,iii)).*statevariable(:,iii))   + aifterm ;
      end
    end

    disp('build objective function')
    sumstatevariable =  sum(sum(sin(subFaList).*(ve*statevariable  + (1-ve) *VIF_scale_fact(1)  * [(SubTimeList+t0qp ).^(alphamean-1).* exp(-(SubTimeList+t0qp )/ betamean) /(betamean^alphamean* gamma(alphamean));zeros(1,(Ntime-1)*nsubstep+1)]  ),2));

    %% 
    % Create an optimization problem using these converted optimization expressions.
    
    disp('create optim prob')
    convprob = optimproblem('Objective',sumstatevariable , "Constraints",stateconstraint);
    myidx = varindex(convprob )
    %% 
    % View the new problem.
    
    %show(convprob)
    problem = prob2struct(convprob,'ObjectiveFunctionName','reducedObjective','ConstraintFunctionName','reducedConstraint');


% MC parameters 
N_in = [50]; % number of theta realizations
size_N_in = size(N_in,2);
N_out = [100]; % number of z realizations 
size_N_out = size(N_out,2); 
N_trials = 1; 


% rng(0); % fix random seed for eval purposes
% note: predictions only apply to linear functions f(x) = Mx + b and z ~ N(f(theta), Sigma_z) 

% optimization variables
%k0 = randn(d_theta,1); % we make k which will set the diagonal transformation M = d iag(k)
k0 =  [flips(:);TR* ones(Ntime-1,1) ];   
pmin =  [flips(:)*0;zeros(Ntime-1,1)   ];     
pmax =  [flips(:)*0+35*pi/180;5*ones(Ntime-1,1) ];
                       % in this particular case only, it happens to be
                       % same dim as theta. 

% solver options
    tolx=1.e-9;
    tolfun=5.e-4;
    maxiter=400;

% mapping \mu_z = f(theta,k) = \mathcal{G}(theta,k) = Mx + b 
% Change this function to create custom \mu_z = f(theta)
%f_theta = (@(theta,k,b) diag(k)*theta+b); % the mean of z is a linear function of theta here.
%b = ones(d_theta,1); % randn(d_theta,1); % some offset 
%% signal model
% [f_theta, f_theta_grad]  = @(theta,x) MIGHQuadHPTofts(theta,x,  problem, myidx,auxvariable,Nspecies,Ntime);
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
%         options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
        tic; % clocking about 4 minutes execution in two dimensions 
        A = []; c = []; Aeq=[]; ceq=[]; lb=[]; ub=[]; nonlcon =[];
        Fx = @(k) DifferentialEntropy_QuadHPT(mu_theta, cov_theta, d_theta, d_z, cov_z, N_in(idx_n_in), N_out(idx_n_out), k, N_trials, problem, myidx, auxvariable,Nspecies,Ntime);
        %myobjfun = Fx(k0);%[myobjfun, myobjfun_Der]= Fx(k0)
        [k_fin, fval, exitflag] = fmincon(Fx , k0, [],[],[],[],pmin,pmax,[],...
        optimset('TolX',tolx,'TolFun',tolfun,'MaxIter', ...
        maxiter,'Display','iter-detailed','Hessian',{'lbfgs',1}, ...
        'GradObj','on','PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' }) ...
            );
        toc; 
%         ent_expected =  0.5*log((2*pi*exp(1))^d_z * det(cov_z)) % in other words the contribution of cov_mu_z is 0 
        ent_measured = fval 
    end
end
