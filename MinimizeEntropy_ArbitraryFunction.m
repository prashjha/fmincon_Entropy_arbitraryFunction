% In this script we use fmincon to minimize the differential entropy for
% some optimization variables 
clc;
clear all;
close all;

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
alphamean  =  [2.0];
alphasttd  =  [.3];
betamean  =  [4.5];
betasttd  =  [.3];
tisinput=[T1pmean; T1pstdd; T1lmean; T1lstdd; kplmean; kplstdd; kvemean; kvestdd;t0mean;t0stdd;alphamean; alphasttd; betamean ; betasttd ];

%% Variable Setup
lnsr2pi = 0.9189385332046727; % log(sqrt(2*pi))
Ntime = 30;
TR = 3;
TimeList = (0:(Ntime-1))*TR;
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
TRList = optimvar('TR',Ntime-1,1,'LowerBound',0, 'UpperBound',5);

statevariable    = optimvar('state',Nspecies,Ntime ,'LowerBound',0);
stateconstraint  = optimconstr(    [Nspecies,Ntime ]);


%% reference solution
opts.TolFun = 1e-09;
opts.TolX = 1e-09;
opts.Display = 'off';
params = struct('t0',[t0mean(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
    'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplmean(1) ;0,0],...
    'TRList',TimeList,'PerfusionTerms',[kvemean(1),0],'volumeFractions',ve,...
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
    % precompute
    gammaa = gamma(alphamean);
    %gampdfdenom = (betamean^alphamean* gamma(alphamean));
    % loop over time
    for jjj = 1:Ntime-1
%      subTR = TRList(jjj) /nsubstep;
%      expAsubTR = [ exp(-subTR*(kplqp + kveqp/ve + 1/T1Pqp)),                   0; (kplqp*exp(-subTR/T1Lqp) - kplqp*exp(-subTR*(kplqp + kveqp/ve + 1/T1Pqp)))/(kplqp + kveqp/ve - 1/T1Lqp + 1/T1Pqp), exp(-subTR/T1Lqp)];
      expATR = [ exp(-TRList(jjj)*(kplqp + kveqp/ve + 1/T1Pqp)),                   0; (kplqp*exp(-TRList(jjj)/T1Lqp) - kplqp*exp(-TRList(jjj)*(kplqp + kveqp/ve + 1/T1Pqp)))/(kplqp + kveqp/ve - 1/T1Lqp + 1/T1Pqp), exp(-TRList(jjj)/T1Lqp)];
      tk   = TRList(jjj)*(jjj  );
      tkm1 = TRList(jjj)*(jjj-1);
      aifterm   =  [ -(VIF_scale_fact(1)*((T1Pqp^2*betamean^2*ve^2*exp(kplqp*tkm1 - kplqp*tk - tk/T1Pqp + tkm1/T1Pqp + t0qp/betamean - tkm1/betamean - (kveqp*tk)/ve + (kveqp*tkm1)/ve)*((tkm1*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve)^2 - (T1Pqp^2*betamean^2*ve^2*exp(t0qp/betamean - tk/betamean)*((tk*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve)^2 + (T1Pqp*betamean*t0qp*ve*exp(t0qp/betamean)*(exp(-tk/betamean) - exp(-(kveqp*tk)/ve)*exp((kveqp*tkm1)/ve)*exp(-kplqp*tk)*exp(kplqp*tkm1)*exp(-tk/T1Pqp)*exp(tkm1/T1Pqp)*exp(-tkm1/betamean)))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve)))/(betamean^alphamean*gammaa); (VIF_scale_fact(1)*kplqp*((T1Pqp^2*betamean^2*ve^2*exp(kplqp*tkm1 - kplqp*tk - tk/T1Pqp + tkm1/T1Pqp + t0qp/betamean - tkm1/betamean - (kveqp*tk)/ve + (kveqp*tkm1)/ve)*((tkm1*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve)^2 - (T1Lqp*betamean*exp(t0qp/betamean)*exp(-tk/betamean)*(T1Lqp*betamean + T1Lqp*tk - betamean*tk))/(T1Lqp - betamean)^2 + (T1Lqp*betamean*t0qp*exp(t0qp/betamean)*(exp(-tk/betamean) - exp(-tk/T1Lqp)*exp(tkm1/T1Lqp)*exp(-tkm1/betamean)))/(T1Lqp - betamean) - (T1Pqp^2*betamean^2*ve^2*exp(t0qp/betamean - tk/betamean)*((tk*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve)^2 + (T1Pqp*betamean*t0qp*ve*exp(t0qp/betamean)*(exp(-tk/betamean) - exp(-(kveqp*tk)/ve)*exp((kveqp*tkm1)/ve)*exp(-kplqp*tk)*exp(kplqp*tkm1)*exp(-tk/T1Pqp)*exp(tkm1/T1Pqp)*exp(-tkm1/betamean)))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*kveqp + T1Pqp*betamean*kplqp*ve) + (T1Lqp*betamean*exp(-tk/T1Lqp)*exp(tkm1/T1Lqp)*exp(t0qp/betamean)*exp(-tkm1/betamean)*(T1Lqp*betamean + T1Lqp*tkm1 - betamean*tkm1))/(T1Lqp - betamean)^2))/(betamean^alphamean*gammaa*(kplqp + kveqp/ve - 1/T1Lqp + 1/T1Pqp)) ];
      %hpstatevariable(:,jjj+1) = expATR*(cos(flips(:,jjj)).*(hpstatevariable(:,jjj))) + x0.kve/ve*aifterm ;
      stateconstraint(:,jjj+1) = statevariable(:,jjj+1) ==  expATR *(cos(FaList(:,jjj)).*statevariable(:,jjj))   + kveqp/ve* aifterm ;
      % for kkk = 1:nsubstep
      %   iii = (jjj-1)*nsubstep + kkk;
      %   aifterm   = - kveqp/ve*A_inv*(eye(2) - expAsubTR)*VIF_scale_fact.*[(SubTimeList(iii)+t0qp ).^(alphamean-1).* exp(-(SubTimeList(iii)+t0qp )/ betamean) /gampdfdenom;0] ...
      %             + A_inv_sq*(expAsubTR-(A*subTR)-eye(2))*kveqp/ve/subTR*VIF_scale_fact.* [(SubTimeList(iii+1)+t0qp ).^(alphamean-1).* exp(-(SubTimeList(iii+1)+t0qp )/ betamean) /gampdfdenom-(SubTimeList(iii)+t0qp ).^(alphamean-1).* exp(-(SubTimeList(iii)+t0qp )/ betamean) /gampdfdenom;0];
      %   stateconstraint(:,iii+1) = statevariable(:,iii+1) ==  expAsubTR *(cos(subFaList(:,iii)).*statevariable(:,iii))   + aifterm ;
      % end
    end

    disp('build objective function')
    sumstatevariable =  sum(sum(sin(FaList).*(ve*statevariable  + (1-ve) *VIF_scale_fact(1)  * [(TimeList+t0qp ).^(alphamean-1).* exp(-(TimeList+t0qp )/ betamean) /(betamean^alphamean* gamma(alphamean));zeros(1,Ntime)]  ),2));

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
N_in = [3000]; % number of theta realizations
size_N_in = size(N_in,2);
N_out = [3000]; % number of z realizations 
size_N_out = size(N_out,2); 
N_trials = 1; 


% rng(0); % fix random seed for eval purposes
% note: predictions only apply to linear functions f(x) = Mx + b and z ~ N(f(theta), Sigma_z) 

% optimization variables
%k0 = randn(d_theta,1); % we make k which will set the diagonal transformation M = diag(k)
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
f_theta  = @(theta,x) MIGHQuadHPTofts(theta,x,  problem, myidx,params,Nspecies,Ntime);
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
        Fx = @(k) DifferentialEntropy_ArbitraryFunction(mu_theta, cov_theta, d_theta, d_z, f_theta, cov_z, N_in(idx_n_in), N_out(idx_n_out), k, N_trials)
        %myobjfun = Fx(k0);%[myobjfun, myobjfun_Der]= Fx(k0)
        [k_fin, fval, exitflag] = fmincon(Fx , k0,[],[],[],[],pmin,pmax,[],...
        optimset('TolX',tolx,'TolFun',tolfun,'MaxIter', ...
        maxiter,'Display','iter-detailed','Hessian',{'lbfgs',1}, ...
        'GradObj','off','PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' }) ...
            );
        toc; 
%         ent_expected =  0.5*log((2*pi*exp(1))^d_z * det(cov_z)) % in other words the contribution of cov_mu_z is 0 
        ent_measured = fval 
    end
end

function objfun =MIGHQuadHPTofts(theta,xopt,problem,myidx,params,Nspecies,Ntime)
   x0= struct();
   x0.kpl = theta(1);
   x0.kve = theta(2);
   x0.t0  = theta(3);
   x0.FaList = reshape(xopt(myidx.FaList),Nspecies,Ntime);
   x0.TR     = xopt(myidx.TR);
   T1Pqp = params.T1s(1);
   T1Lqp = params.T1s(2);
   alphamean = params.gammaPdfA(1);
   betamean = params.gammaPdfB(1);
   disp('compute state')
   tic
   % use analytic integrals
   x0.state  = zeros(Nspecies,Ntime);
   gammaa = gamma(alphamean);
   for jjj = 1:Ntime-1
     expATR = [ exp(-TR*(x0.kpl + x0.kve/ve + 1/T1Pqp)),                   0; (x0.kpl*exp(-TR/T1Lqp) - x0.kpl*exp(-TR*(x0.kpl + x0.kve/ve + 1/T1Pqp)))/(x0.kpl + x0.kve/ve - 1/T1Lqp + 1/T1Pqp), exp(-TR/T1Lqp)];
     tk   = sum(TR(1:jjj));
     tkm1 = sum(TR(1:jjj-1));
     %aifterm   =  eval([ vifone;viftwo]);
     aifterm   =  [ -(jmA0*((T1Pqp^2*betamean^2*ve^2*exp(x0.kpl*tkm1 - x0.kpl*tk - tk/T1Pqp + tkm1/T1Pqp + x0.t0/betamean - tkm1/betamean - (x0.kve*tk)/ve + (x0.kve*tkm1)/ve)*((tkm1*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)^2 - (T1Pqp^2*betamean^2*ve^2*exp(x0.t0/betamean - tk/betamean)*((tk*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)^2 + (T1Pqp*betamean*x0.t0*ve*exp(x0.t0/betamean)*(exp(-tk/betamean) - exp(-(x0.kve*tk)/ve)*exp((x0.kve*tkm1)/ve)*exp(-x0.kpl*tk)*exp(x0.kpl*tkm1)*exp(-tk/T1Pqp)*exp(tkm1/T1Pqp)*exp(-tkm1/betamean)))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)))/(betamean^alphamean*gammaa); (jmA0*x0.kpl*((T1Pqp^2*betamean^2*ve^2*exp(x0.kpl*tkm1 - x0.kpl*tk - tk/T1Pqp + tkm1/T1Pqp + x0.t0/betamean - tkm1/betamean - (x0.kve*tk)/ve + (x0.kve*tkm1)/ve)*((tkm1*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)^2 - (T1Lqp*betamean*exp(x0.t0/betamean)*exp(-tk/betamean)*(T1Lqp*betamean + T1Lqp*tk - betamean*tk))/(T1Lqp - betamean)^2 + (T1Lqp*betamean*x0.t0*exp(x0.t0/betamean)*(exp(-tk/betamean) - exp(-tk/T1Lqp)*exp(tkm1/T1Lqp)*exp(-tkm1/betamean)))/(T1Lqp - betamean) - (T1Pqp^2*betamean^2*ve^2*exp(x0.t0/betamean - tk/betamean)*((tk*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)^2 + (T1Pqp*betamean*x0.t0*ve*exp(x0.t0/betamean)*(exp(-tk/betamean) - exp(-(x0.kve*tk)/ve)*exp((x0.kve*tkm1)/ve)*exp(-x0.kpl*tk)*exp(x0.kpl*tkm1)*exp(-tk/T1Pqp)*exp(tkm1/T1Pqp)*exp(-tkm1/betamean)))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve) + (T1Lqp*betamean*exp(-tk/T1Lqp)*exp(tkm1/T1Lqp)*exp(x0.t0/betamean)*exp(-tkm1/betamean)*(T1Lqp*betamean + T1Lqp*tkm1 - betamean*tkm1))/(T1Lqp - betamean)^2))/(betamean^alphamean*gammaa*(x0.kpl + x0.kve/ve - 1/T1Lqp + 1/T1Pqp)) ];
   
     x0.state(:,jjj+1) = expATR*(cos(flips(:,jjj)).*(x0.state(:,jjj))) + x0.kve/ve*aifterm ;
   end
   toc
   % handle = figure;plot([1:size(x0.state,2)],x0.state(1,:),'k',[1:size(x0.state,2)],x0.state(2,:),'b');saveas(handle,sprintf('state%d',size(x0.state,2)),'png')
   Xfull = [ x0.FaList(:);x0.TR(:);x0.kpl;x0.kve; x0.state(:);x0.t0];
   %extraParamscon = functions(problem.nonlcon).workspace{1}.extraParams;
   %save('extraParams.mat','extraParamscon','Xfull' )
   % [initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = reducedConstraint(Xinit ,extraParamscon );
 
   
   % prevent JIT compiler by clearing the functions
   clear reducedObjective
   clear reducedConstraint 
   disp('compute obj/nonlcon')
   tic;
   [objfun ,mygradient] = problem.objective(Xfull);
   disp(sprintf('objective %f',objfun ))
   %disp('evaluate constraint')
   initConst = struct();
   %profile on
   [initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = problem.nonlcon(Xfull);
   %clear problem.nonlcon
   %myprof = profile('info')
   %save('profiledata','myprof')
   %profile off
   toc;
   %whos
   %disp('constraint complete')
   objectiveGradFA    = mygradient(myidx.FaList);
   objectiveGradTR    = mygradient(myidx.TR);
   objectiveGradState = mygradient(myidx.state);
   jacobianFA    = initConst.ceqGrad(myidx.FaList,:);
   jacobianTR    = initConst.ceqGrad(myidx.TR,:);
   jacobianState = initConst.ceqGrad(myidx.state,:);
   adjointvar =-jacobianState \objectiveGradState ;
   objfun_Der = [objectiveGradFA;objectiveGradTR] +  [jacobianFA;jacobianTR    ] *   adjointvar ;
   
end
