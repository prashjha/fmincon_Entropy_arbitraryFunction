
clear all 
close all

%% reference forward solve timing
A = tril(magic(1e2));
opts.LT = true;
b = ones(size(A,2),1);
tic
x2 = linsolve(A,b,opts);
t1 = toc


%% analytic integrals
syms tau tk tkm1 t0 alphamean betamean kpl kve T1Pqp T1Lqp ve gammaa  jmA0    
vifone= jmA0/(betamean^alphamean *gammaa)*int( exp(-( kpl+kve/ve+1/T1Pqp )*(tk-tau))*exp( -(tau-t0)/betamean)*(tau - t0) ,tau,tkm1,tk)
viftwo= jmA0* kpl/(kpl+kve/ve-1/T1Lqp+1/T1Pqp)/(betamean^alphamean *gammaa)*int( (exp(-(tk-tau)/T1Lqp) - exp(-(tk-tau)*(kpl+kve/ve+1/T1Pqp )) )*exp( -(tau-t0)/betamean)*(tau - t0) ,tau,tkm1,tk)

%% Tissue Parameters
T1pmean = [ 30 ]; % s
T1pstdd = [ 10 ]; % s
T1lmean = [ 25 ]; % s
T1lstdd = [ 10 ]; % s
kplmean = [ .15 ];       % s
kplstdd = [ .03 ];       % s
kvemean = [ 0.05 ];       % s
kvestdd = [ .01  ];       % s
t0mean  = [ 0    ];       % s
t0stdd  = [ 1.3  ] ;       % s
alphamean  =  [2];
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

%% plot gamma
jmA0    = VIF_scale_fact(1);
jmalpha = alphamean(1);
jmbeta  = betamean(1);
jmt0    = t0mean(1);
gpres = [0:Ntime*TR];
jmaif   =    jmA0 * gampdf(gpres + jmt0  , jmalpha , jmbeta);
jmaifad = .5*jmA0 * 1/(betamean^2* gamma(2)) *(gpres+jmt0).^(2-1).* exp(-(gpres+jmt0)/ betamean);
jmaifad2= .2*jmA0 * 1/(betamean^2* gamma(3)) *(gpres+jmt0).^(3-1).* exp(-(gpres+jmt0)/ betamean);
jmaifad3=    jmA0 * 1/(betamean^2* gamma(1.5)) *(gpres+jmt0).^(1.5-1).* exp(-(gpres+jmt0)/ betamean);
figure(4)
plot(gpres,jmaif ,'b', gpres,jmaifad ,'r',gpres,jmaifad2,'g',gpres,jmaifad3,'k')
ylabel('aif')
xlabel('sec')

    
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
tic 
[t_axis,Mxy,Mz] = model.compile(M0.',params);
t2 = toc
handle = figure(1);plot(TR_list,Mz(1,:),'k',TR_list,Mz(2,:),'b'); hold on;

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
    % precompute
    gampdfdenom = (betamean^alphamean* gamma(alphamean));
    % loop over time
    for jjj = 1:Ntime-1
      subTR = TRList(jjj) /nsubstep;
      expAsubTR = [ exp(-subTR*(kplqp + kveqp/ve + 1/T1Pqp)),                   0; (kplqp*exp(-subTR/T1Lqp) - kplqp*exp(-subTR*(kplqp + kveqp/ve + 1/T1Pqp)))/(kplqp + kveqp/ve - 1/T1Lqp + 1/T1Pqp), exp(-subTR/T1Lqp)];
      for kkk = 1:nsubstep
        iii = (jjj-1)*nsubstep + kkk;
        aifterm   = - kveqp/ve*A_inv*(eye(2) - expAsubTR)*VIF_scale_fact.*[(SubTimeList(iii)+t0qp ).^(alphamean-1).* exp(-(SubTimeList(iii)+t0qp )/ betamean) /gampdfdenom;0] ...
                  + A_inv_sq*(expAsubTR-(A*subTR)-eye(2))*kveqp/ve/subTR*VIF_scale_fact.* [(SubTimeList(iii+1)+t0qp ).^(alphamean-1).* exp(-(SubTimeList(iii+1)+t0qp )/ betamean) /gampdfdenom-(SubTimeList(iii)+t0qp ).^(alphamean-1).* exp(-(SubTimeList(iii)+t0qp )/ betamean) /gampdfdenom;0];
        auxvariable(:,iii+1) = expAsubTR*(cos(subFaList(:,iii)).*(auxvariable(:,iii))) + aifterm          ;
        stateconstraint(:,iii+1) = statevariable(:,iii+1) ==  expAsubTR *(cos(subFaList(:,iii)).*statevariable(:,iii))   + aifterm ;
      end
    end


x0.kpl = kplmean ;
x0.kve = kvemean;
x0.t0  = t0mean;
x0.FaList = flips;
x0.TR     = TR* ones(Ntime-1,1);
tic
state  = evaluate(auxvariable ,x0);
t3 = toc
handle = figure(1);plot(TR_list, cos(flips(1,:)).* state(1,1:nsubstep:(Ntime-1)*nsubstep+1) ,'k--',TR_list, cos(flips(2,:)).* state(2,1:nsubstep:(Ntime-1)*nsubstep+1) ,'b--');


%% use analytic integrals
hpstatevariable    = zeros(Nspecies,Ntime);
tic
gammaa = gamma(alphamean);
for jjj = 1:Ntime-1
  expATR = [ exp(-TR*(x0.kpl + x0.kve/ve + 1/T1Pqp)),                   0; (x0.kpl*exp(-TR/T1Lqp) - x0.kpl*exp(-TR*(x0.kpl + x0.kve/ve + 1/T1Pqp)))/(x0.kpl + x0.kve/ve - 1/T1Lqp + 1/T1Pqp), exp(-TR/T1Lqp)];
  tk   = TR*(jjj+1);
  tkm1 = TR*jjj;
  aifterm   =  [ -(jmA0*((T1Pqp^2*betamean^2*ve^2*exp(x0.kpl*tkm1 - x0.kpl*tk - tk/T1Pqp + tkm1/T1Pqp + x0.t0/betamean - tkm1/betamean - (x0.kve*tk)/ve + (x0.kve*tkm1)/ve)*((tkm1*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)^2 - (T1Pqp^2*betamean^2*ve^2*exp(x0.t0/betamean - tk/betamean)*((tk*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)^2 + (T1Pqp*betamean*x0.t0*ve*exp(x0.t0/betamean)*(exp(-tk/betamean) - exp(-(x0.kve*tk)/ve)*exp((x0.kve*tkm1)/ve)*exp(-x0.kpl*tk)*exp(x0.kpl*tkm1)*exp(-tk/T1Pqp)*exp(tkm1/T1Pqp)*exp(-tkm1/betamean)))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)))/(betamean^alphamean*gammaa); (jmA0*x0.kpl*((T1Pqp^2*betamean^2*ve^2*exp(x0.kpl*tkm1 - x0.kpl*tk - tk/T1Pqp + tkm1/T1Pqp + x0.t0/betamean - tkm1/betamean - (x0.kve*tk)/ve + (x0.kve*tkm1)/ve)*((tkm1*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)^2 - (T1Lqp*betamean*exp(x0.t0/betamean)*exp(-tk/betamean)*(T1Lqp*betamean + T1Lqp*tk - betamean*tk))/(T1Lqp - betamean)^2 + (T1Lqp*betamean*x0.t0*exp(x0.t0/betamean)*(exp(-tk/betamean) - exp(-tk/T1Lqp)*exp(tkm1/T1Lqp)*exp(-tkm1/betamean)))/(T1Lqp - betamean) - (T1Pqp^2*betamean^2*ve^2*exp(x0.t0/betamean - tk/betamean)*((tk*(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve))/(T1Pqp*betamean*ve) - 1))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve)^2 + (T1Pqp*betamean*x0.t0*ve*exp(x0.t0/betamean)*(exp(-tk/betamean) - exp(-(x0.kve*tk)/ve)*exp((x0.kve*tkm1)/ve)*exp(-x0.kpl*tk)*exp(x0.kpl*tkm1)*exp(-tk/T1Pqp)*exp(tkm1/T1Pqp)*exp(-tkm1/betamean)))/(betamean*ve - T1Pqp*ve + T1Pqp*betamean*x0.kve + T1Pqp*betamean*x0.kpl*ve) + (T1Lqp*betamean*exp(-tk/T1Lqp)*exp(tkm1/T1Lqp)*exp(x0.t0/betamean)*exp(-tkm1/betamean)*(T1Lqp*betamean + T1Lqp*tkm1 - betamean*tkm1))/(T1Lqp - betamean)^2))/(betamean^alphamean*gammaa*(x0.kpl + x0.kve/ve - 1/T1Lqp + 1/T1Pqp)) ];
  hpstatevariable(:,jjj+1) = expATR*(cos(flips(:,jjj)).*(hpstatevariable(:,jjj))) + x0.kve/ve*aifterm ;
end
t4 = toc
handle = figure(1);plot(TR_list, cos(flips(1,:)).* hpstatevariable(1,:) ,'k.-',TR_list, cos(flips(2,:)).* hpstatevariable(2,:) ,'b.-');
