// ****************************** Variables *******************************
var y x c cs cb hs hb hr l ls lb T ws wb bc bs r dp q q_r bar_hb a;

// y  = output 
// x  = markup 
// cs = consumption savers 
// cb = consumption borrowers 
// c  = total consumption(=y) 
// hs = housing demand savers 
// hb = housing demand borrowers 
// hr = rent housing borrowers 
// ls = labor savers 
// lb = labor borrowers 
// l  = total labor 
// wu = wage borrowers 
// wc = wage savers 
// bc = borrowing 
// bs = borrowing 
// r  = interest rate 
// dp = inflation 
// q  = housing prices 
// q_r = rent prices

// ****************************** Variables exogenas ********************** 
varexo u eps_e;

// ****************************** Parameters ******************************
parameters 
betat,          % betat = borrowers discount factor
beta,           % beta = savers discount factor
J,              % J = weight of housing on the utility function
eta,            % eta and 1/(eta-1) is the labor-supply elasticity.
K,              % K = loan-to-value ratio
gamma,          % gamma = savers labor income share.
XSS,            % XSS = steady state markup
ktil,           % ktil = markup parameter in Phillips curve
rho,            % rho = Taylor Rule, interest rate smoothing parameter
phip,           % phip =  Taylor Rule, inflation parameter
phiy,           % phiy = Taylor Rule, output parameter
tau_h, 
% lambda,
% depre_h, 
% depre_r, 
omega, 
e_h, 
rho_e,
tau_r, 
rho_a, 
SIGMA,
s,
delta;

// ****************************** Calibration *****************************
%%% betat, beta, J, eta, K, gamma, XSS, ktil, rho, phip, phiy, tau_h, 
%%% depre_h, depre_r, lambda, omega_h, e_h, tau_r, rho_a, SIGMA, delta, lambda

beta = 0.99 ;
betat = 0.98;
J = .1 ;
eta =1.2;
K = 0.9 ;
gamma = 0.4;
XSS = 1.2 ;
ktil = 0.0858 ;
rho = 0.99;
phip = 0.5 ;
phiy = 0.1 ;
tau_h = 0.01 ;
% depre_h = 0.075;
% depre_r = 0.01;
omega = 0.2;
% lambda = 0.5;
e_h = 2 ;
tau_r = -0.0185 ;
rho_a = 0.8;
rho_e = 0.95;
SIGMA = 0.003;
delta = 0.029;
s = 0.1;


// ****************************** Model ***********************************
model;
%----------- Savers   ------------------
% (1)
exp(cs)+exp(bs)-exp(q)*(1-tau_h)*(exp(hs(-1)-hs))-
exp(q_r)*exp(hr(-1)-hr)=exp(wb+lb)+exp(r(-1)+bs(-1)-dp)+exp(T)+ exp(q_r+hr);

% (2)
1/exp(cs) = beta*exp(r-dp(+1)-cs(+1)); 

% (3)
ws = (eta-1)*ls + cs;

% (4)
J/exp(hs) = (1-tau_h)*(exp(q-cs)-beta*(exp(q(+1)-cs(+1))));

% (5)
exp(q-cs) = exp(q_r-cs)+beta*exp(q(+1)-cs(+1));

%----------- Borrowers ------------------
% (6)
exp(cb)+exp(r(-1)+bc(-1)-dp)-exp(q)*(1-tau_h)*(exp(hb(-1)-hb))+
exp(q_r)*(1-tau_r)*exp(hr) = exp(wb+lb) + exp(bc);

% (7)
bc = log(K) + q(+1) + hb + dp(+1) - r;

% (8)
% 1/exp(cb) = betat*exp(r-cb(+1)-dp(+1));

% (9)
wb = (eta-1)*lb + cb;

% (10)
(J/exp(bar_hb))*((omega)*exp(bar_hb-hb))^(1/e_h)
=(1-tau_h)*(1/exp(cb))*(exp(q-K*exp(q(+1)-dp(+1)-r)))
-(1-tau_h)*(1/exp(cb(+1)))*betat*(1-K)*exp(q(+1));

% (11)
(J/exp(bar_hb))*((1-omega)*exp(bar_hb-hr))^(1/e_h)
= (1-tau_r)*exp(q_r-cb);

%----------- Market Equilibrium ---------
% (12)
exp(c) = exp(cb) + exp(cs) ;

exp(bs) = exp(bc);

% (14)
exp(l) = exp(ls) + exp(lb) ;

% (15)
1 = hs + hr + hb;

% (16)
bar_hb = (e_h/(e_h-1))*
((omega^(1/e_h))*((e_h-1)/e_h)*hb+
((1-omega)^(1/e_h))*((e_h-1)/e_h)*hr);

%------------- Production ---------------
% (17)
y = a + gamma*ls + (1-gamma)*lb ;

%------------- Labour market ------------
% (18)
ws = gamma + y - x - ls ;

% (19)
wb = 1-gamma + y - x - lb ;

%------------- Phillips curve -----------
% (20)
dp = beta*dp(+1) - ktil*(x-log(XSS));

%------------- Government ---------------
% Taylor rule
% (21)
r = rho*r(-1)+(1-rho)*(1+phip)*dp+
(1-rho)*phiy*(y-y(-1))+(1-rho)*log(1/beta) + s*eps_e;

% Equilibrium government budget constraint
% (22)
T = tau_r*exp(q_r+hr)+
tau_h*exp(q)*exp(hs-hs(-1))+
exp(hb-hb(-1));

%------------- Technology Shocks ---------

% (23)

a = rho_a*a(-1) + u;

end ;

// ****************************** Initial Values **************************
% var y x cs cb c hs hb hr ls lb l T ws wb bc bs r dp q q_r a bar_hb;

initval;
y = 0 ;
x = 0 ;
c = 0;
cs = 0.0 ;
cb = 0.0 ;
hs = 0;
hb = 0 ;
hr = 0;
ls = 0 ;
lb = 0 ;
l = 0 ;
T = 0 ;
ws = 0 ;
wb = 0 ;
bc = 0 ;
bs = 0 ;
r = 0.01 ;
dp = 0 ;
q = 0 ;
q_r = 0 ;
a = 0;
bar_hb = 0;
u = 0 ;
eps_e = 0 ;
% S = 0 ;

end;
resid(1);

options_.solve_tolf = 5e-3;
steady(solve_algo = 2, maxit = 10000);

steady ;

shocks;
var u; stderr SIGMA;
var eps_e; stderr delta;
end;

check ;

% var y x cs cb c hs hb hr ls lb l T ws wb bc bs r dp q q_r a bar_hb;

stoch_simul(order=1,irf=15) y dp r cb q q_r;