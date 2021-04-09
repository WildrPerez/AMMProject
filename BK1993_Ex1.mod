%------------------------------------------------------------------------------------------%
% MODELO de Baxter y King 1993                                                                    

%------------------------------------------------------------------------------------------%
% VARIABLES (12)                                                                        
%------------------------------------------------------------------------------------------%

var
c %consupsion
i %investment                                            
y %product
k % private capital 
n % labor
gb % basic government expenditure (gasto básico)
rk  % renta de alquiler de capital
tr % transferencias
kg % government capital stock
g % total government expenditure
ig % public investment
w % real salary
;

varexo e_gb e_ig; % e_gb: choque al gasto básico y e_ig: choque a la inversión pública

%------------------------------------------------------------------------------------------%
% PARAMETERS ()                                                                         
%------------------------------------------------------------------------------------------%
parameters
A
theta_l
theta_n
beta
delta_k
theta_k
t
theta_g
rho_gb
rho_ig
sigma_gb
sigma_ig
y_ss 
c_ss 
i_ss
w_ss
%r_ss
k_ss
n_ss %(From paper)
gb_ss
tr_ss
kg_ss
g_ss
ig_ss
r_a % annual interest rate (Dato del paper)
gy  % ratio g_ss/y_ss (Dato del paper)
igy % ratio ig_ss/y_ss (Dato del paper)
gby % ratio gb_ss/y_ss (X)
tr_y  % ratio tr_ss/y_ss (X)
kgy  % ratio kg_ss/y_ss (X)
yk  % ratio y_ss/k_ss (X)
iy  % ratio i_ss/y_ss (X)
cy  % ratio c_ss/y_ss (X)
wy  % ratio w_ss/y_ss (X)
;

%------------------------------------------------------------------------------------------%
% CALIBRATION                                                                      
%------------------------------------------------------------------------------------------%
%0)PREFERENCES y RK_SS
%--------------
r_a = 0.065; % Interest annual rate
rk_ss = (1+r_a)^(1/4);
delta = 0.025 % deprec annual rate
delta_k = (1+delta)^(1/4)-1 ; % deprec quaterly rate
beta = (rk_ss + (1-delta_k))^(-1); 
%--------------
%1)FIRMS
%--------------
A = 1; % supuesto propio.
theta_n = 0.58;
theta_k = 0.42;
%--------------
%2) GOVERNMENT
%--------------
theta_g = 0.05; % modelo base rbc (theta_g = 0), modelo con kg aumentador de productividad (theta_g = 0.05)
t = 0.2; % from paper P.320
gy = 0.2; % from paper
igy = 0.05; % from paper
%--------------
%3)SHOCKS
%--------------
rho_gb = 1 ;
sigma_gb = 0.01;
rho_ig = 0.01; 
sigma_ig = 0.01;
%------------------------------------------------------------------------------------------%
%4)STEADY STATE                                                                   
%------------------------------------------------------------------------------------------%
n_ss = 0.1; % from paper P.320
gby = gy - igy; % sum of the government investment and consumption  
tr_y = t - gy; % budget constrain of the government 
kgy = igy*(1/delta_k);
yk = rk_ss/(theta_k); % yk = ratio y_ss/k_ss (X)
iy = delta_k*yk^(-1); % ratio i_ss/y_ss (X)
cy = 1 - iy - gy;
wy = theta_n*(1-t)/n_ss;

y_ss = (A*(yk^(-theta_k))*(kgy^(theta_g))*n_ss^(theta_n))^(-1/(theta_g - theta_n));

%y_ss = A*k^{theta_k)*n^{theta_n)*kg^{theta_g);
%y_ss = A*(k_ss/y_ss)^{theta_k)*(n_ss/y_ss)^{theta_n)*(kg_ss/y_ss)^{theta_g);
%y_ss = A*(k_ss/y_ss)^{theta_k)*(n_ss/y_ss)^{theta_n)*(kg_ss/y_ss)^{theta_g);

gb_ss = gby*y_ss;
tr_ss = tr_y*y_ss;
kg_ss = kgy*y_ss;
k_ss = (yk^(-1))*y_ss;
i_ss = iy*y_ss;
c_ss = cy*y_ss;
w_ss = wy*y_ss;
g_ss = gy*y_ss;
ig_ss = igy*y_ss;

theta_l = (1-n_ss)*wy/cy;
%------------------------------------------------------------------------------------------%
%4)MODEL                                                                   
%------------------------------------------------------------------------------------------%
model;
%================================
% Familias 
%================================
theta_l/(1-n) = w*(1-0.2)/c;
1/c = beta*(1/c(+1))*(rk(+1) + (1-delta_k));
k = (1-delta_k)*k(-1) + i;
%================================
% Firmas
%================================
y = A*((k(-1))^theta_k)*(n^(theta_n))*((kg(-1))^(theta_g));
rk = theta_k*y/k(-1); % deberia contener impuestos?
w = (1-theta_n)*y/n; % deberia contener impuestos?
%================================
% Gobierno
%================================
kg = (1-delta_k)*kg(-1) + ig;
t*y = g + tr;
g = gb + ig;
%================================
% Condición de mercado 
%================================
y = c + i + g;  
%================================
% Fuentes de incertidumbre 
%================================
ln(gb/gb_ss) = rho_gb*ln(gb(-1)/gb_ss) + e_gb;
ln(ig/ig_ss) = rho_ig*ln(ig(-1)/ig_ss) + e_ig;
end;
%------------------------------------------------------------------------------------------%
% VALORES INICIALES
%------------------------------------------------------------------------------------------%
initval;
n = n_ss;
k = k_ss;
i = i_ss;
c = c_ss;
w = w_ss;
rk = rk_ss;
y = y_ss;
kg = kg_ss;
ig = ig_ss;
g = g_ss;
tr = tr_ss;
gb = gb_ss;
end;
resid(1);
steady;

%------------------------------------------------------------------------------------------%
% CHOQUES                                                                          
%------------------------------------------------------------------------------------------%
shocks;
var e_gb = (sigma_gb)^2;
var e_ig = (sigma_ig)^2;
end;

check;

%------------------------------------------------------------------------------------------%
%SIMULACIÓN: 
%------------------------------------------------------------------------------------------%
stoch_simul(order=1,irf=12,irf_plot_threshold=0)y i c gb ig n w rk;

saveas(gcf,'simul_3','pdf');
