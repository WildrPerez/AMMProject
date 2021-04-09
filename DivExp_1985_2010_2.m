% --------------------------------------------------------------------------------------------------
% By César Salinas and Wilder Pérez, May 2018
% --------------------------------------------------------------------------------------------------
clear all
close all
clc


linea = '-----------------------------------------------------------------------------------------------------------';
disp(linea)
fprintf(1, '"Bayesian model averaging" aplicada a la diversificación exportadora.');
fprintf(1, '\n');
disp(linea)
fprintf(1, 'Basado en Fernandez, Ley y Steel (2001), "Benchmark Priors for Bayesian Model Averaging",\n');
fprintf(1, '                                          Journal of Econometrics, 100(2), 381–427.\n');
fprintf(1, '          Fernandez, Ley y Steel (2001), "Model Uncertainty in Cross-Country Growth Regressions",\n');
fprintf(1, '                                          Journal of Applied Econometrics, 16(5), 563-576.\n');
disp(linea), disp(' ')


% Loading data: See the MS Excel file
%--------------------------------------------------------------------------
[A, B] = xlsread('bd_1985_2010_2');


for i = 1:size(A, 2)
    vari = B{1, i};
    data = A(:, i);
    eval([ vari ' = data;']);
end

Yini   = rtfpna5;                                      
%Xini   = [pop5 hc5 rgdpo_pc5 ck5 xr5 lp5 lpc5 ls5 lsc5 yr_sch5 yr_sch_pri5 yr_sch_sec5 pve5 gee5 rqe5 cce5 polity5 peg5 agl5 agmach5 agmachpal5 agvalad5 arabl5 br5 bm5 coos5 cab5 cbfs5 ctps5 ctpsbb5 epe5 ese5 to5 fdi5 gdppc5 ggce5 geet5 gcf5 gds5 gfcf5 gne5 gs5 icp5 igdpd5 mva5 me5 oer5 tr5 nrr5 ts5 ped5 sed5 ses5 sep5 fe5];
%Xnames = char('pop5', 'hc5', 'rgdpo_pc5', 'ck5', 'xr5', 'lp5', 'lpc5', 'ls5', 'lsc5', 'yr_sch5', 'yr_sch_pri5', 'yr_sch_sec5', 'pve5', 'gee5', 'rqe5', 'cce5', 'polity5', 'peg5', 'agl5', 'agmach5', 'agmachpal5', 'agvalad5', 'arabl5', 'br5', 'bm5', 'coos5', 'cab5', 'cbfs5', 'ctps5', 'ctpsbb5', 'epe5', 'ese5', 'to5', 'fdi5', 'gdppc5', 'ggce5', 'geet5', 'gcf5', 'gds5', 'gfcf5', 'gne5', 'gs5', 'icp5', 'igdpd5', 'mva5', 'me5', 'oer5', 'tr5', 'nrr5', 'ts5', 'ped5', 'sed5', 'ses5', 'sep5', 'fe5');
Xini    = [r_na5 r_as5 r_eu5 r_af5 r_sa5 lsc5 lpc5 pop5 hc5 rgdpo_pc5 rgdpna_pc5 xr5 ck5 agl5 arabl5 sep5 ses5 ped5 sed5 polity5 nrr5 trad5 ctps5 bm5 gcf5 ts5 icp5 ggce5 fdi5 peg5];
Xnames  = char('r_na5', 'r_as5', 'r_eu5', 'r_af5', 'r_sa5', 'lsc5' , 'lpc5' , 'pop5', 'hc5', 'rgdpo_pc5', 'rgdpna_pc5', 'xr5', 'ck5', 'agl5', 'arabl5', 'sep5', 'ses5', 'ped5', 'sed5', 'polity5', 'nrr5', 'trad5', 'ctps5', 'bm5', 'gcf5', 'ts5', 'icp5', 'ggce5', 'fdi5', 'peg5');

k   = size(Xini, 2);  % # of regressors
N   = length(iso3n);

% Dealing with missing values 
%--------------------------------------------------------------------------
[D0, N1] = DELMISS([Yini Xini], N); 
yraw = D0(:, 1);     
Xraw = D0(:, 2:k+1);  

[n, K] = size(Xraw);

uno = ones(n ,1);
y = 100*yraw;
Z = Xraw - uno*mean(Xraw);
yM1y = sum( (y - mean(y)).^2 );

% MCO para determinar un modelo inicial
X     = [Z uno];
XXinv = (X'*X)\eye(K + 1);
b     = XXinv*X'*y;
e     = y - X*b;
Vb    = diag(XXinv)*(e'*e)/(n - K);
tb    = abs( b(1:K)./sqrt(Vb(1:K)) );

% Modelo actual
Ms = (tb > 1)';
Ms(K + 1) = true; % Siempre se incluye el intercepto en las regresiones


% Parámetros de simulación
R      = 300000;        % Numero FINAL de repeticiones
burnin = R/5;        % Número de repeticiones "burn in"
iter   = R + burnin;

MOD_   = nan(R, K + 1);
BETA_   = zeros(R, K + 1);
VBETA_  = zeros(R, K + 1); % Solo guardamos las varianzas
logpM_ = nan(R, 1);
SIGMA2_ = nan(R, 1);

logPosts  = -1/eps; % "Aceptamos" la primera iteración.
Acc    = 0;   % Repeticiones donde la realización candidata es aceptada
Xbig   = X;


fprintf('Algoritmo MC^3:\n')
for i = 1:iter
    % Tomamos un número aleatorio entre 0, 1, 2, 3, ..., K
    h = round(K*rand);
    Msold = Ms;
    if h == 0
        % No hay nada que hacer, se mantiene el modelo
        Acc = Acc + 1;
    else
        if Ms(h)            % Si "h" se encuentra en el modelo actual, ...
            Ms(h) = false;  %    ... lo eliminamos
        else                % Si "h" no se encuentra en el modelo actual, ...
            Ms(h) = true ;  %    ... lo agregamos
        end
        
        % Selección del modelo
        k = sum(Ms) - 1; % No contamos al intercepto
        X = Xbig(:, Ms);
        
        % Regresión MCO
        XXinv = (X'*X)\eye(k + 1);
        b     = XXinv*X'*y;
        e     = y - X*b;
        yMy   = e'*e;
        
        % Determinación de g y momentos a posteriori
        g = 1/max([n K^2]);
        S = ( yMy + g*yM1y )/(1 + g);
        
        % log marginal likelihood
        logPost = ( k*log(g) - k*log(1 + g) - (n - 1)*log(S) )/2;
        
        % Aceptación o rechazo
        if log(rand) < min([0 logPost-logPosts])
            Ebeta = b/(1 + g);
            Vbeta = XXinv*S/(n - 2)/(1 + g);
            logPosts = logPost;
            Acc = Acc + 1;
        else
            Ms = Msold;
        end
    end
    % Las repeticiones se almacenan una vez que pasó el periodo "burn in"
    if i > burnin
        if i == burnin + 1
            fprintf(1, '. Fin de las %s repeticiones "Burn in"\nAvance (%s repeticiones): ', num2str(burnin), num2str(R));
        end
        if rem(i, round(R/20)) == 0
            fprintf('%2.0f%% ', (i - burnin)/R*100)
        end
        MOD_(i - burnin,:)    = Ms;
        BETA_(i - burnin,Ms)  = Ebeta;
        VBETA_(i - burnin,Ms) = diag(Vbeta);
        logpM_(i - burnin,1)  = logPosts;
        SIGMA2_(i - burnin,1)  = S/n;
    else
        if i == 1
            fprintf(1, 'Avance (%s repeticiones burnin): ', num2str(burnin));
        else
            if rem(i, round(R/20)) == 0
                fprintf('%2.0f%% ', i/burnin*100)
            end
        end
    end
end
fprintf('\n')
fprintf('Tasa de aceptación (%%) = %s\n\n', num2str(Acc/iter*100))

%% Tamano del modelo: A priori y a posteriori
[n1,x1] = hist(binornd(K, 0.5, iter, 1), 1:40);
[n2,x2] = hist(sum(MOD_,2), 1:40);
f = figure; set(f, 'name', 'Distribución del tamaño del modelo', 'Units', 'normalized', 'Position', [30 15 60 60]/100)
hold on
b1 = bar(x1, 100*n1/sum(n1));
b2 = bar(x2, 100*n2/sum(n2));
hold off
set(b1, 'EdgeColor', 'r', 'FaceColor', 'r');
set(b2, 'EdgeColor', 'b', 'FaceColor', 'b', 'Barwidth', 0.6)
set(gca, 'Xtick', 2:25, 'Xlim', [1.5 25.5], 'Box', 'Off')
ylabel('100 \times Frecuencia')
leg = legend([b1 b2], '\itA priori', '\itA posteriori');
title('\bfDistribución del tamaño del modelo')
set(leg, 'Box', 'Off', 'Location', 'NorthEast')

%%
save _temp.mat BETA_ VBETA_ logpM_ SIGMA2_ MOD_
clear BETA_ VBETA_ logpM_ SIGMA2_

fprintf('Ordenando los modelos')
[MOD, I, J] = unique(MOD_, 'rows');
% MODELOS es una matriz de r x K que contiene los r "únicos" modelos
% visitados por la cadena. "i" es un índice tal que MODELOS = MODS(i, :).
% "j" es un índice tal que MOD = MODELOS(j,:).

clear MOD_
r     = length(I);
fprintf('. La cadena visitó %s modelos únicos.\n', num2str(r))
T     = tabulate(J);
frec  = T(:,2); % Número de veces que los modelos la matriz MODELOS han aparecido

load _temp
% Las realizaciones de los modelos únicos
BETA   = BETA_(I, :);
VBETA  = VBETA_(I, :);
SIGMA2 = SIGMA2_(I, :);

% Probabilidades a posteriori de los modelos
logpM = logpM_(I, :);

clear BETA_ VBETA_ logpM_ SIGMA2_ MOD_ T
delete('_temp.mat')

probA = exp( logpM - max(logpM) ); % Probabilidad "analítica"
probA = probA/sum(probA);

probN = frec/sum(frec);            % Probabilidad "numérica"

%% Resultado en figura
[~, ind] = sortrows(probA, -1);
indparamejor = ind;
pA = probA(ind);
pN = probN(ind);

f = figure; set(f, 'name', 'Probabilidad a posteriori de los 20 mejores modelos', 'Units', 'normalized', 'Position', [30 15 60 60]/100)
b = bar(100*[pA(1:20) pN(1:20)]);
set(b(1), 'EdgeColor', 'r', 'FaceColor', 'r')
set(b(2), 'EdgeColor', 'b', 'FaceColor', 'b')
set(gca, 'Xtick', 1:20, 'Xlim', [0.5 20.5], 'Box', 'Off')
ylabel('100 \times Probabilidad')
leg = legend(b, 'Analítica', 'Numérica');
set(leg, 'Box', 'Off', 'Location', 'NorthEast')
C = corrcoef([probN probA]);
title(['\bfProbabilidad a posteriori de los 20 mejores modelos (\rho = ' num2str(C(1,2)) ')'])

%% Resultados en pantalla
% Momentos a posteriori de beta
postprob = MOD'*probA;
postmean = BETA'*probA;

postvar = VBETA'*probA; + ( ( BETA - repmat(postmean', r, 1) ).^2 )'*probA;
poststd = sqrt(postvar);

% Ordenamos los datos (descendentemente) de acuerdo a la probabilidad de inclusión
[~, ind] = sort(-postprob(:));

linea = repmat('-', 1, 72);
fprintf('\n\nResultados a posteriori (promedio de modelos)\nComparar con Table 1 (p. 596) de Fernandez, Ley y Steel (2001) \n%s\n', linea)
fprintf('Variable                          \t Pr(a = 1) \t  Media \t Ratio t\n%s\n', linea)
for i = 1:K
    j = ind(i);
    if j >= K + 1
        continue
    end
    signo = ' ';
    if postmean(j) < 0
        signo = '-';
    end
    fprintf('%s \t   %4.3f \t %s%6.4f \t %6.4f \n', Xnames(j,:), postprob(j), signo, abs(postmean(j)), abs(postmean(j))/poststd(j))
end
fprintf('%s\n', linea)

for mej = 1:3
    mejor = indparamejor(mej);
    indmejor = MOD(mejor, :);
    linea = repmat('-', 1, 92);
    fprintf('\nResultados a posteriori (mejor modelo %1.0f)\n%s\n',  mej, linea)
    fprintf('Variable                          \t Pr(a = 1) \t  Media \t Ratio t \t  BMA m \t  BMA t\n%s\n', linea)
    for i = 1:K
        j = ind(i);
        if ~indmejor(j)
            continue
        end
        if j >= K + 1
            continue
        end
        signo = ' ';
        if BETA(mejor, j) < 0
            signo = '-';
        end
        fprintf('%s \t   %4.3f \t %s%6.4f \t %6.4f ', Xnames(j,:), postprob(j), signo, abs(BETA(mejor, j)), abs(BETA(mejor, j))/sqrt(VBETA(mejor, j)))
        signo = ' ';
        if postmean(j) < 0
            signo = '-';
        end
        fprintf('\t %s%6.4f \t %6.4f\n', signo, abs(postmean(j)), abs(postmean(j))/poststd(j))
    end
    fprintf('%s\n', linea)
end
