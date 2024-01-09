clc
clear all
close all

%% Caricamento dei dati

stock_data_ADI = read_stock("ADI.OQ.csv");
stock_data_AMZN = read_stock("AMZN.OQ.csv");
stock_data_APTV = read_stock("APTV.N.csv");
stock_data_ITW = read_stock("ITW.N.csv");
stock_data_PFE = read_stock("PFE.N.csv");
stock_data_SPG = read_stock("SPG.N.csv");
stock_data_SPX = read_stock("SPX.csv");

%% Creazione meta dati

stock_names = ["ADI.OQ";"AMZN.OQ";"APTV.N";"ITW.N";"PFE.N";"SPG.N"];
number_stock = 6;

returns = [...
    table2array(stock_data_ADI(:, 4)), ...
    table2array(stock_data_AMZN(:, 4)), ...
    table2array(stock_data_APTV(:, 4)), ...
    table2array(stock_data_ITW(:, 4)), ...
    table2array(stock_data_PFE(:, 4)), ...
    table2array(stock_data_SPG(:, 4)) 
];

market_prices = [...
    table2array(stock_data_SPX(:, 4))];

capitalization = [...
    table2array(stock_data_ADI(:,5)), ...
    table2array(stock_data_AMZN(:,5)), ...
    table2array(stock_data_APTV(:,5)), ...
    table2array(stock_data_ITW(:,5)), ...
    table2array(stock_data_PFE(:,5)), ...
    table2array(stock_data_SPG(:,5))
];

prices = [...
    table2array(stock_data_ADI(:,6)), ...
    table2array(stock_data_AMZN(:,6)), ...
    table2array(stock_data_APTV(:,6)), ...
    table2array(stock_data_ITW(:,6)), ...
    table2array(stock_data_PFE(:,6)), ...
    table2array(stock_data_SPG(:,6))
];

%% Gestione valori Nan

% Come suggerito in classe, tronchiamo i valori Nan presenti nella tabella stock_data_APTV
% Di conseguenza, per preservare la stessa dimensione, tronchiamo i valori delle altre tabelle
% anche se non sono NaN.
% Utilizziamo la funzione rmmissing che elimina le righe che presentano almeno un NaN

returns = rmmissing(returns);
prices = rmmissing(prices);
market_prices  = rmmissing(market_prices);
capitalization = rmmissing(capitalization);

returns = returns./100;

%% PUNTO 2

mu = mean(returns)';

% Shrinkage toward Constant Correlation Approach

k = 0.3;

% Calcolo matrice varianza-covarianza
sigma = cov(returns);
% standard_dev
vola = sqrt(diag(sigma));

correlations = corr(returns);
rho = (sum(sum(correlations))-number_stock) / (number_stock*(number_stock-1));

correlations_CC = ones(size(correlations))*rho + (1-rho)*diag(ones(number_stock,1));
sigma_CC = diag(vola)*correlations_CC*diag(vola);
sigma_SCC = (1-k)*sigma + k*sigma_CC;

% Stima del vettore delle medie con una media esponenziale
lambda = 0.01;

T = size(returns,1);
t = 1:T;
mu_exp = sum(returns.*(exp(-lambda*(T-t)))')'/sum(exp(-lambda*(T-t)));

% Commento sulle correlazioni
disp('Matrice di correlazione tra le azioni:');
disp(correlations);

% Commento sulla matrice di varianza-covarianza storica e sulla matrice shrinkage
disp('Matrice di varianza-covarianza storica:');
disp(sigma);

disp('Matrice di varianza-covarianza con shrinkage toward Constant Correlation Approach:');
disp(sigma_SCC);


%% PUNTO 3

% Definizioni variabili
rf = 0.02/12;

% Calcolo frontiera efficiente con risk-free asset
A = (ones(number_stock,1)')*(sigma_SCC\mu_exp);
B = mu_exp'*(sigma_SCC\mu_exp);
C = ones(number_stock,1)'*(sigma_SCC\ones(number_stock,1));
H = B - 2*A*rf + C*rf^2;

m = linspace(rf,0.05,100);

var_rf = @(x) (x - rf).^2 / H;

% Calcolo frontiera efficiente senza risk-free asset
D=B*C-A^2;

g=(B*(sigma_SCC\ones(number_stock,1))-A*(sigma_SCC\mu_exp))/D;
h=(C*(sigma_SCC\mu_exp)-A*(sigma_SCC\ones(number_stock,1)))/D;

w_min = sigma_SCC\ones(number_stock,1)/C;

m1=linspace(A/C,0.05,100);

var_no_rf=@(m)C/D*(m-A/C).^2+1/C;

% Calcolo portafoglio tangente
w_tang = (sigma_SCC\(mu_exp-rf*ones(number_stock,1))/(ones(number_stock,1)'*(sigma_SCC\(mu_exp-rf*ones(number_stock,1)))));

r_tang = w_tang'*mu_exp;
var_tang = w_tang'*sigma_SCC*w_tang;

% Plot
figure
hold on
grid on
plot(sqrt(var_rf(m)), m, 'LineWidth', 2, 'DisplayName', 'Con Risk-Free Asset')
scatter(0, rf, 50, 'blue', 'filled','DisplayName', 'Risk-Free ')
plot(sqrt(var_no_rf(m1)), m1, 'LineWidth', 2, 'DisplayName', 'Senza Risk-Free Asset', 'Color', 'red')
scatter(sqrt(var_no_rf(A/C)), A/C , 50, 'red', 'filled','DisplayName', 'A/C')
scatter(sqrt(var_tang), r_tang, 50, [150, 50, 200]/255, 'filled', 'DisplayName', 'Portafoglio Tangente')
xlabel('Deviazione Standard')
ylabel('Rendimento Atteso')
title('Frontiere e Portafoglio Tangente')
text(0.02, rf, ' rf ', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'blue');
text(sqrt(var_no_rf(A/C)), A/C, ' A/C', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', 'red');
text(sqrt(var_tang), r_tang, '    Portafoglio Tangente', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', [150, 50, 200]/255);
legend('Location', 'Best')
hold off

%% PUNTO 4

% Punto a: plot delle tre frontiere

x0 = w_min;
m2 = linspace(rf,0.01,100);

% Vincolo 1
Aeq1 = [(mu_exp - rf)'; 1 1 0 0 0 0];

% Ciclo per trovare la frontiera
var_vincolo1_rf = zeros(length(m2), 1);

for i = 1:length(m2)
    beq1 = [m2(i) - rf; 0.5];
    w_1 = (m2(i) - rf) / H * (sigma_SCC \ (mu_exp - rf));
    ww_1 = fmincon(@(ww) ww' * sigma_SCC * ww, x0, [], [], Aeq1, beq1);
    var_vincolo1_rf(i) = ww_1' * sigma_SCC * ww_1;
end

% Vincolo 2
Aeq2 = (mu_exp - rf)';
A2 = - eye(number_stock);
b2 = - 0.1 * ones(number_stock, 1);

% Ciclo per trovare la frontiera
var_vincolo2_rf = zeros(length(m2), 1);

for i = 1:length(m2)
    beq2 = m2(i) - rf;
    w_2 = (m2(i) - rf) / H * (sigma_SCC \ (mu_exp - rf));
    ww_2 = fmincon(@(ww) ww' * sigma_SCC * ww, x0, A2, b2, Aeq2, beq2);
    var_vincolo2_rf(i) = ww_2' * sigma_SCC * ww_2;
end

% Punto b: calcolo dei portafogli con rendimento atteso del 0.5%
target_return = 0.005;

% No Vincolo
Aeq_no = (mu_exp - rf)';
beq_no = target_return - rf;
w_no_target_return = fmincon(@(w) w' * sigma_SCC * w, x0, [], [], Aeq_no, beq_no);

% Vincolo 1
beq1 = [target_return - rf; 0.5];
w_vincolo1_target_return = fmincon(@(ww) ww' * sigma_SCC * ww, x0, [], [], Aeq1, beq1);

% Vincolo 2
beq2 = target_return - rf;
w_vincolo2_target_return = fmincon(@(w) w' * sigma_SCC * w, x0, A2, b2, Aeq2, beq2);

% Calcolo delle varianze dei portafogli ottenuti
var_no_target_return = w_no_target_return' * sigma_SCC * w_no_target_return;
var_vincolo1_target_return = w_vincolo1_target_return' * sigma_SCC * w_vincolo1_target_return;
var_vincolo2_target_return = w_vincolo2_target_return' * sigma_SCC * w_vincolo2_target_return;

%% Plot dei risultati
figure
hold on
grid on
plot(sqrt(var_rf(m2)), m2, 'b', 'LineWidth', 2, 'DisplayName', 'Frontiera no vincoli')
scatter(sqrt(var_no_target_return), target_return, 100, 'b', 'filled', 'DisplayName', 'Target no vincoli')
plot(sqrt(var_vincolo1_rf), m2, 'r', 'LineWidth', 2, 'DisplayName', 'Frontiera vincolo 1')
scatter(sqrt(var_vincolo1_target_return), target_return, 100, 'r', 'filled', 'DisplayName', 'Target vincolo 1')
plot(sqrt(var_vincolo2_rf), m2, 'g', 'LineWidth', 2, 'DisplayName', 'Frontiera vincolo 2')
scatter(sqrt(var_vincolo2_target_return), target_return, 100, 'g', 'filled', 'DisplayName', 'Target vincolo 2')
text(sqrt(var_no_target_return), target_return, ' No Vincoli', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'b');
text(sqrt(var_vincolo1_target_return), target_return, ' Vincolo 1', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'r');
text(sqrt(var_vincolo2_target_return), target_return, ' Vincolo 2', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'g');
legend('Location', 'Best')
hold off

%% PUNTO 5

% riduco i market prices per accoppiare le dimensioni;
market_prices = market_prices(95:239);

% calcolo i market returns;
market_returns = zeros(1,length(market_prices));
for i=2:length(market_prices)
    market_returns(i) = (market_prices(i)-market_prices(i-1))./(market_prices(i-1));
end

% riduzione dimensione per adattamento lm
returns_CAPM = returns(2:145,:)*100;
market_returns_CAPM=market_returns(2:145)*100;
rf_CAPM = rf*100;

% Computo beta e alpha per ciascuno stock procedendo in ordine alfabetico
alpha=zeros(6,1);
beta=zeros(6,1);
pi_alpha=zeros(6,1);
pi_beta=zeros(6,1);
pi_test=ones(6,1);
for i=1:length(stock_names)
    lm = fitlm(market_returns_CAPM-rf_CAPM,returns_CAPM(:,i)-rf_CAPM); 
    LM = table2array(lm.Coefficients);
    alpha(i)=LM(1,1);
    beta(i)=LM(1,2);
    pi_alpha(i)=LM(1,4);
    pi_beta(i)=LM(2,4);
    pi_test(i)=coefTest(lm);
end

% riporto in tabella per pulizia grafica
LM_table = array2table([alpha,pi_alpha,beta,pi_beta,pi_test],'VariableNames',{'Alpha','pvalue di Alpha','Beta','pvalue di Beta','pvalue del modello'},'RowNames',stock_names);

%% PUNTO 6 

% calcolo la total market capitalization
total_market_cap=sum(sum(capitalization));

% e i pesi del portafoglio di mercato
w_market=zeros(6,1);
for i=1:length(stock_names)
    w_market(i)=sum(capitalization(:,i))/total_market_cap;
end

% e gli implicit market excess returns (PI)
PI = sigma_SCC*w_market;

%% PUNTO 7

% calcolo dei rendimenti attesi a Gennaio con metodo Monte Carlo 'trained'
% sulla serie storica dei rendimenti percentuali
jan_returns=zeros(6,1);
simul_jan_returns=zeros(1000,6);
mean_hist = mean(returns)';
standev_hist = std(returns)';
max_hist = max(returns);
min_hist = min(returns);
for j=1:length(stock_names)
    for i=1:1000
        simul_jan_returns(i,j) = mean_hist(j) + standev_hist(j)*((max_hist(j)-min_hist(j))*rand + min_hist(j));
    end
    jan_returns(j)=mean(simul_jan_returns(:,j));
end

% calcolo E[r] e il portafoglio con una prima view determinata tramite il
% metodo Monte Carlo sovrastante, scegliendo il primo e quarto titolo con
% una certainty 

% VIEW 1

% Amazon batterà Simon Property Group della media dei valori stimati con montecarlo
Q1 = (jan_returns(2)-jan_returns(6)) / 2;
Certainty1 = 0.3; 
P1 = [0, 1, 0, 0, 0, -1];

[PI_new1, sigma_new1, w_new1] = calculateView(sigma_SCC, PI, P1, Q1, Certainty1);

vola_new1 = sqrt(w_new1'*sigma_new1*w_new1);
r_new1 = w_new1'*PI_new1;

% VIEW 2

% Pfizer sale seguendo il momentum e le news disponibili al 23 dicembre
Q2 = 0.01; 
Certainty2 = 0.40; 
P2 = [0, 0, 0, 0, 1, 0];

[PI_new2, sigma_new2, w_new2] = calculateView(sigma_SCC, PI, P2, Q2, Certainty2);

vola_new2 = sqrt(w_new2'*sigma_new2*w_new2);
r_new2 = w_new2'*PI_new2;

% VIEW UNITA
Q = [(jan_returns(2)-jan_returns(6)) / 2; 0.01];
Certainty = [0.3 ; 0.4]; 
P = [0, 1, 0, 0, 0, -1;
     0, 0, 0, 0, 1, 0];

[PI_new, sigma_new, w_new] = calculateView(sigma_SCC, PI, P, Q, Certainty);

vola_new = sqrt(w_new'*sigma_new*w_new);
r_new = w_new'*PI_new;
