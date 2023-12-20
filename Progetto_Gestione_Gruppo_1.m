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

capitalization = [...
    table2array(stock_data_ADI(:,5)), ...
    table2array(stock_data_AMZN(:,5)), ...
    table2array(stock_data_APTV(:,5)), ...
    table2array(stock_data_ITW(:,5)), ...
    table2array(stock_data_PFE(:,5)), ...
    table2array(stock_data_SPG(:,5))
];

%% Gestione valori Nan

% Come suggerito in classe, tronchiamo i valori Nan presenti nella tabella stock_data_APTV
% Di conseguenza, per preservare la stessa dimensione, tronchiamo i valori delle altre tabelle
% anche se non sono NaN.
% Utilizziamo la funzione rmmissing che elimina le righe che presentano almeno un NaN

returns = rmmissing(returns);
capitalization = rmmissing(capitalization);

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

% % Commento sulle correlazioni
% disp('Matrice di correlazione tra le azioni:');
% disp(correlations);
% 
% % Commento sulla matrice di varianza-covarianza storica e sulla matrice shrinkage
% disp('Matrice di varianza-covarianza storica:');
% disp(sigma);
% 
% disp('Matrice di varianza-covarianza con shrinkage toward Constant Correlation Approach:');
% disp(sigma_SCC);


%% PUNTO 3

% Definizioni variabili
rf = 0.02/12;

% Calcolo Markovitz con no risk asset
A = (ones(number_stock,1)')*(sigma_SCC\mu_exp);
B = mu_exp'*(sigma_SCC\mu_exp);
C = ones(number_stock,1)'*(sigma_SCC\ones(number_stock,1));
H = B - 2*A*rf + C*rf^2;

% PROVA
D=B*C-A^2;

g=(B*(sigma_SCC\ones(number_stock,1))-A*(sigma_SCC\mu_exp))/D;
h=(C*(sigma_SCC\mu_exp)-A*(sigma_SCC\ones(number_stock,1)))/D;

w_min = sigma_SCC\ones(number_stock,1)/C;

sum(w_min)

m1=linspace(A/C,3,100);


Var_w=@(m)C/D*(m-A/C).^2+1/C;

% Calcolo portafoglio tangente
w_tang = (sigma_SCC\(mu_exp-rf*ones(number_stock,1))/(ones(number_stock,1)'*(sigma_SCC\(mu_exp-rf*ones(number_stock,1)))));

r_tang = w_tang'*mu_exp;
var_tang = w_tang'*sigma_SCC*w_tang;

% Calcolo frontiera efficiente
m = linspace(rf,3,100);

var_w = @(x) (x - rf).^2 / H;

% Plot
figure
hold on
plot(sqrt(var_w(m)), m)
plot(sqrt(Var_w(m1)),m1)
scatter(sqrt(var_tang), r_tang)

%% PUNTO 4
%% PUNTO 5
%% PUNTO 6
%% PUNTO 7
%% PUNTO 8
%% PUNTO 9

