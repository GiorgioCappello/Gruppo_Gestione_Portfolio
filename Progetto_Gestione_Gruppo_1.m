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

sigma = cov(returns);
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
% disp(sigma);
% 
% % Commento sulla matrice di varianza-covarianza storica e sulla matrice shrinkage
% disp('Matrice di varianza-covarianza storica:');
% disp(historical_covariance);
% 
% disp('Matrice di varianza-covarianza con shrinkage toward Constant Correlation Approach:');
% disp(shrinkage_covariance);


%% PUNTO 3
%% PUNTO 4
%% PUNTO 5
%% PUNTO 6
%% PUNTO 7
%% PUNTO 8
%% PUNTO 9

