%% MATLAB script mle_example.m
% MLE estimation for the elasticity problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
clear all; close all; clc

%% generating synthetic data based on the true parameter 
% number of data points	
n = 10;   
% true Young's modulus(MPa)                      
E_true = 200*1e3; 
% stress (MPa)              
T_data = linspace(50, 330, n); 
% strain 
e_data = T_data./(E_true*1e3) + 0.002.*(T_data./250).^0.5;  % %Strain Case 1
% e_data = [0.0009,0.0011,0.0013,0.0015,0.0017,0.0018,0.0019,0.0021,0.0022,0.0023];   %Strain Case 1
% e_data = [0.0009,0.0012,0.0014,0.0016,0.0018,0.0020,0.0022,0.0024,0.0025,0.0027];     %Strain Case 2  

%% compute MLE and Fisher information matrix
Ehat = (sum(T_data))/(sum(e_data))
sigmahat = sqrt((1/n) * sum((T_data/e_data - Ehat).^2))

% In = [ sigmahat^(-2), 0
%     0, n*sigmahat^(-2)-sum(3.*(T_data/e_data - Ehat).^2) * sigmahat^(-4)];

In = [ 1*sigmahat^(-2), 0
    0, 2*n*sigmahat^(-2)];

% inverse of Fisher information
InvIn = inv(In);

% estimation with plus and minus stDev
Ehat_up = Ehat + 1.645*sqrt(InvIn(1,1))
Ehat_lo = Ehat - 1.645*sqrt(InvIn(1,1))

%% plot data and model prediction
e_model_hat = T_data./(Ehat);
e_model_up  = T_data./(Ehat_up);
e_model_lo  = T_data./(Ehat_lo);

figure
plot(e_data, T_data, 'ob', 'MarkerSize',10)
hold on
plot(e_model_up, T_data, '--b', 'linewidth', 3)
plot(e_model_hat, T_data, 'k', 'linewidth', 3)
plot(e_model_lo, T_data, '--r', 'linewidth', 3)
xlabel('Strain'); ylabel('Stress');
ylim([0,350])
legend('data', 'upper bound', 'MLE', 'lower bound')
set(gca,'FontSize',24)

InvIn