%% MATLAB script for
% scatterplots of the additive model
%	Y = c1X1 + c2X2		
clear all; close all; clc	

%  coefficients
c1 = 2;	
c2 = 1;	
sigma1 = 1;
sigma2 = 3;

% number od samples
N = 100000;

% drawing samples from parameters and construct sample matrix M
x1 = normrnd(0,sigma1,[N,1]);
x2 = normrnd(0,sigma2,[N,1]);	

M = [x1, x2];

% compute vector of model outputs
Y = c1*x1 + c2*x2;

%% scater plots
figure
plot(M(:,1), Y, '*b')
xlim([-4 4]), grid on
axis square,xlabel('x_1'),ylabel('y')		
set(gca,'FontSize',24)	
print('x1y','-dpng')	
	
figure
plot(M(:,2), Y, '*b')
xlim([-10 10]), grid on		
axis square,xlabel('x_2'),ylabel('y')		
set(gca,'FontSize',24)	
print('x2y','-dpng')
