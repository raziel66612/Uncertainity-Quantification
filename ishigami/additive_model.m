function [y] = additive_model(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIVE FUNCTION
% Author: Danial Faghihi, University at Buffalo
%          
% INPUTS:
% xx = [theta1, theta2, theta3]
%
% OUTPUTS:
% Y = theta1 + theta2 + theta3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta1 = theta(1);
theta2 = theta(2);
theta3 = theta(3);
y = theta1 + theta2 + theta3;
end