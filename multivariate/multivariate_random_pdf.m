mean=[0 0];
cov=[0.8 1;1 1.4];
rng('default') % random number generator
X=mvnrnd(mean,cov,500);
% [X1,X2]=meshgrid(X(:,1),X(:,2));
% Y=[X1,X2];

pdf=mvnpdf(X,mean,cov);
% pdf=reshape([pdf],[length(X)],[length(X)]);

scatter3(X(:,1),X(:,2),pdf)

% x1=linspace(-10,10,350)
% x2=linspace(-10,10,350)
[x1,x2] = meshgrid(-10:0.5:10 , -10:0.5:10)
mesh(x1,x2,pdf)
