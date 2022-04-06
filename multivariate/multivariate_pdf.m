%---------------------------Question 4 ------------------------------------
mean=[0 0];
cov=[0.8 1;1 1.4];

x1=linspace(-10,10,250);
x2=linspace(-10,10,250);
[X1,X2] = meshgrid(x1, x2);
X=[X1(:) X2(:)];

pdf= mvnpdf(X,mean,cov); %generating probability distribution function
% size(pdf)

pdf=reshape(pdf,length(x2),length(x1));
% size(pdf)

figure
surf(x1,x2,pdf)
caxis([min(pdf(:))-1*range(pdf(:)),max(pdf(:))]);
axis([-3 3 -5 5 0 0.5])
xlabel('x1')
ylabel('x2')
zlabel('Probability Density(pdf)')
 
figure
plot(x1,pdf)  %marginal density
xlabel('x2')
ylabel('Marginal density')
grid on
grid minor

% plot of 500 random samples(x1 Vs x2)
pi_plot=mvnrnd(mean,cov,500);
figure
plot(pi_plot(:,1),pi_plot(:,2),'+r'); 
xlabel('x1')
ylabel('x2')
grid on
grid minor
axis([-4 4 -5 5])

%compute Kernal Density Function
[f,xi]= ksdensity(pi_plot(:,2));
figure
plot(xi,f);
ylabel('Density')
xlabel('x2')

