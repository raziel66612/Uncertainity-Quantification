function [y] = ishigami(xx,a,b)

x1=xx(1);
x2=xx(2);
x3=xx(3);

term1=sin(x1);
term2=a*(sin(x2))^2;
term3=b*(x3^4)*sin(x1);

y=term1 + term2 + term3;
end
