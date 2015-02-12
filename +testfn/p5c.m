function [ g ] = p5c( x,y )
g=zeros(4,1);
g(1)=x(1)-x(2);
g(2)=-x(1)+x(2);
g(3)=y(1)-y(2);
g(4)=-y(1)+y(2);
end