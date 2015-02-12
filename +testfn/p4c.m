function [ g ] = p4c( x,y )
g=zeros(3,1);
g(1)=x(1)-x(2);
g(2)=y(1)-y(2);
g(3)=-y(1)+y(2);
end