function [ g ] = p1c( x,y )
g=zeros(2,1);
g(1)=y-x*(x+6.28);
g(2)=y-x*(x-6.28);
end