function [ g ] = p2c( x,y )
g=zeros(2,1);
g(1)=-(x-5)^2-(y-3)^2+4;
g(2)=(x-5)^2+(y-3)^2-16;
end