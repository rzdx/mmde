function [ O ] = cafit( confnnm,x,y )
g=feval(confnnm,x,y);
O=abs(max(g));
end