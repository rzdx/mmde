function [ O ] = dfrandi( L,N,varargin )
rjset=[];
ci=0;
inii=1;
if ~isempty(varargin)
    if length(varargin{1})>=2
        for i=inii:length(varargin{1})
            ci=ci+1;
            rjset(ci)=varargin{1}(i);
        end
        inii=2;
    end
    if isempty(varargin{1})
        inii=2;
    end
end
for i=inii:length(varargin)
    ci=ci+1;
    rjset(ci)=varargin{i};
end

i=0;
while i<L
    i=i+1;
    t=floor(1+N*rand);
    if ~ismember(t,rjset)
        O(i)=t;
        rjset(length(rjset)+1)=t;
    else
        i=i-1;
    end
end
end

% L = number of random numbers
% N = max number of randi()
% varargin = list of never-appear elements
% O = 1~N exclude [varargin]