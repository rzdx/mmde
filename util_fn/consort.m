function [ nX,nevX,nXidx ] = consort( sd,X,evX,minmax,e )
% recurrsive function 
% sd: current level
% X: pop.
% evX: pop. fitness
% minmax: 'min' or 'max' for each level
% e: error tolerance threshold (suggest:1e-8) (i.e. deem as 0: no violation occured)
stv=evX(sd,:);
if strcmp(minmax{sd},'min')
    [srtstv,idx]=sort(stv);
elseif strcmp(minmax{sd},'max')
    [srtstv,idx]=sort(stv,'descend');
else
    error('minmax_indetermined');
end
X=X(:,idx);
evX=evX(:,idx);

if sd<=1 || size(evX,2)==1
    nX=X;
    nevX=evX;
    nXidx=idx;
    return;
end

df=diff(srtstv);
df(df<e)=0;
seg=find(df);
seg=[0,seg,size(evX,2)];
nX=[];
nevX=[];
nXidx=[];
for i=1:length(seg)-1
    segX=X(:,seg(i)+1:seg(i+1));
    segevX=evX(:,seg(i)+1:seg(i+1));
    segidx=idx(:,seg(i)+1:seg(i+1));
    [tnX,tnevX,tnXidx]=consort(sd-1,segX,segevX,minmax,e);
    nX=[nX,tnX];
    nevX=[nevX,tnevX];
    nXidx=[nXidx,segidx(tnXidx)];
end
end