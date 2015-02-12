function [ bestx,bestFv,O ] = mmbest1bin( fnnm,bound,algo_para,terminate_cond,lv2_para )
Xminmax={'min','min'};
Yminmax={'max','min'};
e=0;
fitfnnm=fnnm.fitfnnm;
consfnnm=[];
cafnnm=fnnm.cafnnm;
defnnm2=fnnm.defnnm2;
LB=bound.LB;
UB=bound.UB;
NP=algo_para.NP;
F=algo_para.F;
CR=algo_para.CR;
maxFEs=terminate_cond.maxFEs;
maxStag=inf;
TolX=0;
TolFv=0;

NP2=lv2_para.algo_para.NP;

if isfield(fnnm,'consfnnm')
    consfnnm=fnnm.consfnnm;
end
if isfield(terminate_cond,'maxStag')
    maxStag=terminate_cond.maxStag;
end
if isfield(terminate_cond,'TolX')
    TolX=terminate_cond.TolX;
end
if isfield(terminate_cond,'TolFv')
    TolFv=terminate_cond.TolFv;
end

D=size(LB,1);
caxArch=zeros(D,NP);
FEs=0;
% --------------------initialization
if isfield(algo_para,'X')
    X=algo_para.X;
    Y=algo_para.Y;
else
    X=zeros(D,NP);
    Y=zeros(D,NP2,NP);
    for j=1:D
        X(j,:)=LB(j)+(UB(j)-LB(j))*rand(1,NP);
    end
    for i=1:NP
        for j=1:D
            Y(j,:,i)=LB(j)+(UB(j)-LB(j))*rand(1,NP,1);
        end
    end
end
% --------------------constraint active
for i=1:NP
    lv2_para.algo_para.x=X(:,i);
    [cax,~,O]=feval(defnnm2,cafnnm,lv2_para.bound,lv2_para.algo_para,lv2_para.terminate_cond);
    FEs=FEs+O.FEs;
    caxArch(:,i)=cax;
end
% --------------------evaluation function value (Y)
Yfit=zeros(1,NP2,NP);
for i=1:NP
    x=X(:,i);
    parfor k=1:NP2
        Yfit(1,k,i)=feval(fitfnnm,x,Y(:,k,i));
        FEs=FEs+1;
    end
end
% --------------------evaluation constraint violation (Y)
Ycfit=zeros(1,NP2,NP);
for i=1:NP
    x=X(:,i);
    parfor k=1:NP2
        g=feval(consfnnm,x,Y(:,k,i));
        Ycfit(1,k,i)=sum(g(g>0));
        FEs=FEs+1;
    end
end
% --------------------sort (Y)
for i=1:NP
    srtYfit=[Yfit(1,:,i);Ycfit(1,:,i)];
    srtY=Y(:,:,i);
    [srtY,srtYfit]=consort(numel(Yminmax),srtY,srtYfit,Yminmax,e);
    Y(:,:,i)=srtY;
    Yfit(1,:,i)=srtYfit(1,:,i);
    Ycfit(1,:,i)=srtYfit(2,:,i);
end
% --------------------sort (X)
Xfit=Yfit(1,1,:);
Xcfit=Ycfit(1,1,:);
srtXfit=[Xfit;Xcfit];
[X,srtXfit,srtXidx]=consort(numel(Xminmax),X,srtXfit,Xminmax,e);
Xfit=srtXfit(1,:);
Xcfit=srtXfit(2,:);
Y=Y(:,:,srtXidx);
Yfit=Yfit(:,:,srtXidx);
Ycfit=Ycfit(:,:,srtXidx);
% --------------------get bestx
bestx=X(:,1);
bestFv=Xfit(1);
% --------------------vertical selection // applicable only when NP=NP2
if NP==NP2
    vYfit=zeros(1,NP,NP);
    vYcfit=zeros(1,NP,NP);
    for i=1:NP
        x=X(:,i);
        parfor k=1:NP
            if i==k
                continue;
            end
            vYfit(1,i,k)=feval(fitfnnm,x,Y(:,i,k));
            g=feval(consfnnm,x,Y(:,i,k));
            vYcfit(1,i,k)=sum(g(g>0));
            FEs=FEs+2;
        end
    end
else
    disp('---verical selection not used---');
end
% --------------------horizontal evolution
for i=1:NP
    x=X(:,i);
    parfor k=1:NP2
        feval();
    end
end
% --------------------loop start
stdX=std(X,0,2);
G=0;
terminate_state=0;
while 1
    % --------------------terminate_state check
    if FEs>maxFEs
        terminate_state=1;
    end
    if Stags>maxStag
        terminate_state=2;
    end
    if stdX<TolX
        terminate_state=3;
    end
    if std(Xfit,0,2)<TolFv
        terminate_state=4;
    end
    % --------------------termination
    if terminate_state
        break;
    end
    G=G+1;
    % --------------------mutation
    V=zeros(D,NP);
    for i=1:NP
        r=dfrandi(2,NP,i);
        V(:,i)=bestx+F*(X(:,r(1))-X(:,r(2)));
    end
    % --------------------crossover
    U=zeros(D,NP);
    for i=1:NP
        randj=randi(D);
        for j=1:D
            if rand<=CR || randj==j
                U(j,i)=V(j,i);
            else
                U(j,i)=X(j,i);
            end
        end
    end
    % --------------------out of bound handling: bounce back then align
    for i=1:NP
        for j=1:D
            if U(j,i)<LB(j)
                U(j,i)=2*LB(j)-U(j,i);
            end
            if U(j,i)>UB(j)
                U(j,i)=2*UB(j)-U(j,i);
            end
            if U(j,i)<LB(j)
                U(j,i)=LB(j);
            end
            if U(j,i)>UB(j)
                U(j,i)=UB(j);
            end
        end
    end
    % --------------------selection
    Ufit=zeros(1,NP);
    for i=1:NP
        Ufit(:,i)=feval(fitfnnm,U(:,i));
        FEs=FEs+1;
    end
    
    for i=1:NP
        if Ufit(i)<=Xfit(i)
            X(:,i)=U(:,i);
            Xfit(i)=Ufit(i);
            suc_ctG(G)=suc_ctG(G)+1;
        end
    end
    % --------------------update bestx & bestFv
    [~,minidx]=min(Xfit);
    bestx=X(:,minidx);
    bestFv=Xfit(minidx);
    stdX=std(X,0,2);
    
    bestFvG(G)=bestFv;
    stdXG(G)=stdX;
end
% --------------------loop end
O.bestFvG=bestFvG;
O.stdXG=stdXG;
O.suc_ctG=suc_ctG;
O.terminate_state=terminate_state;
end