function [ bestx,bestFv,O ] = best1bin( fnnm,bound,algo_para,terminate_cond )
fitfnnm=fnnm.fitfnnm;
consfnnm=[];
LB=bound.LB;
UB=bound.UB;
NP=algo_para.NP;
F=algo_para.F;
CR=algo_para.CR;
maxFEs=terminate_cond.maxFEs;
maxStag=inf;
TolX=0;
TolFv=0;

if isfield(algo_para,'x')
    x=algo_para.x;
end
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
FEs=0;
% --------------------initialization
if isfield(algo_para,'X')
    X=algo_para.X;
else
    X=zeros(D,NP);
    for j=1:D
        X(j,:)=LB(j)+(UB(j)-LB(j))*rand(1,NP);
    end
end
% --------------------evaluation
Xfit=zeros(1,NP);
if exist('x','var')
    parfor i=1:NP
        Xfit(i)=feval(fitfnnm,x,X(:,i));
    end
else
    parfor i=1:NP
        Xfit(i)=feval(fitfnnm,X(:,i));
    end
end
FEs=FEs+NP;
% --------------------find best x
[~,minidx]=min(Xfit);
bestx=X(:,minidx);
bestFv=Xfit(minidx);
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
    % --------------------mutation
    G=G+1;
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
    if exist('x','var')
        parfor i=1:NP
            Ufit(:,i)=feval(fitfnnm,x,U(:,i));
        end
    else
        parfor i=1:NP
            Ufit(:,i)=feval(fitfnnm,U(:,i));
        end
    end
    FEs=FEs+NP;
    
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
O.bestFvG=bestFvG;
O.stdXG=stdXG;
O.suc_ctG=suc_ctG;
O.terminate_state=terminate_state;
end