function [ bestx,bestv,O ] = rand1bin( fitfnnm,LB,UB,maxFEs,options )
D=options.D;
NP=options.NP;
F=options.F;
CR=options.CR;
maxStag=options.maxStag;
consfnnm=options.consfnnm;

FEs=0;
% --------------------initialization
X=zeros(D,NP);
for d=1:D
    X(d,:)=LB(d)+(UB(d)-LB(d)).*rand(1,NP);
end
% --------------------evaluation
Xfit=zeros(1,NP);
for i=1:NP
    Xfit(i)=feval(fitfnnm,X(:,i));
    FEs=FEs+1;
end

G=0;
terminate_cond=0;
while 1
    % --------------------terminate_cond check
    if FEs>maxFEs % maxFEs
        terminate_cond=1;
    elseif Stag>maxStag % maxStag
        terminate_cond=2;
    end
    
    if terminate_cond
       break; 
    end
    % --------------------mutation
    G=G+1;
    V=zeros(D,NP);
    for i=1:NP
        r=dfrandi(3,NP,i);
        V(:,i)=X(:,r(1))+F*(X(:,r(2))-X(:,r(3)));
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
    % --------------------update bestx & bestv
    [~,minidx]=min(Xfit);
    bestx=X(:,minidx);
    bestv=Xfit(minidx);

    bestvG(G)=bestv;
end
O.bestvG=bestvG;
O.suc_ctG=suc_ctG;
O.terminate_cond=terminate_cond;
end