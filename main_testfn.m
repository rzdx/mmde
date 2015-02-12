rng('default')
mpool
clear
clc
add_path

defdnm='de.';
defnnm='mmbest1bin';
defnnm2='best1bin';
evofnnm=[defdnm,defnnm];
evofnnm2=[defdnm,defnnm2];

XUBs=[3.14,6,5,1,1];
XLBs=[-3.14,0,-5,-1,-1];
YUBs=[3.14,8,5,1,1];
YLBs=[-3.14,2,-5,-1,-1];
Ds=[1,1,1,2,2];

NP=30;
F=0.7;
CR=0.9;
maxFEs=1e7;
maxStag=20;
TolX=0;
TolFv=0;

NP2=30;
F2=0.7;
CR2=0.9;
maxFEs2=NP2*200;
maxStag2=20;
TolX2=0;
TolFv2=0;

suc_thrd=1e-6;

runmax=2;

funmin=1;
funmax=5;

funlenmax=5;

tstfdnm='testfn.';
fitfnnms=cell(funlenmax,1);
confnnms=cell(funlenmax,1);
optsols=cell(funlenmax,1);
for i=funmin:funmax
    fitfnnms{i}=[tstfdnm,'p',num2str(i)];
    confnnms{i}=[tstfdnm,'p',num2str(i),'c'];
    optsols{i}=feval([tstfdnm,'optvalues'],i);
end

algo_para.NP=NP;
algo_para.F=F;
algo_para.CR=CR;
terminate_cond.maxFEs=maxFEs;
terminate_cond.maxStag=maxStag;
terminate_cond.TolX=TolX;
terminate_cond.TolFv=TolFv;

lv2_para.algo_para.NP=NP2;
lv2_para.algo_para.F=F2;
lv2_para.algo_para.CR=CR2;
lv2_para.terminate_cond.maxFEs=maxFEs2;
lv2_para.terminate_cond.maxStag=maxStag2;
lv2_para.terminate_cond.TolX=TolX2;
lv2_para.terminate_cond.TolFv=TolFv2;

dstrst=zeros(funlenmax,runmax);
Srst=zeros(funlenmax,runmax);
bestXY=cell(funlenmax,runmax);
bestV=zeros(funlenmax,runmax);
tic
for rnum=1:runmax
    disp(['rnum:',num2str(rnum)]);
    for fnnum=funmin:funmax
        disp(['fnnum:',num2str(fnnum)]);
        
        fnnm.fitfnnm=fitfnnms{fnnum};
        fnnm.consfnnm=confnnms{fnnum};
        fnnm.defnnm2=evofnnm2;
        bound.UB=XUBs(fnnum)*ones(Ds(fnnum),1);
        bound.LB=XLBs(fnnum)*ones(Ds(fnnum),1);
        lv2_para.bound.UB=YUBs(fnnum)*ones(Ds(fnnum),1);
        lv2_para.bound.LB=YLBs(fnnum)*ones(Ds(fnnum),1);
        
        [bestx,besty,bestFv,O]=feval(evofnnm,fnnm,bound,algo_para,terminate_cond,lv2_para);
        
        bestXY{fnnum,rnum}=[bestx,besty];
        bestV(fnnum,rnum)=bestFv;
        
        optdst=zeros(size(optsols{fnnum},2),1);
        for di=1:size(optsols{fnnum},2)
            optdst(di)=norm(optsols{fnnum}(:,di)-[bestx;besty]);
        end
        dstrst(fnnum,rnum)=min(optdst);
        if dstrst(fnnum,rnum)<=suc_thrd
            Srst(fnnum,rnum)=1;
        else
            Srst(fnnum,rnum)=0;
        end
    end
end
toc

rownm=cell(funlenmax,1);
colnm={'mean','std','success_rate'};
datav=zeros(length(rownm),length(colnm));
for dr=funmin:funmax
    datav(dr,1)=mean(dstrst(dr,:));
    datav(dr,2)=std(dstrst(dr,:),1);
    datav(dr,3)=sum(Srst(dr,:))/runmax;
    rownm{dr}=['kmiltering-',defnnm,'_p',num2str(dr)];
end
T={datav,rownm,colnm};
shT={T};

fldn=['tf_result-',datestr(now,'yyyymmddHHMMSS')];
sfldn={'meanstd','conv_gragh'};
fio.nfolds(fldn,sfldn);

save([fio.addslash(1,fldn,sfldn{1}),'T.mat'], 'T');
save([fio.addslash(1,fldn,sfldn{1}),'dstrst.mat'], 'dstrst');
save([fio.addslash(1,fldn,sfldn{1}),'bestXY.mat'], 'bestXY');
save([fio.addslash(1,fldn,sfldn{1}),'bestV.mat'], 'bestV');

XT={[fio.addslash(1,fldn,sfldn{1}),'mmde_rst.xls'],0,shT};
tio.xlswt(XT);

disp('---OVER---');