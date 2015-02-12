rng('default')
mpool
clear
clc
add_path

defdnm='de.';
defnnm='best1bin'; %{rand1bin,best1bin}
evofnnm=[defdnm,defnnm];

kmO=kmproblem(1);
XUBs=kmO.ub;
XLBs=kmO.lb;
YUBs=XUBs;
YLBs=XLBs;
Ds=kmO.D;
NP=30;
F=0.7;
CR=0.9;
Gmax=200;
maxFEs=1e7;

runmax=2;

funmin=1;
funmax=1; 
funlenmax=1;

fitfnnms=cell(funlenmax,1);
confnnms=cell(funlenmax,1);
for i=funmin:funmax
    fitfnnms{i}='error_nozeta'; % D=[4,4] : [KA11,KA12;KA22,KC12] , [SA11,SA12;SA22,SC12]
    confnnms{i}='ConstraintViolation';
end

bestXY=cell(funlenmax,runmax);
bestV=zeros(funlenmax,runmax);
dstrst=zeros(funlenmax,runmax);

tic
for rnum=1:runmax
    disp(['rnum:',num2str(rnum)]);
    for fnnum=funmin:funmax
        disp(['fnnum:',num2str(fnnum)]);
        
        solverOptions1.dimensionFactor = 30;
        solverOptions1.F = F;
        solverOptions1.CR = CR;
        solverOptions1.TolX = 0;
        solverOptions1.TolFun = 0;
        solverOptions1.RecordPoint = 1000;
        solverOptions1.nonlcon = confnnms{fnnum};
        solverOptions1.innerMaxIter = 200;
        solverOptions1.migrateFactor = 0.7;
        
        solverOptions2.dimensionFactor = 30;
        solverOptions2.F = F;
        solverOptions2.CR = CR;
        solverOptions2.TolX = 0;
        solverOptions2.TolFun = 0;
        
        ub1 = XUBs;
        lb1 = XLBs;
        ub2 = ub1;
        lb2 = lb1;
        
        [bestx,besty,bestv,out]=mmdeb1b_pce(fitfnnms{fnnum},maxFEs,lb1,ub1,lb2,ub2,solverOptions1,solverOptions2);
        
        bestXY{fnnum,rnum}=[bestx,besty];
        bestV(fnnum,rnum)=bestv;
        
        dstrst(fnnum,rnum)=bestv;
    end
end
toc

rownm=cell(funlenmax,1);
colnm={'mean','std'};
datav=zeros(length(rownm),length(colnm));
for dr=funmin:funmax
    datav(dr,1)=mean(dstrst(dr,:));
    datav(dr,2)=std(dstrst(dr,:),1);
    rownm{dr}=['kmfiltering-',defnnm,'_p',num2str(dr)];
end
T={datav,rownm,colnm};
shT={T};

fldn=['km_result-',datestr(now, 'yyyymmddHHMMSS')];
sfldn={'meanstd','filter_gragh'};
fio.nfolds(fldn,sfldn);

save([fio.addslash(1,fldn,sfldn{1}),'T.mat'], 'T');
save([fio.addslash(1,fldn,sfldn{1}),'dstrst.mat'], 'dstrst');
save([fio.addslash(1,fldn,sfldn{1}),'bestXY.mat'], 'bestXY');
save([fio.addslash(1,fldn,sfldn{1}),'bestV.mat'], 'bestV');

XT={[fio.addslash(1,fldn,sfldn{1}),'kmfiltering_rst.xls'],0,shT};
tio.xlswt(XT);

disp('---OVER---');