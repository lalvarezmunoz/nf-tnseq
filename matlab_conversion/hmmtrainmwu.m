function [hmmmatrix,stateconfidence,trainingmatrix] = hmmtrainmwu(experiment,controlbootstraps,trainingstates,uniquenames,uniqueindices,mwuthr,denominator1,denominator2,denominator3,q1,q2,q3,q4)
%Full function to train and export hmm matrices
%version _dc3 denotes using the conditional discretization version three
%that has far more parameters.
%   Detailed explanation goes here
%trainingstates referes to the originalessential versus non essentialcalls
[~,n]=size(controlbootstraps);
hmmmatrix=zeros(length(experiment),n);
trainingmatrix=zeros(length(experiment),n);
for i=1:n;
    %step1 is to discretize
    [discta] = discretizeconditional_3(experiment,controlbootstraps(:,i),denominator1,denominator2,denominator3,q1,q2,q3,q4);
    mwutrain=trainingstates;
    aucmwutrain=trainingstates;
    for j=1:length(uniquenames);
        ratio=sum(experiment(uniqueindices(j,1):uniqueindices(j,2),1))/sum(controlbootstraps(uniqueindices(j,1):uniqueindices(j,2),i));
        labels=[];
        scores1=[];
        scores2=[];
        labels(1:((uniqueindices(j,2)-uniqueindices(j,1))+1),1)=0;
        labels(((uniqueindices(j,2)-uniqueindices(j,1))+2):(2*((uniqueindices(j,2)-uniqueindices(j,1))+1)),1)=1;
        scores1(1:((uniqueindices(j,2)-uniqueindices(j,1))+1),1)=experiment(uniqueindices(j,1):uniqueindices(j,2),1);
        scores2(1:((uniqueindices(j,2)-uniqueindices(j,1))+1),1)=controlbootstraps(uniqueindices(j,1):uniqueindices(j,2),1);
        [p,h] = ranksum(scores1,scores2);
        if ratio > 1;
            if p<mwuthr;
                mwutrain(uniqueindices(j,1):uniqueindices(j,2),1)=3;
            end
        end
        if ratio < 1;
            if p<mwuthr;
                mwutrain(uniqueindices(j,1):uniqueindices(j,2),1)=4;
            end
        end
    end
    [truse,em]=hmmestimate(discta,mwutrain);
    hmm=hmmviterbi(discta,truse,em);
    hmmmatrix(:,i)=hmm;
    trainingmatrix(:,i)=mwutrain;
end
j=1;
stateconfidence=zeros(length(experiment),max(hmm));
for j=1:length(experiment);
    for k=1:n;
        stateconfidence(j,hmmmatrix(j,k))=stateconfidence(j,hmmmatrix(j,k))+1;
    end
end
stateconfidence=stateconfidence./n;
end
