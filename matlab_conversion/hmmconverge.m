function [outmat,a_0,b_0] = hmmconverge(genometa,initialtrans,initialemission)
i=1;
[m,n]=size(genometa);
outmat=zeros(m,100);
a=initialtrans;
b=initialemission;
j=1;
while i > 0
    a_0=a;
    b_0=b;
    [states]=hmmviterbi(genometa,a_0,b_0);
    [a,b]=hmmestimate(genometa,states);
    outmat(:,j)=states;
    i=(a(1,1)-a_0(1,1));
    j=j+1;
end

