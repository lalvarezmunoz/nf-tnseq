function [sensitivityanalysis] =  checkseqdepth(taseqcontrol)
%
%
totalsequences=sum(taseqcontrol);
tahit=taseqcontrol(taseqcontrol~=0);
fractionhit=length(tahit)/length(taseqcontrol);
taprop=taseqcontrol./totalsequences;
y=1;
sensitivityanalysis=[];
for i=1:100000:totalsequences
    x=mnrnd(i,taprop,10);
    x=x';
    a=logical(x);
    x=sum(a);
    z=mean(x);
    sensitivityanalysis(y,1)=i;
    sensitivityanalysis(y,2)=z;
    y=y+1;
end

scatter(sensitivityanalysis(:,1),sensitivityanalysis(:,2))
end