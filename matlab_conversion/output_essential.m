function [outputessential,callstats] = oputput_essential(finalessential,uniquenames,uniqueindices)
%Goes through and finds all annotations that have essential
%regions  You get a 1 for non essential, 2 for essential 
outputessential=zeros(length(uniqueindices),1);
callstats=zeros(3,1);
for i=1:length(uniqueindices);
    annotcallsess=finalessential(uniqueindices(i,1):uniqueindices(i,2),1);
    len_ce = length(annotcallsess);
    x=round(len_ce*0.20);
    if x==0;
        x=1;
    end
    min_ce=min(annotcallsess(x:(len_ce-x),1));
    max_ce=max(annotcallsess(x:(len_ce-x),1));
    outputessential(i,1)=1;
    if max_ce > 1;
        outputessential(i,1)=3;
        if min_ce > 1;
            outputessential(i,1)=2;
        end;
    end
    
    callstats(outputessential(i,1),1)=callstats(outputessential(i,1),1)+1;
    
end

