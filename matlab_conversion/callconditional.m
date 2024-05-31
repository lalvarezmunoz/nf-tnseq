function [outputcond,callstats] = callconditional(stateconfidence,uniquenames,uniqueindices,minsl,maxsl,minenr,maxenr,propignore)
%Finds cond essential and enriched regions. 
outputcond=zeros(length(uniqueindices),1);
callstats=zeros(5,1);
for i=1:length(uniqueindices)
    annotcallsl=stateconfidence(uniqueindices(i,1):uniqueindices(i,2),4);
    annotcallenr=stateconfidence(uniqueindices(i,1):uniqueindices(i,2),3);
    len_ce = length(annotcallsl);
    len_enr = length(annotcallenr);
    x=round(len_ce*propignore);
    if x==0;
        x=1;
    end
    min_ce=min(annotcallsl(x:(len_ce-x),1));
    max_ce=max(annotcallsl(x:(len_ce-x),1));
    min_enr=min(annotcallenr(x:(len_enr-x),1));
    max_enr=max(annotcallenr(x:(len_enr-x),1));
    outputcond(i,1)=5;
    if max_ce >= maxsl;
        outputcond(i,1)=1;
        if min_ce >= minsl;
            outputcond(i,1)=2;
        end;
    end
    if max_enr >= maxenr;
        outputcond(i,1)=3;
        if min_enr >= minenr;
            outputcond(i,1)=4;
        end;
    end
    callstats(outputcond(i,1),1)=callstats(outputcond(i,1),1)+1;
    
end


