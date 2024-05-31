function [essentialregions,essentiallist,bootstats,essentialpvals] = Slidingwindow(taseq,bootstraps,windowsize,uniquenames,uniqueindices,threshold)
%first generate nbootstraps of size windowsize
%during this score all bootstraps for % TAs that are not empty and total
%tas hit
boots=randi(length(taseq),windowsize,bootstraps);
bootstats=zeros(2,bootstraps);
essentialregions=zeros(length(taseq),1);
essentiallist=zeros(length(uniquenames),1);
essentialpvals=zeros(length(taseq),1);
for i=1:bootstraps;
    taboots=taseq(boots(:,i),1);
    totalta=sum(taboots);
    numbernothit=length(taboots(taboots==0));
    numberhit=length(taboots(taboots~=0));
    fractionhit=numberhit/(numberhit+numbernothit);
    bootstats(1,i)=totalta;
    bootstats(2,i)=fractionhit;
end
i=1;
for i=1:length(taseq);
    windowval=0;
    if (length(taseq)-i) > windowsize;
        windowval=taseq(i:(i+windowsize)-1);
        sumwindowval=sum(windowval);
        [probofgettincurrentwindow] = ksdensity(bootstats(1,:),sumwindowval,'function','cdf');
        essentialpvals(i,1)=probofgettincurrentwindow;
        if probofgettincurrentwindow < threshold;
            essentialregions(i:i+windowsize,1)=1;
      
        end
    end
        
    if (length(taseq)-i) < windowsize;
        windowval=vertcat((taseq(i:i+(length(taseq))-i)),taseq((i-(length(taseq)-i)):(i-1)));
        sumwindowval=sum(windowval);
        [probofgettincurrentwindow] = ksdensity(bootstats(1,:),sumwindowval,'function','cdf');
         essentialpvals(i,1)=probofgettincurrentwindow;
        if probofgettincurrentwindow < threshold;
            essentialregions(i:length(taseq),1)=1;
        end
    end
end
i=1;
for i=1:length(uniquenames);
   essentialfeatures= essentialregions(uniqueindices(i,1):uniqueindices(i,2),1);
   if sum(essentialfeatures)==length(essentialfeatures);
       essentiallist(i,1)=1;
   end
end
