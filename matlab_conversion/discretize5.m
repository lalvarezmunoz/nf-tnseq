function [discta] = discretize5(TAarray,zero,one,sevfive,whisker)
for i=1:length(TAarray)
    if TAarray(i,1)==zero
        discta(i,1)=1;
    end
    if TAarray(i,1)>0
        if TAarray(i,1)<=one
        discta(i,1)=2;
        end
    end
    if TAarray(i,1)>one
        if TAarray(i,1)<sevfive
        discta(i,1)=3;
        end
    end
    if TAarray(i,1)>(sevfive-1)
        if TAarray(i,1)<whisker
        discta(i,1)=4;
        end
    end
    if TAarray(i,1)>(whisker-1)
        discta(i,1)=5;
    end
    
end
end

    

