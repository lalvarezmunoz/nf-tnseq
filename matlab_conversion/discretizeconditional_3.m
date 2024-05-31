
function [discta] = discretizeconditional_3(Experiment,Control,denominator1,denominator2,denominator3,q1,q2,q3,q4)
discta=zeros(length(Control),1);
for i=1:length(Experiment);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%No Reads%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if Experiment(i,1)==0;
        if Control(i,1)==0;
            discta(i,1)=1;
        end
  end
  %%%%%%%%%%%%%%%%%%%%%%Quantitative differences%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if log2(Experiment(i,1)/Control(i,1))>0;
        if log2(Experiment(i,1)/Control(i,1))<=q3;
        discta(i,1)=2;
        end
  end
  if log2(Experiment(i,1)/Control(i,1))>q3;
        if log2(Experiment(i,1)/Control(i,1))<=q4;
        discta(i,1)=3;
        end
  end
  if log2(Experiment(i,1)/Control(i,1))>=q4;
      discta(i,1)=4;
  end
  if log2(Experiment(i,1)/Control(i,1))<0;
        if log2(Experiment(i,1)/Control(i,1))>q1;
        discta(i,1)=5;
        end
  end
  if log2(Experiment(i,1)/Control(i,1))<=q1;
        if log2(Experiment(i,1)/Control(i,1))>q2;
        discta(i,1)=6;
        end
  end
  if log2(Experiment(i,1)/Control(i,1))<=q2;
        discta(i,1)=7;
  end
  %%%%%%%%%%%%%%%%%Experiment not equal to zero%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if Experiment(i,1)>0;
      if Experiment(i,1)<=denominator1;
        if Control(i,1)==0;
           discta(i,1)=8;
        end
      end
  end
  if Experiment(i,1)>denominator1;
      if Experiment(i,1)<=denominator2;
        if Control(i,1)==0;
           discta(i,1)=9;
        end
      end
  end
  if Experiment(i,1)>denominator2;
      if Experiment(i,1)<=denominator3;
        if Control(i,1)==0;
           discta(i,1)=10;
        end
      end
  end
  if Experiment(i,1)>denominator3;
        if Control(i,1)==0;
           discta(i,1)=11;
        end
  end
  
  %%%%%%%Control not equal to zero Experiment equal to zero%%%%%%%%%%%%%%%%
   if Control(i,1)>0;
      if Control(i,1)<=denominator1;
        if Experiment(i,1)==0;
           discta(i,1)=12;
        end
      end
  end
  if Control(i,1)>denominator1;
      if Control(i,1)<=denominator2;
        if Experiment(i,1)==0;
           discta(i,1)=13;
        end
      end
  end
  if Control(i,1)>denominator2;
      if Control(i,1)<=denominator3;
        if Experiment(i,1)==0
           discta(i,1)=14;
        end
      end
  end
  if Control(i,1)>denominator3;
        if Experiment(i,1)==0;
           discta(i,1)=15;
        end
  end
end
  
end  

