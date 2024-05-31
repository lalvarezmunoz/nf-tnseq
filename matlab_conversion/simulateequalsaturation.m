function [bootstrapcontrol]=simulateequalsaturation(taproportion,reads,inputreads,boots)
%taproportion is the gap in ta density between libraries
%reads is the number of reads to simulate.  Reads should equal the total
%number of reads in the experimental condition that you want to compare to
%input reads is the vector of control input
inputreadssum=sum(inputreads);
inputproportion=inputreads./inputreadssum;
inputproportiontanorm=inputproportion.*taproportion;
inputproportiontanorm(length(inputproportiontanorm)+1,1)=1-taproportion;
multinominputsample=mnrnd(reads,inputproportiontanorm,boots);
multinominputsample=multinominputsample';
multisum=sum(multinominputsample);
difference=multisum-multinominputsample(length(multinominputsample),1);
correction=repmat((reads./difference),length(multinominputsample),1);
correctedinput=multinominputsample.*correction;
bootstrapcontrol=correctedinput(1:length(inputreads),:);
end



