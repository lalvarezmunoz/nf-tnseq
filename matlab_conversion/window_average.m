function [normalized_reads] = window_average(TAsites, Total_reads, windowsize, genome_length)
%% Sets up the genomic window that you want with the first window starting at the first nt

num_windows=genome_length/windowsize; %calculates how many windows you will create
windowfix=fix(num_windows); %rounds down the number of window
leftover_window_start=((windowfix)*windowsize)+1;
leftover_window_end=((windowfix)*windowsize) + (genome_length-(windowfix*windowsize)); %uses the reamining reads to define

nonzeroindices=find(Total_reads>0); % finds where sites with reads are in the list from TAsites or Tn5sites
TAsites_nonzero=TAsites(nonzeroindices);
%zeroindices=find(Total_reads==0); % find where the sites without reads are in the list
%TAsites_zero=TAsites(zerodindices);
reads_nonzero=Total_reads(nonzeroindices,1); % pulls out the reads at each insertion when they are >0) 

chr_ave=(sum(Total_reads))/length(reads_nonzero); % caulcuates chromosomal read average at all sites with reads

%% Iterates through the windows to calculate the window average and scaling factor, and adjusts the reads accordingly
TA_window_indice=zeros(length(nonzeroindices),1);
TA_window_reads=zeros((windowfix),1);
counter=zeros((windowfix),1);

for x=1:windowfix;
    for p=1:length(nonzeroindices);
        if TAsites_nonzero(p) >= ((x*windowsize)-windowsize+1) & TAsites_nonzero(p) <= (x*windowsize);
            TA_window_indice(p) = x;
            counter(x)=counter(x)+1;
            TA_window_reads([x])=(reads_nonzero(p)+TA_window_reads([x])); % running total of all reads that fit in each window
        end
    end
end

mean_reads=TA_window_reads./counter; %for main windows, calculate the reads in each window
scale=chr_ave./mean_reads; %the scale factor is to be multiplied to the reads in each window

corrected_reads=zeros(length(TA_window_indice),1);
normalized_reads=zeros(length(Total_reads),1);

for x=1:length(TA_window_indice);
    if TA_window_indice(x) > 0;
        window_scalar=scale(TA_window_indice(x));
        corrected_reads(x)=round(reads_nonzero(x)*window_scalar); % scale each read to the correction factor calcuated for the window
        normalized_reads(nonzeroindices(x))=corrected_reads(x);
    end
end

%% Performs the same functions for the last leftover sequence

TA_window_indice=zeros(length(nonzeroindices),1); %resent TA_window_indice so you don't calculate the past window data again
TA_window_reads= zeros((windowfix+1),1);
counter=zeros((windowfix+1),1);

for p=1:length(nonzeroindices);
    if TAsites_nonzero(p) >= leftover_window_start & TAsites_nonzero(p) <= leftover_window_end;
        TA_window_indice(p) = windowfix+1;
        counter(windowfix+1)=counter(windowfix+1)+1;
        TA_window_reads([windowfix+1])=(reads_nonzero(p)+TA_window_reads([windowfix+1])); % running total of all reads that fit in each window
    end
end

mean_reads=TA_window_reads./counter;
scale=chr_ave./mean_reads;

for x=1:length(TA_window_indice);
    if TA_window_indice(x) > 0;
        window_scalar=scale(TA_window_indice(x));
        corrected_reads(x)=round(reads_nonzero(x)*window_scalar);
        normalized_reads(nonzeroindices(x))=corrected_reads(x);
    end
end
