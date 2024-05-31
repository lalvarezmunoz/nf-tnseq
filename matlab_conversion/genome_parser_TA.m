%% Parse out the GTF file into appropriate fields
% FIrst, import GTF as a table, renaming the fields
% "Chrom" (col 1),"Type" (col3),"Start" (col4),"End"(5),"Name"(9)
% You will run this as many times as you have chromosomes, substituting the
% chromosome designation (and replacing the variable names) as you go

%% This section pulls the positions of all the entries related to genes (CDS) on each chromosome
function [Chr1_names, Chr1_starts, Chr1_ends, TAsites, TAid] = genome_parser_TA (GTFtable, chr_search_term, chrom_end, fasta,motif)
GTF_file=GTFtable;
Chromo=GTF_file.Chrom; %pull out the chromosome information
Chr_pos=strcmp(Chromo,chr_search_term); %creates a vector of logicals where very cell is yes (1) or no(0) on being on chrI
Type=GTF_file.Type; %pulls out the CDS,exon,attribute type data
CDS=strcmp('CDS',Type);
 
y=[]; %empty vector
for x=1:length(Chr_pos); %iterates through every row
    if Chr_pos(x) == 0
        z = 0;
        y=[y;z]; %if not on chromosome you want, calls it false and puts it in the Y vector
    else
        if strcmp('CDS',Type(x)) == 1 %if on chromosome you want, sees if type is CDS
            z = 1; % if CDS, then call it true
        else
            z = 0; % if not CDS, then call it false
        end
    y=[y;z]; %deposit the data into the y vector
    end
end
Chr_CDSpos=logical(y); %transform y vector in true/false logical

%% This section will take out each gene on the chromosome in question and their associated values that we want

Start_chr1=GTF_file.Start(Chr_CDSpos,1); %pulls out start sites for all cells with CDS (in that chromosome) in that row
End_chr1=GTF_file.End(Chr_CDSpos,1); %pulls out end sites for all cells with CDS (in that chromosome) in that row
Name_chr1=GTF_file.Name(Chr_CDSpos,1); %pulls out names/IDs for all cells with CDS (in that chromosome) in that row

gene_names={};
y=strfind(Name_chr1,'gene_id "'); %reports the position of where gene_id starts=usually 1
z=strfind(Name_chr1,'"');
for x=1:length(Name_chr1)
    searchlen=y{x}(1)+length('gene_id "');
    a=z{x}(find(z{x}>searchlen,1)); % finds which value in the vector corresponds to the first " downstream of 'gene_id "'.  Gives the term (1,2,3,4) NOT position
    b=Name_chr1{x}(searchlen:a-1); % finds the gene name between the quotation marks
    gene_names=[gene_names;b]; % deposits the names into a new vector; this is the list of gene names for the organism
end

[start_sort,sort_index]=sort(Start_chr1); %sorts by start site (start_sort); in sort_index, the first value of start_sort used to be at that position in Start_chr1

%% This section creates the IG (intergenic) elements between the annotated genes on the chromosome
IGstart=[];
IGend=[];
IGname={};

IGname=[IGname;'IG_1']; %designates the intergenic region before the first gene
IGstart=[IGstart;1];
IGend=[IGend;start_sort(1)-1];
prev_end=End_chr1(sort_index(1));

for x=2:length(start_sort);
    if start_sort(x) > prev_end % means we have at least 1 nt space between genes
        IGname=[IGname;strcat('IG_',gene_names(sort_index(x)))];
        IGstart=[IGstart;prev_end+1];
        IGend=[IGend;start_sort(x)-1];
        prev_end=End_chr1(sort_index(x));
    end
end
prev_end=End_chr1(sort_index(end));
if prev_end+1 < chrom_end
    IGstart=[IGstart;prev_end+1];
    IGend=[IGend;chrom_end];
    IGname=[IGname;'IG_chrm_end'];
end

% Now to add the IGs to the Genes and re-sort

c=[Start_chr1;IGstart];
[Chr1_starts,final_sort_index]=sort(c);
d=[gene_names;IGname];
e=[End_chr1;IGend];
Chr1_ends=[];
Chr1_names={};

for x=1:length(Chr1_starts)
    ends=e(final_sort_index(x));
    names=d(final_sort_index(x));
    Chr1_ends=[Chr1_ends;ends];
    Chr1_names=[Chr1_names;names];
end
% use some clear functions to clean up workspace?
%% This section next finds all the TA sites from a FASTA sequence, reports them and links them with which locus they are in.

g=fastaread(fasta); %reads in fasta file.  Have to type this out in command window (copy URL), no options selected
% fasta data is located in 2 values:  g.Header = chromosome info;
% g.Sequence is the full sequence in 1 single string

TAsites=strfind(g.Sequence,motif); % finds TA sites in the FASTA file
TAsites=TAsites'; %Transposes the TA sites into a column

% Maps each TA to the locus in which it's located

TA_index=zeros(length(TAsites),1); %creates an array of 0's the length to the TA site list

for x=1:length(TAsites); %For every TA site, find the start and end combination it's between; runs ~7 min for 140000 TA and 5600 loci
    for p=1:length(Chr1_starts); 
        if any(TAsites(x)>=Chr1_starts(p) & TAsites(x)<=Chr1_ends(p)) == 1;
            TA_index([x])=[p]; % if the TA is between the starts and ends, replace the appropriate cell with the position of the locus
            break %finds first value that fits criteria, reports it and stops the loop to save time
        end
    end
end

TAid={}; % creates a file in which every TAsite is now linked to its locus name (4-5 min on 4GB machine) 
for x=1:length(TAsites);
    TAid=[TAid;Chr1_names(TA_index(x))];
end
