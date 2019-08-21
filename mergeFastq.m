function [FQ_merged]=mergeFastq (varargin)

% [FQ]=mergeFastq (dirString,filePair,chQ)
% mergeFASTq merges R2 and R1 .fastq files to result in a merged sequence
% in R1 direction. Index read file is ignored.
% Inputs: i. dirstring: director location for input files; ii. fileString:
% if only specific subset of files have to be merged; and iii. minO: minimum
% overlap length for the merged sequence.
% Main output in: i. Result, containing merged
% sequence and additional information needed to output fastq file. Here, Quality
% Score is rewritten as 'J' if the same sequence is read in both R1 and R2. If
% there is a mismatch, sequence and quality score are based on better QS.
% Additional outputs: overCt if more than one read may correspond to a single R1
% read or underCt, if the paired read in R2 is absent

overlap=1;
minO=10;
overlapLen=25;  % minimum of 25 base overlap between sequences
sm = [1 -4 -4 -4; -4 1 -4 -4; -4 -4 1 -4; -4 -4 -4 1];  % scoring matrix
LenSeq1=1;
LenSeq2=1;
UMILen1=6;
UMILen2=6;

dirString='D:\NGS\E1';
filePair={{'RS-03300320_Lectin_S1_L001_R1_001.fastq'} {'RS-03300320_Lectin_S1_L001_R2_001.fastq'}};
% Qmin=0.05;      % 1:20 chance of having an error (to control Read quality)
% chQ=acceptableQ(Qmin);

if nargin>0
    dirString=varargin{1};
end
if nargin>1
    filePair=varargin{2};
end
if nargin>2
    chQ=varargin{3};
end

%[filePair]=filesToMerge (dirString,fileString);  % read file names, eliminating index reads, and pairing R1 and R2

for i=1:size(filePair,1)
    FQ1 = fastqread(filePair{i,1});
    FQ2 = fastqread(filePair{i,2});
    ProdFQ1={};
    ProdFQ2={};
    for k=1:length(FQ1)
        idx=findstr(FQ1(k).Header,':');
        ProdFQ1(k)={[FQ1(k).Header(idx(4)+1:idx(5)-1),FQ1(k).Header(idx(5)+1:idx(6)-1),FQ1(k).Header(idx(6)+1:idx(7)-2)]};
    end
    for k=1:length(FQ2)
        idx=strfind(FQ2(k).Header,':');
        ProdFQ2(k)={[FQ2(k).Header(idx(4)+1:idx(5)-1),FQ2(k).Header(idx(5)+1:idx(6)-1),FQ2(k).Header(idx(6)+1:idx(7)-2)]};
    end
    [~, idxMatch] = ismember(ProdFQ2, ProdFQ1);    % for every sequence in read 2 what is the index in read 1
    filePair{i,1},sum(idxMatch==0)                 % output if any read 2 did not have a match in read 1
    
    for k=1:length(FQ2)
        Seq1=FQ1(idxMatch(k)).Sequence;
        LenSeq1(length(Seq1)>LenSeq1)=length(Seq1);  % set length of Seq 1 read
        Qual1=FQ1(k).Quality;
%        [Seq1,Qual1]=cleanUP(Seq1,Qual1,chQ);
        Seq2=seqrcomplement(FQ2(k).Sequence);
        LenSeq2(length(Seq2)>LenSeq2)=length(Seq2);  % set length of Seq 2 read
        Qual2=reverse(FQ2(k).Quality);
%        [Seq2,Qual2]=cleanUP(Seq2,Qual2,chQ);
        merge=1;
        if ((overlap==1)&&(length(Seq1)>LenSeq1-5)&&(length(Seq2)>LenSeq2-5))                           % done only if overlap may occur and read lengths are complete
            loc1=[];
            loc2=[];
            if contains(Seq1,'N')
                loc1=strfind(Seq1,'N');
                if length(Seq1)-loc1(end)>minO
                    Seq1=Seq1(loc1(end)+1:end);
                else
                    merge=0;
                end
            end
            if contains(Seq2,'N')
                loc2=strfind(Seq2,'N');
                if loc2(1)>minO
                    Seq2=Seq2(1:loc2(1)-1);
                else
                    merge=0;
                end
          end
            if (merge~=0)
                struct = localalign(Seq1, Seq2, 'alpha','nt','scoringmatrix',sm,'gapopen', 100,'numaln', 3);  % Top 3 hits returned with scoring done with custom matrix; large gapopen penalty
                merge=0;
                if isempty(loc1)&&isempty(loc2)
                for j=1:length(struct.Alignment)
                    if ((length(Seq1)-struct.Stop(j,1))<6)&&(struct.Start(j,2)<6)&&(length(struct.Alignment{j})>overlapLen)  % within 6 bases and have overlap of 25
                        Seq_merged=[Seq1(1:struct.Stop(j,1)),Seq2(struct.Stop(j,2)+1:end)];
                        FQ1(k).Sequence=Seq_merged(UMILen1+1:end-UMILen2);
                        FQ1(k).Quality=[Qual1(UMILen1+1:struct.Stop(j,1)),Qual2(struct.Stop(j,2)+1:end-UMILen2)];
                        merge=struct.Stop(j,2);
                        FQ1(k).merge=merge;
                        FQ1(k).UMI=[Seq_merged(1:UMILen1),Seq_merged(end-UMILen2+1:end)];
                        FQ1(k).Seq1=Seq1;
                        FQ1(k).Seq2=Seq2;
                        FQ1(k).struct=struct.Alignment{j};
                    end
                end
                else 
                    FQ1(k).merge=0;  % place holder for now
                    FQ1(k).UMI='';
                  end
            end
        end
        
        if (overlap==0 || (merge==0)||(overlap==1 && (merge==1)))  % simply find reverse complement and merge
            Seq_merged=[Seq1,Seq2];
            FQ1(k).Sequence=Seq_merged(UMILen1+1:end-UMILen2);
            FQ1(k).Quality=[Qual1,Qual2];
            FQ1(k).merge=0;
            FQ1(k).UMI=[Seq_merged(1:UMILen1),Seq_merged(end-UMILen2+1:end)];
        end
    end
    FQ_merged{i}=FQ1;
end
end

% function [Seq,Qual]=mergeSeq(Seq1,Seq2,Qual1,Qual2)
% Seq=[Seq1,Seq2];
% Qual=[Qual1,Qual2];
% end

function [Seq,Qual]=cleanUP(Seq,Qual,chQ)
ismL=ismember(Qual,chQ);
if any(ismL==0)
    Seq(ismL==0)='N';
end
end

function chQ=acceptableQ(P)
ch1=-10*log10(P)+33;
count=0;
for i=floor(ch1):73
    count=count+1;
    chQ(count)=char(i);
end
end