function amplicon
Qmin=0.05;      % 1:20 chance of having an error (to control Read quality)
chQ=acceptableQ(Qmin);
Indel=[];
NoIndel=[];
File={};
gCSV='C:\Users\neel\Desktop\India\molther\Multiplex PCR NGS analysis input table_071919.csv'; % read sgRNA for each gene from .csv
[FPN,FSeq,RPN,RSeq,GSeq,Strand,PCRprod,checkSeq,checkLoc]=getData(gCSV);
dirString='C:\Users\neel\Desktop\India\molther\RPCI071819';   % read fastq data file
fileString='*Dox*R*.fastq';
filePair=filesToMerge (dirString,fileString);
%   FQ=mergeFastq (dirString,{filePair{i,:}},chQ);
FQ=mergeFastq (dirString,filePair,chQ);
for i=1:length(FQ)
    [noindel_ct,indel_ct,FQ{i}]=indel(FPN,FSeq,RSeq,GSeq,PCRprod,checkSeq,checkLoc,Strand,FQ{i});
    File{i}=filePair(i,1);
    Indel=[Indel;indel_ct'];
    NoIndel=[NoIndel;noindel_ct'];
    [filePair(i,1),'total',length(FQ{i}),'indel',sum(Indel(i,:)),'noIndel',sum(NoIndel(i,:))]
    FQ_temp=FQ{i};
    cell_Indel={FQ_temp.Indel};
    cell_IndelLen={FQ_temp.IndelLen};
    arrIndel=cell2mat(cell_Indel);
    arrIndelLen=cell2mat(cell_IndelLen);
    ['changes to GSeq',{sum(arrIndel>0)},'Length change F-Rprimer',sum(abs(arrIndelLen)>0)]
end
    Percent=100*Indel./(Indel+NoIndel);
end

function [Fwd_Primer_Name,FSeq,Rev_Primer_Name,RSeq,GSeq,Strand,amplicon,checkSeq,checkLoc]=getData(gCSV)
fileID = fopen(gCSV);
Data=textscan(fileID,'%s %s %s %s %s %s %s','headerLines', 1,'Delimiter',',');
fclose(fileID)
Fwd_Primer_Name=Data{1,1};
FSeq=lower(Data{1,2});
Rev_Primer_Name=Data{1,3};
RSeq=lower(Data{1,4});
RSeq=cellfun(@(x) seqrcomplement(x), RSeq,'un',0);
GSeq=lower(Data{1,5});
Strand=Data{1,6};
amplicon=lower(Data{1,7});
Pos=[];
checkLoc=[];
for i=1:length(Fwd_Primer_Name)
    Pos(i)=strfind(amplicon{i},GSeq(i));
    if (Strand{i}=='+')
        checkSeq{i}=amplicon{i}(Pos(i)-20:Pos(i)-1);
        checkLoc=[checkLoc;Pos(i)-20,Pos(i)-1];
    else
        checkSeq{i}=amplicon{i}(Pos(i)+20:Pos(i)+39);
        checkLoc=[checkLoc;Pos(i)+20,Pos(i)+39];
    end
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

function [noindel_ct,indel_ct,FQ]=indel(FPN,FSeq,RSeq,GSeq,PCRprod,checkSeq,checkLoc,Strand,FQ)

condition=3;   % based on sgRNA

count=zeros(length(FPN),1);
countLen=zeros(length(FPN),1);
countBoth=zeros(length(FPN),1);
UMIcount_noindel=zeros(length(FPN),1);
UMI_noindel={};
UMIcount_indel=zeros(length(FPN),1);
UMI_indel={};
UMIcount_noindel_Len=zeros(length(FPN),1);
UMI_noindel_Len={};
UMIcount_indel_Len=zeros(length(FPN),1);
UMI_indel_Len={};
UMIcount_noindel_both=zeros(length(FPN),1);
UMI_noindel_both={};
UMIcount_indel_both=zeros(length(FPN),1);
UMI_indel_both={};
indelStore={};
for i=1:length(FQ)
    if (FQ(i).merge~=0)
        Seq=lower(FQ(i).Sequence);
        y=cellfun(@(x) strfind(Seq,x), FSeq,'un',0);
        idx1=find(~cellfun(@isempty,y));
        y=cellfun(@(x) strfind(Seq,x), RSeq,'un',0);
        idx2=find(~cellfun(@isempty,y));
        if (~isempty(idx1)&&~isempty(idx2))
            if (length(idx1)~=1) || (length(idx2)~=1)
                FQ(i).Indel=2;
            elseif idx1~=idx2
                FQ(i).Indel=2;
            else
                theoLen=strfind(PCRprod{idx1},RSeq(idx1))-strfind(PCRprod{idx1},FSeq(idx1));
                exptLen=strfind(Seq,RSeq{idx1})-strfind(Seq,FSeq{idx1});
                if (length(theoLen)>1) || (length(exptLen)>1)
                    FQ(i).Indel=2;
                else
                    FQ(i).IndelLen=exptLen-theoLen;
                    if (FQ(i).IndelLen==0)
                        UMIcount_noindel_Len(idx1)=UMIcount_noindel_Len(idx1)+1;
                        UMI_noindel_Len{idx1,UMIcount_noindel_Len(idx1)}=FQ(i).UMI;
                    else
                        UMIcount_indel_Len(idx1)=UMIcount_indel_Len(idx1)+1;
                        UMI_indel_Len{idx1,UMIcount_indel_Len(idx1)}=FQ(i).UMI;
                    end
                    countLen(idx1)=countLen(idx1)+1;
                    if(length(Seq)>checkLoc(idx1,2)+length(GSeq{idx1}))
                        %                    if strcmp(checkSeq{idx1},Seq(checkLoc(idx1,1):checkLoc(idx1,2)))
                        count(idx1)=count(idx1)+1;
                        guideTarget=GSeq{idx1};
                        FQ(i).match=idx1;
                        if (Strand{idx1}=='+')
                            guideActual=Seq(checkLoc(idx1,2)+1:checkLoc(idx1,2)+length(guideTarget));
                        else
                            guideActual=Seq(checkLoc(idx1,1)-20:checkLoc(idx1,1)-21+length(guideTarget));
                        end
                        if strcmp(guideTarget,guideActual)  % no indel
                            UMIcount_noindel(idx1)=UMIcount_noindel(idx1)+1;
                            UMI_noindel{idx1,UMIcount_noindel(idx1)}=FQ(i).UMI;
                            FQ(i).Indel=0;
                        else
                            UMIcount_indel(idx1)=UMIcount_indel(idx1)+1;
                            UMI_indel{idx1,UMIcount_indel(idx1)}=FQ(i).UMI;
                            indelStore{idx1,UMIcount_indel(idx1)}=guideActual;
                            FQ(i).Indel=1;
                            % indel=length(Seq)-PCRprod(idx1);
                        end
                        %                    end
                    if ((FQ(i).IndelLen~=0)&&(FQ(i).Indel==1))
                            UMIcount_indel_both(idx1)=UMIcount_indel_both(idx1)+1;
                            UMI_indel_both{idx1,UMIcount_indel_both(idx1)}=FQ(i).UMI;
                    else
                            UMIcount_noindel_both(idx1)=UMIcount_noindel_both(idx1)+1;
                            UMI_noindel_both{idx1,UMIcount_noindel_both(idx1)}=FQ(i).UMI;
                    end
                    end
                end
            end
        end
    end
end
switch condition
    case 1
        [noindel_ct,indel_ct]=countIndel(UMI_indel,UMI_noindel,FPN);
    case 2
        [noindel_ct,indel_ct]=countIndel(UMI_indel_Len,UMI_noindel_Len,FPN);
    case 3
        [noindel_ct,indel_ct]=countIndel(UMI_indel_both,UMI_noindel_both,FPN);
end
end

function [noindel_ct,indel_ct]=countIndel(UMI_indel,UMI_noindel,FPN)
noindel_ct=zeros(length(FPN),1);
indel_ct=zeros(length(FPN),1);
for i=1:size(UMI_indel,1)
    y=UMI_indel(i,:);
    y=y(~cellfun('isempty',y));
    y=unique(y);
    indel_ct(i)=length(y);
end
for i=1:size(UMI_noindel,1)
    y=UMI_noindel(i,:);
    y=y(~cellfun('isempty',y));
    y=unique(y);
    noindel_ct(i)=length(y);
end
end
