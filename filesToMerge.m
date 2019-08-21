function [filePair]=filesToMerge (varargin)
%
%  filePair=filesToMerge (dirString,fileString)
% readFQfiles reads all fastq files in the specified directory
% input 1 is source directory, and
% input 2 is text that is part of the file name  like: '*.merged*.fastq'.
% Outputs includes cell array with the file names and the individual fastq
% data in a structure
dirString='D:\NGS\E1';
fileString='*.fastq';
if nargin>0
    dirString=varargin{1};
end
if nargin>1
    fileString=varargin{2};
end

cd (dirString);
if isempty(fileString)
    listing=[dir('*.merged*.fastq')];        % only read merged files
else
    listing=[dir(fileString)];
end
fileName={listing.name};

for i=size(fileName,2):-1:1                 % if we want to get rind of index files
    if ~isempty(strfind(fileName{i},'_I1_'))
        fileName{i}=[];
    end
end

pair=0;
filePair={};
for i=1:size(fileName,2)                 % if we want to get rind of index files
    if strfind(fileName{i},'_R1_')
        searchFor=strrep(fileName{i},'_R1_','_R2_');
        fileNo=find(strcmp(fileName,searchFor));
        if ~isempty(fileNo)
            pair=pair+1;
            filePair{pair,1}= fileName{i};
            filePair{pair,2}= fileName{fileNo};
        end
    end
end
end