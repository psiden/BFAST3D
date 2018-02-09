%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Read files like condition_key.txt, task_contrasts.txt or 
%               events.tsv and output as cell-matrix of character arrays
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-03-23
% REVISED:      2017-11-24
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cellMat = readTaskFile(filePath)
%%

fileID = fopen(filePath);

tline = fgets(fileID);
i = 1;
while ischar(tline)
%     disp(tline)
    spaceInd = sort([strfind(tline,' '),strfind(tline,char(9))]);
    cellMat{i,1} = tline(1:spaceInd(1)-1);
    cellMat{i,2} = tline(spaceInd(1)+1:spaceInd(2)-1);    
    if length(spaceInd) == 2
      cellMat{i,3} = tline(spaceInd(2)+1:end);
    else
      cellMat{i,3} = tline(spaceInd(2)+1:spaceInd(3)-1);
    end
    tline = fgets(fileID);
    i = i + 1;
end


