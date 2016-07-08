%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Read files like condition_key.txt or task_contrasts.txt and
%               output as cell-matrix of character arrays
%
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-03-23
% REVISED:      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cellMat = readTaskFile(filePath)
%%

fileID = fopen(filePath);

tline = fgets(fileID);
i = 1;
while ischar(tline)
%     disp(tline)
    spaceInd = strfind(tline,' ');
    cellMat{i,1} = tline(1:spaceInd(1)-1);
    cellMat{i,2} = tline(spaceInd(1)+1:spaceInd(2)-1);
    cellMat{i,3} = tline(spaceInd(2)+1:end);
    tline = fgets(fileID);
    i = i + 1;
end


