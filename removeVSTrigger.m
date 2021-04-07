% Author: Sharon John
% Description: Removes the "Trigger" keyword from events in EEGLab's 
% data structure for VS task.

function [data] = removeVSTrigger(data)
    for i = 1:length(data)
        tempType = str2double(cell2mat(regexp(data(i).type,"\d|\d\d|\d\d\d",'match')));
        data(i).type = tempType;
    end
end