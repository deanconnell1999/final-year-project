% Author: Karl Crofton
% Description: Removes the "Trigger" keyword from events in EEGLab's 
% data structure and correcting VPA trigger values.

function [data] = removeVPATrigger(data)
    cor = 1;
    for i = 1:length(data)
        tempType = str2num(cell2mat(regexp(data(i).type,"\d",'match')));
        if(tempType == 9 || tempType == 8)
            cor = 0;
            
            break;
        end
    end
    
    for i = 1:length(data)
        tempType = str2num(cell2mat(regexp(data(i).type,"\d",'match')));
        if(cor == 1)
            data(i).type = str2num(cell2mat(regexp(data(i).type,"\d",'match')));
        end
        
        if(cor == 0)
            if(tempType == 9 || tempType == 8 || tempType == 7 || tempType == 6)
                data(i).type = tempType - 2;
            else
                data(i).type = tempType;
            end
        end
    end
end
