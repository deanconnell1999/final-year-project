% Author: Sharon John
% Description: Index values in EEG.data, corresponding to each 
% event in EEG.event
function [indx_st] = getEventIndex(eventSeconds,timeOfEEG)
    clear indx_st;
    
    indx_st = zeros(1,length(eventSeconds));
    for i = 1:length(eventSeconds)
        indx_st(i) = find((round(eventSeconds(1,i),1) == round(timeOfEEG,1)),1);
    end
end