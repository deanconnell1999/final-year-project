% Author: Sharon John
% Description: Creates time vector of timeOfEEG, the time (s) 
% corresponding to each index in EEG.data
function [timeOfEEG] = createTimeVector(t,maxTime)
    Index_length = 1:length(t);
    timeOfEEG = ((Index_length(:))/length(t))*maxTime;
    timeOfEEG = timeOfEEG';
end