% Author: Sharon John
% Description: Takes data, channel number, event numbers, time of data 
% (in seconds) and event times (in seconds). Returns a struct containing
% ERPs for all five correct reponse types (flc,lc,rc,frc,ac) for a patient
% on a particular channel in VS task.

function [VSerps, avg_len] = getVSERP(data,channels,events,timeOfEEG,eventSeconds)
    
    %Gathering question and response stimulus arrays from data.
    curr = 0;
    clear eve;
    eve.flc(1:2,1) = 0;
    eve.lc(1:2,1) = 0;
    eve.rc(1:2,1) = 0;
    eve.frc(1:2,1) = 0;
    eve.ac(1:2,1) = 0;
    map = ["flc" "lc" "rc" "frc" "ac"];
    for i = 1:2:length(events)
        curr = curr+1;
        eve.ques(curr) = events(1,i).type;
        eve.resp(curr) = events(1,i+1).type;
        if(eve.resp(curr) == 10)
            if (eve.ques(curr) == 128) 
                currType = map(5);
            else
                currType = map(eve.ques(curr));
            end
            eve.(currType)(1,length(eve.(currType)(1,:))) = i;
            eve.(currType)(2,length(eve.(currType)(2,:))) = i+1;
            eve.(currType)(:,length(eve.(currType)(1,:))+1) = 0;
        end
    end
    eve.flc(:,end) = [];
    eve.lc(:,end) = [];
    eve.rc(:,end) = [];
    eve.frc(:,end) = [];
    eve.ac(:,end) = [];
    eve = rmfield(eve,'resp');
    eve = rmfield(eve,'ques');
    
    header = fieldnames(eve);
    
    for k = 1:numel(header)
        biggest = 0;
        indx_value = zeros(1,length(eve.(header{k})));
        diffs = zeros(1,length(eve.(header{k})));
        
        for i =1:length(eve.(header{k}))
            indx_value(1,i) = find(round(eventSeconds(eve.(header{k})(1,i)),1) == round(timeOfEEG,1),1);
            indx_value(2,i) = find(round(eventSeconds(eve.(header{k})(2,i)),1) == round(timeOfEEG,1),1);
            diffs(i) = length(indx_value(1,i):indx_value(2,i));
            if(biggest<diffs(i))
                biggest = diffs(i);
            end
        end
        
        % Takes the time of the event and minus' 200ms before and gets the index
        indicesIndexBefore = find(round(eventSeconds(eve.(header{k})(1,1)),1) == round(timeOfEEG,1),1) - find(round(eventSeconds(eve.(header{k})(1,1))-0.2,1) == round(timeOfEEG,1),1);
        % Takes the time of the event and minus' 2500ms after and gets the index
        indicesIndexAfter = find(round(eventSeconds(eve.(header{k})(1,1)) + 2.5,1) == round(timeOfEEG,1),1) - find(round(eventSeconds(eve.(header{k})(1,1)),1) == round(timeOfEEG,1),1);
        
        clear frame; % Frame needs to be cleared every loop
        
        for i = 1:length(eve.(header{k}))
            % snippet: question - 200ms : question + 2500ms
            snippet = (indx_value(1,i) - indicesIndexBefore):(indx_value(1,i) + indicesIndexAfter);
            frame(i,:) = data(channels,snippet);
        end
        
        for i=1:length(frame)
            vs_erp(i) = sum(frame(:,i))/nnz(frame(:,i));
        end
        
        VSerps.(map(k)) = vs_erp;
    end
end