% Author: Karl Crofton
% Description: Takes data, channel number, event numbers, time of data 
% (in seconds) and event times (in seconds). Returns a struct containing
% ERPs for all four correct reponse types (fcc,fic,tcc,tic) for a patient
% on a particular channel in VPA task.

function [VPAerps] = getVPAERP(data,chans,events,timeOfEEG,eventSeconds)
    for i = 1:length(events)
        if(events(i).type > 1)
            c = i;
            break;
        end
    end
    
    curr = 0;
    clear eve; % DELETE AFTER USE (Or 'eve' won't be accurate)
    eve.fcc(1:2,1) = 0;
    eve.tcc(1:2,1) = 0;
    eve.fic(1:2,1) = 0;
    eve.tic(1:2,1) = 0;
    map = ["fcc" "fic" "tcc" "tic"];
    for i = c:2:length(events)
        curr = curr+1;
        eve.ques(curr) = events(1,i).type;
        eve.resp(curr) = events(1,i+1).type;
        if(eve.resp(curr) == 10)
            currType = map(eve.ques(curr)-3);
            eve.(currType)(1,length(eve.(currType)(1,:))) = i;
            eve.(currType)(2,length(eve.(currType)(2,:))) = i+1;
            eve.(currType)(:,length(eve.(currType)(1,:))+1) = 0;
        end
    end
    eve.fcc(:,end) = [];
    eve.tcc(:,end) = [];
    eve.fic(:,end) = [];
    eve.tic(:,end) = [];
    
    eve = rmfield(eve,'resp');
    eve = rmfield(eve,'ques');
    
    fields = fieldnames(eve)
    
    for k = 1:numel(fields)
        biggest = 0;
        indx_st = zeros(1,length(eve.(fields{k})));
        diffs = zeros(1,length(eve.(fields{k})));
        
        for i =1:length(eve.(fields{k}))
            indx_st(1,i) = find(round(eventSeconds(eve.(fields{k})(1,i)),1) == round(timeOfEEG,1),1);
            indx_st(2,i) = find(round(eventSeconds(eve.(fields{k})(2,i)),1) == round(timeOfEEG,1),1);
            diffs(i) = length(indx_st(1,i):indx_st(2,i));
            if(biggest<diffs(i))
                biggest = diffs(i);
            end
        end
        
        % Find the index value of event - 200ms and event + 2500ms
        indicesBefore = find(round(eventSeconds(eve.(fields{k})(1,1)),1) == round(timeOfEEG,1),1) - find(round(eventSeconds(eve.(fields{k})(1,1))-0.2,1) == round(timeOfEEG,1),1);
        indicesAfter = find(round(eventSeconds(eve.(fields{k})(1,1)) + 2.5,1) == round(timeOfEEG,1),1) - find(round(eventSeconds(eve.(fields{k})(1,1)),1) == round(timeOfEEG,1),1);
        
        clear frame; % Frame needs to be cleared every loop
        
        for i = 1:length(eve.(fields{k}))
            % snippet: question - 200ms : question + 2500ms
            snippet = (indx_st(1,i) - indicesBefore):(indx_st(1,i) + indicesAfter);
            frame(i,:) = data(chans,snippet);
        end
        
        for i=1:length(frame)
            vpa_erp(i) = sum(frame(:,i))/nnz(frame(:,i));
        end
        
        VPAerps.(map(k)) = vpa_erp;
    end
end