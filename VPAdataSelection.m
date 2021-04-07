% Author: Karl Crofton
% Main file to process for VPA task, then export ERP and ERP Power Data.
% Instructions to run:
% - Run EEGLab and import desired files. For Each file, run this script.
% - Things you may wish to alter: the channels in the "chans" vector, The
%                                 folder for which you want to save the
%                                 file.
% NOTE: When you want to save export to a certain folder. Change line 76's
%       folder location.

clear patient_erp; %clearing patient_erp from buffer so we can start fresh
patient_erp.fcc(1,:) = 0;
patient_erp.fic(1,:) = 0;
patient_erp.tcc(1,:) = 0;
patient_erp.tic(1,:) = 0; %initialising the erps for each trigger value tcc, fcc etc.

% Remove string "Trigger-" from events and keeping event trigger number:
if(~(isa(EEG.event(1).type,'double')))
    EEG.event = removeVPATrigger(EEG.event); 
end

% Find patient number with regex (saving file name as this variable)
patient_erp.id = regexp(regexp(EEG.comments, '\d\d |\\\d\d\\', 'match'),'\d\d','match');

% Conversion from Samples to Seconds - we divide the latencies by sample rate
% to get time in seconds
clear e_secs %e_secs needs to be cleared (else old variable in workspace will interfere with new)
e_secs(1,:) = (cell2mat({EEG.event(:).latency})-1)/EEG.srate; % !! the first row of matrix e_secs is... EEG.event is reshaped into col vec

% produceTimeVect := Creates timeOfEEG, the time (s) corresponding to each
%                    index in EEG.data
timeOfEEG = createTimeVector(EEG.times,EEG.xmax); %.times: a vector the times data is recorded, .xmax: final time (in different decimal point)

% getIndexOfEvent:= Produces indx_st, the index values in EEG.data, 
%                     corresponding to each event in EEG.event
indx_st = getEventIndex(e_secs,timeOfEEG);

% Channel list (CHANGE IF USING FOR OTHER TESTS AND EXPERIMENTS)
chans = [14,18]; % PO3, PO4

%-------------- FOR TASK DATA ONLY -------------------
% VPA:
% chans(1) = PO3 data, chans(2) = PO4 data, these need to be changed to
% get ERPs for specific electrode.
eve = getVPAERP(EEG.data,chans(2),EEG.event,timeOfEEG,e_secs);

% % Gather and plot the ERP of tc, ti, fc & fi (for correct and wrong answers) : 
% 
% % From Question to Response Time
clear fct indx_st;
[vpa_erp] = getVPAERP(EEG.data,chans(2),EEG.event,timeOfEEG,e_secs);

fields = fieldnames(patient_erp);
map = ["fcc" "fic" "tcc" "tic"];

figure(1000);
for i = 1:(length(fieldnames(patient_erp))-1)
    subplot(2,2,i);
    t = (1:length(vpa_erp.(fields{i})))*1/1024;
    plot(t,vpa_erp.(fields{i}));
    title(strcat(patient_erp.id{1}," - ", upper(fields{i})));
    xlabel('Time (s)');
    ylabel('Voltage (mV)');
    
    % PATIENT ERP ASSIGNMENT
    patient_erp.(fields{i}) = vpa_erp.(fields{i});
end

for i = 1:4
    patient_erp.(strcat(map(i),"pow")) = periodogram(patient_erp.(map(i)),[],[],1024);
    patient_erp.(strcat(map(i),"pow")) = log10(patient_erp.(strcat(map(i),"pow")));
end

% ----------- Exporting the Data -----------------
% id = cell2mat(patient_erp.id{1});
% Change to desired folder
% id = strcat('your location',id,'.mat');
% save(id,'patient_erp');

figure(300);
periodogram(vpa_erp.fcc,rectwin(length(vpa_erp.fcc)),length(vpa_erp.fcc),EEG.srate);

% Clearing Temporary Variables:
clear curr fields c i biggest avg_len k sample_size f t; 