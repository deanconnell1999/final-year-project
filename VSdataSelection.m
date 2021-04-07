% Author: Sharon John
% Main file to process for VS task, then export ERP and ERP Power Data.
% Instructions to run:
% - Run EEGLab and import desired files. For Each file, run this script.
% - Things you may wish to alter: the channels in the "channels" vector, 
%                                 The folder for which you want to save 
%                                 the file.
% NOTE: When you want to save export to a certain folder. Change line 84's
%       folder location.

%Trigger removal
if(~(isa(EEG.event(1).type,'double')))
    EEG.event = removeVSTrigger(EEG.event);
end

% Converting latencies to seconds
% time vector of events:
% Filtering

clear patient_erp;  %creating class for patient erp data 
patient_erp.flc(1,:) = 0;
patient_erp.lc(1,:) = 0;
patient_erp.rc(1,:) = 0;
patient_erp.frc(1,:) = 0;
patient_erp.ac(1,:) = 0;

patient_erp.id = regexp(regexp(EEG.comments, '\d\d |\\\d\d\\', 'match'),'\d\d','match');
% Conversion from Samples to Seconds
clear event_secs
%dividing the latency by sample rate to get to seconds 
event_secs(1,:) = (cell2mat({EEG.event(:).latency})-1)/EEG.srate; 

% createTimeVect := Creates timeOfEEGData, the time (s) corresponding to each index in EEG.data
timeOfEEG = createTimeVector(EEG.times,EEG.xmax);

% Creates indx_st := Index values in EEG.data, corresponding to each event in EEG.event
indx_st = getEventIndex(event_secs,timeOfEEG);
%--------------------------------------------------------------------------
%Above 2 statements are together so you are producing a time vector and comparing that with second line to get an index  
%--------------------------------------------------------------------------

% Channel Data Selection:
channels = [7,23]; % Getting channels P7 and P8.

% Visual Search Task ERP:
% channels(1) = P7 data, channels(2) = P8 data, these need to be changed to
% get ERPs for specific electrode.
eve = getVSERP(EEG.data,channels(2),EEG.event,timeOfEEG,event_secs); 
%need to change the channels to get the output for the channels at P7, P8.

% % Gather and plot the ERP of flc, lc, rc, frc, & ac (for correct and wrong answers) : 
% % From Question to Response Time
clear fct indx_value;
[vs_erp] = getVSERP(EEG.data,channels(2),EEG.event,timeOfEEG,event_secs);

%Getting the ERP data for each patient at each channel and plotting it
header = fieldnames(patient_erp);
map = ["flc" "lc" "rc" "frc" "ac"];

figure(1000);
for i = 1:(length(fieldnames(patient_erp))-1)
    subplot(2,3,i);
    t = (1:length(vs_erp.(header{i})))*1/1024;
    plot(t,vs_erp.(header{i}));
    title(strcat(patient_erp.id{1}," - ", upper(header{i})));
    xlabel('Time (s)');
    ylabel('Voltage (mV)');
    
    % Patient ERP Assignment
    patient_erp.(header{i}) = vs_erp.(header{i});
end

%Mapping the periodogram and log10 to the new vector that you create on 
%responsepow (e.g. flcmap) - just for exporting so you can put it into 
%python

for i =1:5
    patient_erp.(strcat(map(i), "map")) = periodogram(patient_erp.(map(i)), [], [], 1024);
    patient_erp.(strcat(map(i), "map")) = log10(patient_erp.(strcat(map(i), "map")));
end 

% ----------- Exporting the Data -----------------
id = cell2mat(patient_erp.id{1});
% Change to desired folder
id = strcat('your location',num2str(id), '.mat');
save(id,'patient_erp');

figure(300);
periodogram(vs_erp.flc,rectwin(length(vs_erp.flc)),length(vs_erp.flc),EEG.srate);

% Clearing Temporary Variables:
clear curr fields c i biggest avg_len k sample_size f t; 
