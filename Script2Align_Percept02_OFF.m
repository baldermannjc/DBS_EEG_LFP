%% loading data
clear; close all; clc

% add folder with lfp_resample.m to matlab paths 

% in this example Percept recordings were started earlier and ended earlier
% than the EEG recordings. This script is (as of now) specifically build
% for that case. In other cases, a part of the LFP signal needs to be cut
% (probably easier) or the EEG channels need zero-padding.

eegfile = 'Percept_02_OFF';
eegpath = '/Users/thomas/Desktop/Fingertapping/01_Data/01_EEG/01_Raw/';
savepath = '/Users/thomas/Desktop/Fingertapping/01_Data/01_EEG/02_LFP_included/';

% load LFP data
lfp = load('/Users/thomas/Desktop/Fingertapping/01_Data/02_LFP/Percept_02_LFP.mat');

% import EEG data and high-pass filter
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab('nogui');
EEG = pop_loadbv(eegpath, [eegfile '.vhdr'],[],[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'Percept_02_OFF');
EEG = eeg_checkset(EEG);

EEG = pop_eegfiltnew(EEG, 1, []);
EEG = eeg_checkset(EEG);

%% LFP data at .BrainSenseTimeDomain

% index
xindx = 1;

% sampling rate
Fs = lfp.json.BrainSenseTimeDomain(1).SampleRateInHz;

% in sets 33 and 34 there was finger tapping in the left and right
% hemispheres
for i = 33:34
    
    % load time domain LFP data
    lfp_raw = lfp.json.BrainSenseTimeDomain(i).TimeDomainData;
    
    % time index (only needed for plotting with different sampleing rates)
    lfp_time = (1:length(lfp_raw))*1/Fs;
    
    % channel name & name parts
    chaname = lfp.json.BrainSenseTimeDomain(i).Channel;
    [chan.one,chan.two] = strtok(chaname,'_');
    [chan.three,chan.four] = strtok(chan.two,'_');
    
    % gain (??) (is apparently not used after that)
    gain = lfp.json.BrainSenseTimeDomain(i).Gain;
    
    eval(['lfp' lower(chan.four) ' = lfp_raw;']);
        
%     figure
%     plot(lfp_time, lfp_raw,'LineWidth',2)
%     title([num2str(i) ' ' chaname ' --> total length =' num2str(lfp_time(end))]);

end

% resample LFP to EEG sampling rate (1000 Hz)
LFP_right = lfp_resample(lfp_right,250,1000);
LFP_left = lfp_resample(lfp_left,250,1000);

% plot first filter artifacts
figure
subplot(1,2,1)
plot(EEG.data(17,1:2*10^4)); hold on
plot(LFP_left(1:2*10^4))
set(gca,'ylim',[-700 700])
title('Left')
subplot(1,2,2)
plot(EEG.data(17,1:2*10^4)); hold on
plot(LFP_right(1:2*10^4))
set(gca,'ylim',[-700 700])
title('Right')

% align LFP and EEG artifacts

lfp_shift = zeros(1,.266*10^4);

figure
plot(EEG.data(17,.4*10^4:1.5*10^4)); hold on
plot([lfp_shift LFP_right(.4*10^4:1.5*10^4)])
set(gca,'ylim',[-700 700])

figure
plot(EEG.data(17,.4*10^4:1.5*10^4)); hold on
plot([lfp_shift LFP_left(.4*10^4:1.5*10^4)])
set(gca,'ylim',[-700 700])

% plot whole recording to check were to start plotting next
figure
plot(EEG.data(17,:)); hold on
plot(LFP_right)

% check artifact at the end of the recording
figure
subplot(1,2,1)
plot(EEG.data(17,3.6*10^5:end)); hold on
plot([lfp_shift LFP_right(3.6*10^5:end)])
set(gca,'ylim',[-700 700])
title('Left')
subplot(1,2,2)
plot(EEG.data(17,3.6*10^5:end)); hold on
plot([lfp_shift LFP_left(3.6*10^5:end)])
set(gca,'ylim',[-700 700])
title('Right')

% add missing data points to LFP (zero-padding): this assumes Percept
% recording was started earlier than EEG recording
LFP_right = [lfp_shift LFP_right]; %add zeros to shift LFP
LFP_left = [lfp_shift LFP_left];

% this assumes LFP recording was stopped earlier than the EEG recording. 
if size(LFP_right,2) < size(EEG.data,2) % check if LFP shorter than EEG
    lfp_pad = zeros(1,size(EEG.data,2)- size(LFP_right,2));
    
    LFP_right = [LFP_right lfp_pad]; %zero-padding to have same array size
    LFP_left = [LFP_left lfp_pad];
end

% append two EEG channels and fill with LFP data
EEG.data(end+1,:) = LFP_right;
EEG.nbchan = size(EEG.data,1);
EEG.chanlocs(end+1).labels = 'LFP_right';

EEG.data(end+1,:) = LFP_left;
EEG.nbchan = size(EEG.data,1);
EEG.chanlocs(end+1).labels = 'LFP_left';

% save EEGLAB dataset
pop_saveset( EEG, 'filename', [savepath  eegfile '.set']);

% %topo of stimulation artifact
% figure
% topoplot(EEG.data(:,1.2*10^4),EEG.chanlocs,'conv','off','electrodes','off')
% set(gca,'clim',[-500 150])
% colorbar



