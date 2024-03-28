function [LFP] = lfp_resample(data, old_freq, new_freq)

% % old sampling rate
% old_freq = 250;
% 
% % new sampling rate
% new_freq = 1000;
 
% Default cutoff frequency (pi rad / smp)
fc = 0.9; 

% Default transition band width (pi rad / smp)
df = 0.2;

% finding the best ratio
[p,q] = rat(new_freq/old_freq, 1e-12);

%boundaries
bounds = [1 size(data,1) + 1]; % Add initial and final boundary event

tmpeeglab = [];

for index1 = 1:size(data,2)
    
    sigtmp = data;
    
    if index1 == 1
        tmpres = [];
        indices = 1;
        for ind = 1:length(bounds)-1
            tmpres  = [ tmpres; lfp_myresample( double( sigtmp(bounds(ind):bounds(ind+1)-1,:)), p, q, fc, df ) ];
            indices = [ indices size(tmpres,1)+1 ];    
        end
        
    else
        for ind = 1:length(bounds)-1
            tmpres(indices(ind):indices(ind+1)-1,:) = lfp_myresample( double( sigtmp(bounds(ind):bounds(ind+1)-1,:) ), p, q, fc, df);
        end
    end 
    
    tmpeeglab(index1,:, :) = tmpres;
    
end

LFP = tmpeeglab;

% LFP.data = tmpeeglab;
% LFP.pnts = size(LFP.data,2);
% LFP.srate = old_freq*p/q;
% LFP.xmin = round(min((1:length(data)*1)/old_freq));
% LFP.xmax = LFP.xmin + (LFP.pnts-1)/LFP.srate;
% LFP.times = linspace(LFP.xmin*1000,LFP.xmax*1000,LFP.pnts);

end

% resample with signal processing toolbox
% -----------------------------------
function tmpeeglab = lfp_myresample(data, p, q, fc, df)

if length(data) < 2
    tmpeeglab = data;
    return;
end

% Conservative custom anti-aliasing FIR filter, see bug 1757
nyq = 1 / max([p q]);
fc = fc * nyq; % Anti-aliasing filter cutoff frequency
df = df * nyq; % Anti-aliasing filter transition band width
m = pop_firwsord('kaiser', 2, df, 0.002); % Anti-aliasing filter kernel
b = firws(m, fc, windows('kaiser', m + 1, 5)); % Anti-aliasing filter kernel
b = p * b; % Normalize filter kernel to inserted zeros

% Padding, see bug 1017
nPad = ceil((m / 2) / q) * q; % Datapoints to pad, round to integer multiple of q for unpadding
startPad = repmat(data(1, :), [nPad 1]);
endPad = repmat(data(end, :), [nPad 1]);

% Resampling
tmpeeglab = resample([startPad; data; endPad], p, q, b(:));

% Remove padding
nPad = nPad * p / q; % # datapoints to unpad
tmpeeglab = tmpeeglab(nPad + 1:end - nPad, :); % Remove padded data

end