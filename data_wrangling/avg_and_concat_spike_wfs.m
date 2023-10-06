function [averaged_wfs] = avg_and_concat_spike_wfs(wf_wire1,wf_wire2,...
    wf_wire3,wf_wire4)

%average spike waveform shapes and concatenate across wires; right now,
%this piece of code only works for tetrode recordings, and assumes a shape
%of n_time_samples x n_spikes


%take the average per wire
avg_wf1 = mean(wf_wire1,2);
avg_wf2 = mean(wf_wire2,2);
avg_wf3 = mean(wf_wire3,2);
avg_wf4 = mean(wf_wire4,2);


%concatenate each averaged waveform
averaged_wfs = [avg_wf1;avg_wf2;avg_wf3;avg_wf4];


end