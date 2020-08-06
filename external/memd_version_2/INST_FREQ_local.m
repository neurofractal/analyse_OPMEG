function [instAmp,instFreq] = INST_FREQ_local(data)

% This function computes Hilbert-Huang spectrum using Hilbert transform,
% the instantaneous frequency and the instantaneous amplitude for each IMF. 
% 
% INPUT
% data : [N X M] IMF matrix 
%        (N : the number of IMF, M : data length)
%
% OUTPUT
% instAmp : [N X M] instantaneous amplitude matrix
%        (N : the number of IMF, M : data length)
% instFreq : [N X M] instantaneous frequency matrix
%        (N : the number of IMF, M : data length)


fs=1;
ts=1/fs;

dimension=size(data);

% initialise
instAmp = zeros(dimension(1),dimension(2));
instFreq = zeros(dimension(1),dimension(2));

for k=1:dimension(1)

    % Calculate Hilbert Transform
    % ```````````````````````````
    h=hilbert(data(k,:));


    % Instantatious Amplitude
    % ```````````````````````
    instAmp_temp = abs(h);
    instAmp(k,:)=instAmp_temp(:);
    
    % Instantanious Frequency
    % ```````````````````````
    phi = unwrap(angle(h));
    PHI(k,:)=(angle(h));
    instFreq_temp = (phi(3:end) - phi(1:end-2))/(2*ts);
    instFreq_temp = [instFreq_temp(1) instFreq_temp instFreq_temp(end)];
    instFreq(k,:)=instFreq_temp(:)/(2*pi);
    
end

instAmp(1)=instAmp(2);
instAmp(end)=instAmp(end-1);

instFreq=instFreq.*(double(instFreq>0));
instFreq = instFreq/fs;



