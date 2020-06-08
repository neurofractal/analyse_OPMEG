function [filt] = ft_preproc_dftfilter_NA(pow, freq, dat, Fs, lowerFreqLim)

% FT_PREPROC_DFTFILTER reduces power line noise (50 or 60Hz) via two 
% alternative methods:
% A) DFT filter (Flreplace = 'zero') or
% B) Spectrum Interpolation (Flreplace = 'neighbour').
%
% A) The DFT filter applies a notch filter to the data to remove the 50Hz
% or 60Hz line noise components ('zeroing'). This is done by fitting a sine 
% and cosine at the specified frequency to the data and subsequently 
% subtracting the estimated components. The longer the data is, the sharper 
% the spectral notch will be that is removed from the data.
% Preferably the data should have a length that is a multiple of the
% oscillation period of the line noise (i.e. 20ms for 50Hz noise). If the
% data is of different lenght, then only the first N complete periods are
% used to estimate the line noise. The estimate is subtracted from the
% complete data.
%
% B) Alternatively line noise is reduced via spectrum interpolation
% (Leske & Dalal, 2019, NeuroImage 189,
%  doi: 10.1016/j.neuroimage.2019.01.026)
% The signal is:
% I)   transformed into the frequency domain via a discrete Fourier 
%       transform (DFT), 
% II)  the line noise component (e.g. 50Hz, Flwidth = 1 (±1Hz): 49-51Hz) is 
%       interpolated in the amplitude spectrum by replacing the amplitude 
%       of this frequency bin by the mean of the adjacent frequency bins 
%       ('neighbours', e.g. 49Hz and 51Hz). 
%       Neighwidth defines frequencies considered for the mean (e.g. 
%       Neighwidth = 2 (±2Hz) implies 47-49 Hz and 51-53 Hz). 
%       The original phase information of the noise frequency bin is
%       retained.
% III) the signal is transformed back into the time domain via inverse DFT
%       (iDFT).
% If Fline is a vector (e.g. [50 100 150]), harmonics are also considered. 
% Preferably the data should be continuous or consist of long data segments
% (several seconds) to avoid edge effects. If the sampling rate and the
% data length are such, that a full cycle of the line noise and the harmonics
% fit in the data and if the line noise is stationary (e.g. no variations
% in amplitude or frequency), then spectrum interpolation can also be 
% applied to short trials. But it should be used with caution and checked 
% for edge effects.
%
% Use as
%   [filt] = ft_preproc_dftfilter(dat, Fsample, Fline, varargin)
% where
%   dat             data matrix (Nchans X Ntime)
%   Fsample         sampling frequency in Hz
%   Fline           line noise frequency (and harmonics)
%
% Additional input arguments come as key-value pairs:
%
%   Flreplace       'zero' or 'neighbour', method used to reduce line noise, 'zero' implies DFT filter, 'neighbour' implies spectrum interpolation  
%   Flwidth         bandwidth of line noise frequencies, applies to spectrum interpolation, in Hz
%   Neighwidth      width of frequencies neighbouring line noise frequencies, applies to spectrum interpolation (Flreplace = 'neighbour'), in Hz 
%
% The line frequency should be specified as a single number for the DFT filter.
% If omitted, a European default of 50Hz will be assumed
%
% See also PREPROC

% Undocumented option:
%   Fline can be a vector, in which case the regression is done for all
%   frequencies in a single shot. Prerequisite is that the requested
%   frequencies all fit with an integer number of cycles in the data.
%
% Copyright (C) 2003, Pascal Fries
% Copyright (C) 2003-2015, Robert Oostenveld
% Copyright (C) 2016, Sabine Leske 
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$
%
% Edited by Nicholas Alexander to customised/uneven peak widths.

% determine the size of the data
[~, nsamples] = size(dat);

% Calculate PSD
% cfgPSD              = [];
% cfgPSD.channel      = vertcat(ft_channelselection_opm('TAN',dat));
% cfgPSD.trial_length = 3;
% cfgPSD.method       = 'tim';
% cfgPSD.foi          = [1 200];
% cfgPSD.plot         = 'yes';
% [pow, freq, ~]      = ft_opm_psd(cfgPSD,dat);


% Important to treat all sensors/epochs the same so median across both.
powMed              = median(pow(:,:,:),3);
powMed              = median(pow(:,:),2);

% Identify peaks and remove those below the lowerFreqLim
[peakVal,peakLoc]   = findpeaks(powMed,freq,'MinPeakProminence',10);
peaksToRemove       = (peakLoc < lowerFreqLim);
peakVal(peaksToRemove) = [];
peakLoc(peaksToRemove) = [];


% Find negative peaks for peak edges.
[negPeakVal,negPeakLoc] = findpeaks(-powMed,freq);
lowerPeakLoc            = zeros(size(peakLoc));
upperPeakLoc            = zeros(size(peakLoc));
lowerPeakVal            = zeros(size(peakLoc));
upperPeakVal            = zeros(size(peakLoc));
for peakIdx = 1:length(peakLoc)
    % Find the nearest negPeak to the left of peakLoc
    lowerPeakIdx            = find(negPeakLoc < peakLoc(peakIdx),1,'last');
    % And to the right.
    upperPeakIdx            = find(negPeakLoc > peakLoc(peakIdx),1,'first');
    lowerPeakLoc(peakIdx)   = negPeakLoc(lowerPeakIdx);
    upperPeakLoc(peakIdx)   = negPeakLoc(upperPeakIdx);
    lowerPeakVal(peakIdx)   = -negPeakVal(lowerPeakIdx);
    upperPeakVal(peakIdx)   = -negPeakVal(upperPeakIdx);
end


% ensure to be a column  vector
Fl = peakLoc(:);

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

% Method B): Spectrum Interpolation
F1lower = lowerPeakLoc(:);
F1upper = upperPeakLoc(:);
Neighwidth = ones(size(Fl)).*4;

% error message if periodicity of the interference frequency doesn't match the DFT length 
n = round(floor(nsamples .* Fl./Fs) * Fs./Fl);
if n ~= nsamples 
   ft_error('Spectrum interpolation requires that the data length fits complete cycles of the powerline frequency, e.g., exact multiples of 20 ms for a 50 Hz line frequency (sampling rate of 1000 Hz).');
end

% frequencies to interpolate
f2int = zeros(size([F1upper,F1lower]));
for i = 1:length(Fl)
    f2int(i,:) = [F1lower(i) F1upper(i)];
end
% frequencies used for interpolation
for i = 1:length(Neighwidth)
    f4int(i,:) = [f2int(i,1)-Neighwidth(i) f2int(i,:) f2int(i,2)+Neighwidth(i)];
end

data_fft = fft(dat,nsamples,2); % calculate fft to obtain spectrum that will be interpolated
frq = Fs*linspace(0,1,nsamples+1);

% interpolate 50Hz (and harmonics) amplitude in spectrum
for i = 1:length(Fl)
    smpl2int = nearest(frq,f2int(i,1)):nearest(frq,f2int(i,2)); % samples of frequencies that will be interpolated
    smpl4int = [(nearest(frq,f4int(i,1)):nearest(frq,f4int(i,2))-1),(nearest(frq,f4int(i,3))+1:nearest(frq,f4int(i,4)))]; % samples of neighbouring frequencies used to calculate the mean

    % new amplitude is calculated as the mean of the neighbouring frequencies
    mns4int= bsxfun(@times, ones(size(data_fft(:,smpl2int))), mean(abs(data_fft(:,smpl4int)),2));

    % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
    data_fft(:,smpl2int) = bsxfun(@times, exp(bsxfun(@times,angle(data_fft(:,smpl2int)),1i)), mns4int); 
end

% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
% to ensure a real valued signal after iFFT
filt = ifft(data_fft,[],2,'symmetric');
   
