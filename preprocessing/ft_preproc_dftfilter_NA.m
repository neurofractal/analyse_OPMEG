function [filt] = ft_preproc_dftfilter_NA(cfg, data)

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

%% Extract info from input structure
timeSeries          = data.trial{1};
nsamples            = data.sampleinfo(end) - data.sampleinfo(1) - 1;
samplingFreq        = data.fsample;
fftFreq             = samplingFreq * linspace(0, 1, nsamples + 1);

%% Identify peaks in power in PSD
% Calculate PSD in smaller epochs. 
cfgPSD              = [];
cfgPSD.channel      = 'all';
cfgPSD.trial_length = 3;
cfgPSD.method       = 'tim';
cfgPSD.foi          = cfg.foi;
cfgPSD.plot         = 'no';
[pow,freq,~]        = ft_opm_psd(cfgPSD,data);

% Use those to produce a smoother PSD output for peak identification.
powMed              = median(pow(:,:,:),3);
powMed              = median(pow(:,:),2);

clear pow cfgPSD % tidy up

% The user will be asked if they are satisfied with peak identification.
changePeaks         = true;
minPeakProminence   = 16; % Starting point for MinPeakPriminence method.
while changePeaks
    % Identify peaks based on prominence. Extract halfheight width.
    [peakVal,peakFreq,halfHeightWidth]  = findpeaks(powMed,freq,'MinPeakProminence',minPeakProminence,'WidthReference','halfheight','MinPeakDistance',cfg.minPeakDistance);
    
    % Remove peaks outside foi.
    peaksToRemove           = (peakFreq < cfg.foi(1) | peakFreq > cfg.foi(2));
    peakVal(peaksToRemove)  = [];
    peakFreq(peaksToRemove) = [];

    % Find the halfheight frequencies and interp onto fftFreq and get
    % samples (samples in freq space)
    threeStdDevFftFreqRange = zeros(2, length(peakFreq));
    threeStdDevFreqRange    = zeros(2, length(peakFreq));
    neighbourFreqRange      = zeros(2, length(peakFreq));
    neighbourFftFreqRange   = zeros(2, length(peakFreq));
    wholePeakPow            = cell(size(peakFreq));
    wholePeakFreq           = cell(size(peakFreq));
    for peakIdx = 1:length(peakFreq)        
        % halfHeightWidth (or FWHM) = 2.355 std.dev. = 2sqrtln2 std.dev.
        stdDevWidth                        = halfHeightWidth(peakIdx) / 2.355;
        threeStdDevFreqRange(1)    = peakFreq(peakIdx) - (5 * stdDevWidth);
        threeStdDevFreqRange(2)    = peakFreq(peakIdx) + (5 * stdDevWidth);
        
        neighbourFreqRange(1,peakIdx)    = peakFreq(peakIdx) - cfg.Neighwidth;
        neighbourFreqRange(2,peakIdx)    = peakFreq(peakIdx) + cfg.Neighwidth;
        
        % interpret onto fftFreq
        % PeakFreq
        tmp                         = abs(fftFreq - peakFreq(peakIdx));
        [~, ~, tmp2]                = unique(tmp);
        tmp3                        = find(tmp2 == 1, 1);
        peakFreq(peakIdx)           = fftFreq(tmp3); % Take the lower freq if disputed.

        for i = 1:2
            tmp                             = abs(fftFreq - threeStdDevFreqRange(i));
            [~, ~, tmp2]                    = unique(tmp);
            tmp3                            = find(tmp2 == 1, 1);
            threeStdDevFftFreqRange(i,peakIdx) = fftFreq(tmp3);
            
            tmp                             = abs(fftFreq - neighbourFreqRange(i));
            [~, ~, tmp2]                    = unique(tmp);
            tmp3                            = find(tmp2 == 1, 1);
            neighbourFftFreqRange(i,peakIdx) = fftFreq(tmp3);
        end
        
        % For the plot, get the line for each peak
        startIdx                = nearest(freq,neighbourFftFreqRange(1,peakIdx));
        endIdx                  = nearest(freq,neighbourFftFreqRange(2,peakIdx));
        wholePeakPow{peakIdx}   = powMed(startIdx:endIdx);
        wholePeakFreq{peakIdx}  = freq(startIdx:endIdx);
    end
    
    %% Plot result for verification
    hold on
    verificationPlot    = gcf;
    
    % Plot the peaks
    for peakIdx = 1:length(wholePeakPow)
        plot(wholePeakFreq{peakIdx},wholePeakPow{peakIdx},'LineWidth',2,'Color','b');
        scatter(peakFreq(peakIdx),peakVal(peakIdx),'r');
    end

    % Plot the PSD (from ft_plot_PSD)
    plot(freq,powMed,'LineWidth',1,'Color','k');
    
    % Setup axes
    ax                  = gca;
    ax.FontSize         = 16;
    ax.TickLength       = [0.02 0.02];
    set(ax,'yscale','log')
    labY                = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
    ylabel(labY,'interpreter','latex','FontSize',30)
    xlabel('Frequency (Hz)','FontSize',30)
    xlim([cfg.foi(1), cfg.foi(end)]);

    
    % Ask if the user wants to change the number of peaks
    diffPeaks           = 3;
    while ~any(eq(diffPeaks,[0 1 2]))
        diffPeaks = input('To continue, enter 0\nFor less peaks, enter 1\nFor more peaks, enter 2:\n');
    end
    
    % There is probably a better way, but this works for now.
    if diffPeaks == 1
        minPeakProminence = minPeakProminence * 2;
        close(gcf);
    elseif diffPeaks == 2
        minPeakProminence = minPeakProminence / 2;
        close(gcf);
    else
        changePeaks = 0;
    end
end

% Tidy up
close(verificationPlot)
clear tmp* *Tmp peakIdx peaksToRemove diffPeaks changePeaks ax

%% Checks from interp method still apply.
% preprocessing fails on channels that contain NaN
if any(isnan(timeSeries(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'Time series data contains NaN values');
end

% error message if periodicity of the interference frequency doesn't match the DFT length 
% This error should never be able to happen because I have interpolated
% earlier
n = round(floor(nsamples .* peakFreq./samplingFreq) * samplingFreq./peakFreq);
if n ~= nsamples 
   ft_error('Spectrum interpolation requires that the data length fits complete cycles of the frequencies being interpolated, e.g., exact multiples of 20 ms for a 50 Hz line frequency (sampling rate of 1000 Hz).');
end

%% Remove Gaussian component from amplitude in spectrum
% Run the fft
fftData = fft(timeSeries,nsamples,2); % calculate fft to obtain spectrum that will be interpolated

for peakIdx = 1:length(peakFreq)
    % Find the indices for the frequencies of the peak.
    threeStdDevFreqIndices  = nearest(fftFreq,threeStdDevFftFreqRange(1,peakIdx)):nearest(fftFreq,threeStdDevFftFreqRange(2,peakIdx));
    neighbourFreqIndices    = nearest(fftFreq,neighbourFftFreqRange(1,peakIdx)):nearest(fftFreq,neighbourFftFreqRange(2,peakIdx));
    
    % Guesses for fit (all channels)
    % Gauss centre freq
    centre      = mean(threeStdDevFreqIndices);
    
    % Gauss Width
    widthSamples = length(threeStdDevFreqIndices);
    
    % Define function for a Gaussian on a slope.
    gaussWithSlope = @(A, x0, s, b, c, x)...
        (A*exp(-2*((x-x0)/s).^2)) + (b*x) + c;
    
    % Define function for Gaussian (without slope).
    gauss = @(A, x0, s, x)...
        (A*exp(-2*((x-x0)/s).^2));
    
    slope = @(b, c, x)...
        (b*x) + c;
    
    
    % Fit channels independently (may need to revise this)
    for chanIdx = 1:length(fftData(:,1))
        % Get the data being fit
    	dataToFit   = abs(fftData(chanIdx,neighbourFreqIndices));

        % Guesses for fit (channel level)
        % Floor
        fitFloor    = mean(abs(fftData(chanIdx,threeStdDevFreqIndices)),2);

        % Guass peak
        peak        = max(abs(fftData(chanIdx,threeStdDevFreqIndices)),[],2);

        % Slope
        slopeC       = (abs(fftData(chanIdx,neighbourFreqIndices(end))) - abs(fftData(chanIdx,neighbourFreqIndices(1)))) / (threeStdDevFreqIndices(end) - threeStdDevFreqIndices(1));

        % Constant
        constant    = abs(fftData(chanIdx,neighbourFreqIndices(1))) - (slopeC *  neighbourFreqIndices(1));

        % Fit a gauss with slope to the data.
        fittedGaussWithSlope    = fit(neighbourFreqIndices', dataToFit', gaussWithSlope,...
            'StartPoint', [(peak - fitFloor), centre, widthSamples, slopeC, constant]);

        % Find the Gauss component.
        onlyGauss   = gauss(fittedGaussWithSlope.A,fittedGaussWithSlope.x0, fittedGaussWithSlope.s, neighbourFreqIndices);
        
        % Get the width of 2 std.dev of the fitted Gauss
        tmpStartIdx             = nearest(neighbourFreqIndices,fittedGaussWithSlope.x0 - (5 * fittedGaussWithSlope.s));
        tmpEndIdx               = nearest(neighbourFreqIndices,fittedGaussWithSlope.x0 + (5 * fittedGaussWithSlope.s));

        fittedWidthFreqIndices  = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

        % Get just the slope
        onlySlope   = slope(fittedGaussWithSlope.b,fittedGaussWithSlope.c, fittedWidthFreqIndices);
        
        % Debugging:
        % Use the output parameters to find the fit
        fitted      = gaussWithSlope(fittedGaussWithSlope.A, fittedGaussWithSlope.x0, fittedGaussWithSlope.s, fittedGaussWithSlope.b, fittedGaussWithSlope.c, neighbourFreqIndices);

        % Best guess (something is wrong with the guess, but it is okay for
        % now).
        bestGuess   = gaussWithSlope((peak - fitFloor), centre, widthSamples, slopeC, constant, neighbourFreqIndices);

        % Plots
        hold on
%         plot(neighbourFreqIndices,bestGuess);
        plot(neighbourFreqIndices,fitted);

        plot(neighbourFreqIndices,dataToFit);
%         plot(fittedWidthFreqIndices,onlySlope);
        
        % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
        fftData(chanIdx,fittedWidthFreqIndices) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,fittedWidthFreqIndices)),1i)), onlySlope); 
    end
end

% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
% to ensure a real valued signal after iFFT
filt = ifft(fftData,[],2,'symmetric');
   

