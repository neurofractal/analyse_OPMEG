function [filt] = ft_preproc_dft_remove_gauss(cfg, data)
% Function to peaks in spectrum, model their shape and remove them before
% transforming back into time domain. 
% 
% Function is based on Search Results ft_preproc_dftfilter. In particular,
% the spectrum interpolation method (Leske & Dalal, 2019, NeuroImage 189,
% doi: 10.1016/j.neuroimage.2019.01.026).
% 
% The signal is:
% 1)   transformed into the frequency domain via a discrete Fourier 
%       transform (DFT), 
% 2)   Peaks in the specified frequency range are identified and some
%       characteristics extracted about them.
% 3)   Each peak is modelled as a Guassian on a slope. The Gaussian part is
%       replaced by the slope
% 4)   The signal is transformed back into the time domain via inverse DFT
%       (iDFT).
%
% Use as
%   [filt] = ft_preproc_dft_remove_gause(cfg, data)
% where
%   data            Fieldtrip data structure with continuous data.
%   cfg             Structure containing the following:
%   cfg.foi                 = [start end], freq range to remove peaks between e.g. [40 60];
%   cfg.Neighwidth          = width (in Hz) of peaks to evaluate, e.g. 2;
%   cfg.minPeakDistance     = minimum distance between peaks. minimum of 2 Neighwidth recommended.
%   cfg.method              = 'leaveSlope' or 'removeGauss'

%% Extract info from input structure
timeSeries          = data.trial{1};
nsamples            = length(timeSeries);
samplingFreq        = data.fsample;
fftFreq             = samplingFreq * linspace(0, 1, nsamples);

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
powMedLog           = log(powMed);
clear pow cfgPSD % tidy up

% The user will be asked if they are satisfied with peak identification.
changePeaks         = true;
minPeakProminence   = 16; % Starting point for MinPeakPriminence method.
while changePeaks
    noPeaks     = true;
    while noPeaks
        % Identify peaks based on prominence. Extract halfheight width.
        [peakVal,peakFreq,halfHeightWidth]  = findpeaks(powMedLog,freq,'MinPeakProminence',minPeakProminence,'WidthReference','halfheight','MinPeakDistance',cfg.minPeakDistance);
        % Remove peaks outside foi.
        peaksToRemove           = (peakFreq < cfg.foi(1) | peakFreq > cfg.foi(2));
        peakVal(peaksToRemove)  = [];
        peakFreq(peaksToRemove) = [];

        if isempty(peakVal)
            minPeakProminence = minPeakProminence / 2;
        else
            noPeaks = false;
        end
    end

    % Find the halfheight frequencies and interp onto fftFreq and get
    % samples (samples in freq space)
    threeStdDevFftFreqRange = zeros(2, length(peakFreq));
    neighbourFftFreqRange   = zeros(2, length(peakFreq));
    peakFftFreq             = zeros(1, length(peakFreq));
    wholePeakPow            = cell(size(peakFreq));
    wholePeakFreq           = cell(size(peakFreq));
    for peakIdx = 1:length(peakFreq)        
        % halfHeightWidth (or FWHM) = 2.355 std.dev. = 2sqrtln2 std.dev.
        stdDevWidth                = halfHeightWidth(peakIdx) / 2.355;
        threeStdDevFreqRange(1)    = peakFreq(peakIdx) - (5 * stdDevWidth);
        threeStdDevFreqRange(2)    = peakFreq(peakIdx) + (5 * stdDevWidth);
        
        neighbourFreqRange(1)    = peakFreq(peakIdx) - cfg.Neighwidth;
        neighbourFreqRange(2)    = peakFreq(peakIdx) + cfg.Neighwidth;
        
        % interpret onto fftFreq
        % PeakFreq
        tmp                         = abs(fftFreq - peakFreq(peakIdx));
        [~, ~, tmp2]                = unique(tmp);
        tmp3                        = find(tmp2 == 1, 1);
        peakFftFreq(peakIdx)        = fftFreq(tmp3(1)); % Take the lower freq if disputed.

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
        scatter(peakFreq(peakIdx),exp(peakVal(peakIdx)),'r');
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
% 
% %% Checks from interp method still apply.
% % preprocessing fails on channels that contain NaN
% if any(isnan(timeSeries(:)))
%   ft_warning('FieldTrip:dataContainsNaN', 'Time series data contains NaN values');
% end
% 
% % error message if periodicity of the interference frequency doesn't match the DFT length 
% % This error should never be able to happen because I have interpolated
% % earlier
% n = round(floor(nsamples .* peakFftFreq./samplingFreq) * samplingFreq./peakFftFreq);
% if n ~= nsamples 
%    ft_error('Spectrum interpolation requires that the data length fits complete cycles of the frequencies being interpolated, e.g., exact multiples of 20 ms for a 50 Hz line frequency (sampling rate of 1000 Hz).');
% end

%% Remove Gaussian component from amplitude in spectrum
% Run the fft
fftData = fft(timeSeries,nsamples,2); % calculate fft to obtain spectrum that will be interpolated

for peakIdx = 1:length(peakFftFreq)
    % Find the indices for the frequencies of the peak.
    threeStdDevFreqIndices  = nearest(fftFreq,threeStdDevFftFreqRange(1,peakIdx)):nearest(fftFreq,threeStdDevFftFreqRange(2,peakIdx));
    neighbourFreqIndices    = nearest(fftFreq,neighbourFftFreqRange(1,peakIdx)):nearest(fftFreq,neighbourFftFreqRange(2,peakIdx));
    
    % Guesses for fit (all channels)
    % Gauss Width
    stdDevSamples = length(threeStdDevFreqIndices) / 6;
    
    % Define function for a Gaussian on a slope.
    gaussWithSlope = @(A, x0, s, b, c, x)...
        (A*exp(-2*((x-x0)/s).^2)) + (b*x) + c;
    
    % Define function for Gaussian (without slope).
    gauss = @(A, x0, s, x)...
        (A*exp(-2*((x-x0)/s).^2));
    
    slope = @(b, c, x)...
        (b*x) + c;
    
    lorentzianWithSlope = @(A, x0, g, b, c, x)...
        ((A * (g^2))./((g^2) + ((x - x0).^2))) + (b*x) + c;
    
    twoLortenzianWithSlope = @(A1, A2, x01, x02, g1, g2, b, c, x)...
        ((A1 * (g1^2))./((g1^2) + ((x - x01).^2))) + ((A2 * (g2^2))./((g2^2) + ((x - x02).^2))) + (b*x) + c;
    
    lorentzian = @(A, x0, g, x)...
        ((A * (g^2))./((g.^2) + ((x - x0).^2)));

    twoLorentzian = @(A1, A2, x01, x02, g1, g2, x)...
        ((A1 * (g1^2))./((g1^2) + ((x - x01).^2))) + ((A2 * (g2^2))./((g2^2) + ((x - x02).^2)));
    
    % Fit channels independently (may need to revise this)
    for chanIdx = 1:length(fftData(:,1))
        % Get the data being fit
    	neighbourData   = log(abs(fftData(chanIdx,neighbourFreqIndices)));
        threeStdDevData = log(abs(fftData(chanIdx,threeStdDevFreqIndices)));
%         dataToFit   = sgolayfilt(dataToFit,10,21);
        
        % Guesses for fit (channel level)
        % Floor
        fitFloor    = mean(threeStdDevData,2);

        % Guass peak
        [peak,peakLoc]    = max(threeStdDevData,[],2);
        
        % Gauss centre freq
        centre = threeStdDevFreqIndices(peakLoc);

        % Slope
        slopeC       = (neighbourData(end) - neighbourData(1)) / (threeStdDevFreqIndices(end) - threeStdDevFreqIndices(1));

        % Constant
        constant    = neighbourData(1) - (slopeC *  neighbourFreqIndices(1));
        
        switch cfg.method
            case 'leaveSlopeMinusGauss'
                % Fit a gauss with slope to the data.
                fittedModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlope,...
                                            'StartPoint', [(peak - fitFloor), centre, stdDevSamples, slopeC, constant]);
                
                % Get the width of 3 std.dev of the fitted Gauss
                tmpStartIdx             = nearest(neighbourFreqIndices,fittedModel.x0 - (3 * fittedModel.s));
                tmpEndIdx               = nearest(neighbourFreqIndices,fittedModel.x0 + (3 * fittedModel.s));
                indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                
                % Get just the slope
                replacementData         = slope(fittedModel.b,fittedModel.c, indicesToReplace);
                
                
            case 'removeGauss'
                % Fit a gauss with slope to the data.
                fittedModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlope,...
                                            'StartPoint', [(peak - fitFloor), centre, stdDevSamples, slopeC, constant]);
                
                % Get the width of 3 std.dev of the fitted Gauss
                tmpStartIdx             = nearest(neighbourFreqIndices,fittedModel.x0 - (3 * fittedModel.s));
                tmpEndIdx               = nearest(neighbourFreqIndices,fittedModel.x0 + (3 * fittedModel.s));
                indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                % Find the Gauss component.
                onlyGauss               = gauss(fittedModel.A,fittedModel.x0, fittedModel.s, indicesToReplace);
                
                % Remove the Gauss from the original data.
                replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyGauss;
                
                
                
            case 'removeLorentzian'
                % Fit a Lorentzian with slope
                fittedModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlope,...
                    'StartPoint', [(peak - fitFloor), centre, stdDevSamples, slopeC, constant]);
        
                % Get the width of 6 gamma of fitted lorentzian.
                tmpStartIdx             = nearest(neighbourFreqIndices,fittedModel.x0 - (6 * fittedModel.g));
                tmpEndIdx               = nearest(neighbourFreqIndices,fittedModel.x0 + (6 * fittedModel.g));

                indicesToReplace  = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                
                if length(indicesToReplace) > length(neighbourData)
                    ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                end
                
                % Get the lorentzian component.
                onlyLorentz = lorentzian(fittedModel.A, fittedModel.x0, fittedModel.g, indicesToReplace);

                % Remove it from the original data.
                replacementData        = neighbourData(tmpStartIdx:tmpEndIdx) - onlyLorentz;
                
%                 fitted2     = lorentzianWithSlope(fittedLorentzianWithSlope.A, fittedLorentzianWithSlope.x0, fittedLorentzianWithSlope.g,fittedLorentzianWithSlope.b,fittedLorentzianWithSlope.c, neighbourFreqIndices);
            case 'leaveSlopeMinusLorentzian'
                % Fit a lortenzian with slope to the data.
                fittedModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlope,...
                    'StartPoint', [(peak - fitFloor), centre, stdDevSamples, slopeC, constant]);
               
                % Get the width of 3 std.dev of the fitted Gauss
                tmpStartIdx             = nearest(neighbourFreqIndices,fittedModel.x0 - (6 * fittedModel.g));
                tmpEndIdx               = nearest(neighbourFreqIndices,fittedModel.x0 + (6 * fittedModel.g));
                indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                
                % Get just the slope
                replacementData         = slope(fittedModel.b,fittedModel.c, indicesToReplace);
        end
        
        % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
        fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), exp(replacementData));


        
%         fittedTwoLortenzianWithSlope = fit(neighbourFreqIndices', neighbourData', twoLortenzianWithSlope,...
%             'StartPoint', [(peak - fitFloor), (16.87 - fitFloor), centre, 112671, stdDevSamples, stdDevSamples/3, slopeC, constant]);
        
        
       
        
       
        
        % Debugging:
        % Use the output parameters to find the fit

        
%         fitted3     = twoLortenzianWithSlope(fittedTwoLortenzianWithSlope.A1, fittedTwoLortenzianWithSlope.A2,fittedTwoLortenzianWithSlope.x01, fittedTwoLortenzianWithSlope.x02,fittedTwoLortenzianWithSlope.g1,fittedTwoLortenzianWithSlope.g2,fittedTwoLortenzianWithSlope.b,fittedTwoLortenzianWithSlope.c,neighbourFreqIndices);
        % Best guess (something is wrong with the guess, but it is okay for
        % now).
%         bestGuess   = gaussWithSlope((peak - fitFloor), centre, stdDevSamples, slopeC, constant, neighbourFreqIndices);
        
%         
%         onlyTwoLorentz = twoLorentzian(fittedTwoLortenzianWithSlope.A1, fittedTwoLortenzianWithSlope.A2, fittedTwoLortenzianWithSlope.x01, fittedTwoLortenzianWithSlope.x02, fittedTwoLortenzianWithSlope.g1, fittedTwoLortenzianWithSlope.g2, neighbourFreqIndices);
% 
%         
%         
%         dataMinusTwoLorentz = neighbourData - onlyTwoLorentz;
% 
%         % Plots
%         hold on
%         plot(neighbourFreqIndices,bestGuess);
%         plot(neighbourFreqIndices,dataMinusGauss);
%         plot(neighbourFreqIndices,fitted);
%         plot(neighbourFreqIndices,fitted2);
%         plot(fittedWidthFreqIndices,onlySlope);
%         plot(neighbourFreqIndices,dataMinusLorentz);
%         plot(neighbourFreqIndices,fitted3);
%         plot(neighbourFreqIndices,dataMinusTwoLorentz);

%         close all
        
        
    end
end

% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
% to ensure a real valued signal after iFFT
filteredTimeseries = ifft(fftData,[],2,'symmetric');

% Put it back into original structure.
filt            = data;
filt.trial{1}   = filteredTimeseries;

