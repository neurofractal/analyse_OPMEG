function [filt, filtFreq] = ft_dft_modelling_filter_auto(cfg, data)
% Function to find peaks in spectrum, model their shape and remove them 
% before transforming back into time domain. 
% 
% This function takes inspiration from the spectrum interpolation method 
% (Leske & Dalal, 2019, NeuroImage 189, 10.1016/j.neuroimage.2019.01.026).
% The aim of this function is similar, but intends to automatically detect
% peaks and their distributions so as to minimise the data removed. 
% 
% Overview of function
% 1)   Data transformed into the frequency domain via using Welch method.
% 2)   Peaks in the specified frequency range are identified and some
%       properties extracted about them.
% 3)   Either independently or as a pair, estimates are made about the
%       distribution of peaks. 
% 4)   For each channel, a better fit is made with some channel-wide
%       constraints (centre and max width (g) of peak)
% 5)   Depending on user input, these peak fits are either (1) subtracted 
%       from the original fft power output or (2) used to determine the
%       edges of each peak and to replace data within those edges with a
%       fitted slope. 
% 4)   The signal is transformed back into the time domain via inverse DFT
%       (iDFT). 
%
% Use as
%   [filt, filtFreq] = ft_preproc_dft_remove_gause(cfg, data)
%   where
%   data                    = Fieldtrip data structure with continuous
%                               data. Usually raw data.
%   cfg                     = Structure containing the following:
%   cfg.foi                 = [start end], freq range to remove peaks 
%                               between e.g. [40 60];
%   cfg.Neighwidth          = width (in Hz) of peaks to evaluate, e.g. 2;
%   cfg.minPeakDistance     = minimum distance in Hz between peaks.
%   cfg.method              = 'leaveSlope' will replace data modelled on
%                               the slope beneath the peak. 'subtractPeak'
%                               will subtract the model of the peak from
%                               the PSD. 'chopPeak' will interpolate across
%                               the edges of the peak.
% 	cfg.peakShape           = 'Gaussian' or 'Lorentzian'. Determines the
%                               distribution used to model the peak. 
%                               Lorentzian usually works best.
%   cfg.independentPeaks    = 'yes' or 'no'. Whether to model peaks 
%                               independently or together. Note, at the
%                               moment only two peaks may be modelled
%                               together. 
%   cfg.log                 = 'yes' or 'no'. Whether to log the PSD before
%                               fitting. May help with visualisation. Works
%                               best with 'leaveSlope' when 'yes' and best
%                               with 'removePeak' when 'no', for the most
%                               part.
%   cfg.strength            = multiplier of g, the variable that determines 
%                             distribution width. Higher strength will
%                             remove more of the peak, but increase the
%                             width of data lost.
%   cfg.windowLength        = Number of seconds defining the length of the
%                               window used in the Welch estimate. Adjust 
%                               until the estimate appears smooth. e.g. 20.
%                               Set to [] for pwelch default.
%   cfg.numPeaks            = Number of peaks to be modelled and removed.
% 
%   cfg.centreWeight        = weighting factor for the centre fitting.
% Copyright (C) 2020 Wellcome Trust Centre for Neuroimaging
% 
% Authors: Nicholas Alexander      (n.alexander@ucl.ac.uk) 
%          Oliver Alexander
%__________________________________________________________________________


%% Extract info from input structure
timeSeries          = data.trial{1};
nsamples            = length(timeSeries);
samplingFreq        = data.fsample;
fftFreq             = samplingFreq * linspace(0, 1, nsamples);
fftFreq             = fftFreq(1:nsamples/2+1);

% Run the fft
fftData             = fft(timeSeries,nsamples,2);
fftPow              = abs(fftData);

% To get a better estimation, use the Welch method
welchPow            = zeros(length(data.label),length(fftFreq));
windowLength        = samplingFreq*cfg.windowLength;
for chanIdx = 1:length(data.label)
    if chanIdx == 1
        [welchTmp, welchFreq]       = pwelch(timeSeries(chanIdx,:),windowLength,0,[],samplingFreq);
        welchPow(chanIdx,:)         = interp1(welchFreq,welchTmp,fftFreq);
    else
        welchTmp                    = pwelch(timeSeries(chanIdx,:),windowLength,0,[],samplingFreq);
        welchPow(chanIdx,:)         = interp1(welchFreq,welchTmp,fftFreq);
    end
end

% Put Welch in the same scale as the fft
window          = hanning(nsamples);
scaledWelchPow  = sqrt(welchPow*samplingFreq*(window'*window));
avgWelchPow     = median(scaledWelchPow,1);

% Tidy up
clear timeSeries samplingFreq nsamples welch* window*

% Select pow data for the foi.
foiIdx            	= ~(fftFreq < cfg.foi(1) | fftFreq > cfg.foi(2));
if strcmp(cfg.log,'yes')
    avgFftPowFOI           = log(avgWelchPow(foiIdx));
elseif strcmp(cfg.log,'no')
    avgFftPowFOI           = avgWelchPow(foiIdx);
end
freqFOI             = fftFreq(foiIdx);

%% Identify the requested number of peaks.
% Find all identifyable peaks first
if strcmp(cfg.independentPeaks,'yes')
    [~,~,~,allPeakProm]   = findpeaks(avgFftPowFOI,freqFOI,'MinPeakDistance',cfg.minPeakDistance);
elseif strcmp(cfg.independentPeaks,'no')
    [~,~,~,allPeakProm]   = findpeaks(avgFftPowFOI,freqFOI);
end

if isempty(allPeakProm)
    error('No peaks detected');
else
    allPeakProm         = sort(allPeakProm,'descend');
    minPeakProminence   = allPeakProm(cfg.numPeaks);
end

% Identify peaks based on prominence. Extract halfheight width.
if strcmp(cfg.independentPeaks,'yes')
    [~,peakFreq,peakWidth,peakProminence]  = findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'MinPeakDistance',cfg.minPeakDistance,'WidthReference','halfprom');
elseif strcmp(cfg.independentPeaks,'no')
    [~,peakFreq,peakWidth,peakProminence]  = findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'WidthReference','halfprom');
end


% Tidy up
clear minPeakProminence diffPeaks noPeaks chanIdx...
    initialPeakPromList freqFOI foiIdx smoothPowFOI current*...
    allPeakProm peakSortIdx peakCount avgFftPowFOI freqFOI...
    avgWelchPow


%% Define functions for fitting
% Define function for Gaussian/Lorentzian with and without
% slopes.
justSlope       = @(b, c, x)...
                    (b.*x) + c;
% Gaussian
distName{1}         = 'Gaussian';
peakWithSlope{1}   = @(A, x0, g, b, c, x)...
                    (A*exp(-2*((x-x0)/g).^2)) + (b*x) + c;
justPeak{1}     = @(A, x0, g, x)...
                    (A*exp(-2*((x-x0)/g).^2));

% Lorentzian
distName{2}         = 'Lorentzian';

peakWithSlope{2}   = @(A, x0, g, b, c, x)...
                    ((A .* (g.^2))./((g.^2) + ((x - x0).^2))) + (b.*x) + c;
justPeak{2}     = @(A, x0, g, x)...
                    ((A .* (g.^2))./((g.^2) + ((x - x0).^2)));

twoPeaksWithSlope{2} = @(A1, x01, g1, A2, x02, g2, b, c, x)...
                            ((A1 .* (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2))) + (b.*x) + c;

justTwoPeaks{2}      = @(A1, x01, g1, A2, x02, g2, x)...
                            ((A1 * (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2)));

% Pseudo voigt
voigtWithSlope  = @(m,A,x0,g,b,c,x)...
                        (1-m).*(A*exp(-2*((x-x0)/g).^2)) + m.*((A .* (g.^2))./((g.^2) + ((x - x0).^2))) + (b.*x + c);

justVoigt       = @(m,A,x0,g,x)...
                        (1-m).*(A*exp(-2*((x-x0)/g).^2)) + m.*((A .* (g.^2))./((g.^2) + ((x - x0).^2)));

%% Remove Gaussian component from amplitude in spectrum
switch cfg.independentPeaks
    case 'yes'
        filtFreq        = zeros(length(peakFreq),2);
        for peakIdx = 1:length(peakFreq)
            %% Find the neighbourhood data around the peak
            neighbourLowerBound         = peakFreq(peakIdx) - (cfg.Neighwidth);
            neighbourUpperBound         = peakFreq(peakIdx) + (cfg.Neighwidth);
            neighbourFreqIndicesBound   = nearest(fftFreq,[neighbourLowerBound,neighbourUpperBound]);
            neighbourFreqIndices        = neighbourFreqIndicesBound(1):neighbourFreqIndicesBound(end);
            
            if strcmp(cfg.log,'yes')
                % Get the chan avg data for those bounds.
                neighbourData               = log(scaledWelchPow(:,neighbourFreqIndices));
                fftNeighbourData            = log(fftPow(:,neighbourFreqIndices));
            elseif strcmp(cfg.log,'no')
                % Get the chan avg data for those bounds.
                neighbourData               = scaledWelchPow(:,neighbourFreqIndices);
                fftNeighbourData            = fftPow(:,neighbourFreqIndices);
            end

            %% First fit the peak with slope to each channel independently
            % Guesses for fit (peak level)
            % Amplitude
            A               = peakProminence(peakIdx);
            
            % Centre
            x0              = nearest(fftFreq,peakFreq(peakIdx));
            
            % Gamma (assumes fftFreq starts at 0)
            g               = nearest(fftFreq,peakWidth(peakIdx));
            
            
            % Get a better estimate of g and x0
            for chanIdx = 1:length(fftData(:,1))
                % Guesses for fit (channel level)
                neighbourChanData   = neighbourData(chanIdx,:);
                
                % Slope
                quarterLength   = round(length(neighbourChanData)/4);
                endVals         = mean(neighbourChanData((end - quarterLength):end));
                startVals       = mean(neighbourChanData(1:1+quarterLength));
                eighthLength    = round(quarterLength/2);
                b               = (endVals - startVals) / (neighbourFreqIndices(end - eighthLength) - neighbourFreqIndices(eighthLength));
                
                % Constant - What is c equal to for the slope go through the middle?
                middleValue     = mean([endVals,startVals]);
                c               = middleValue - (b * x0);
                
                % Weighting for the fit
                w               = ones(length(scaledWelchPow),1);
                w(x0-eighthLength:x0+eighthLength)    = cfg.centreWeight;
                w               = w(neighbourFreqIndicesBound(1):neighbourFreqIndicesBound(end));
                
%                 % Fit a peak with slope to the channel data.
%                 for fitIdx = 1:length(peakWithSlope)
%                     peakFitOptions              = fitoptions('method','Nonlinear', 'normal', 'off','weight',w);
%                     peakFitOptions.StartPoint   = [A, x0, g, b, c];
%                     peakFitOptions.Lower        = [0, x0-g, 0, -inf, -inf];
%                     peakFitOptions.Upper        = [inf, x0+g, inf, inf, inf];
%                     peakWithSlopeFitType        = fittype(peakWithSlope{fitIdx});
%                     [fittedModel,goodness]      = fit(neighbourFreqIndices', neighbourChanData', peakWithSlopeFitType, peakFitOptions);
%                     
%                     rmse(fitIdx)                = goodness.rmse;
%                     fitCoeffValue   = coeffvalues(fittedModel);
%                     fitCoeffNames   = coeffnames(fittedModel);
%                     
%                     % Dynamic naming seems to be the best way...
%                     for coeffIdx = 1:length(fitCoeffNames)
%                         fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx,fitIdx) = fitCoeffValue(coeffIdx);
%                     end
%                 end
                peakFitOptions              = fitoptions('method','Nonlinear', 'normal', 'off','weight',w);
                peakFitOptions.StartPoint   = [0.5, A, x0, g, b, c];
                peakFitOptions.Lower        = [-inf,0, x0-g, 0, -inf, -inf];
                peakFitOptions.Upper        = [inf,inf, x0+g, inf, inf, inf];
                peakWithSlopeFitType        = fittype(voigtWithSlope);
                [fittedModel,goodness]      = fit(neighbourFreqIndices', neighbourChanData', peakWithSlopeFitType, peakFitOptions);
%                     
                rmse(chanIdx)               = goodness.rmse;
                fitCoeffValue               = coeffvalues(fittedModel);
                fitCoeffNames               = coeffnames(fittedModel);

                % Dynamic naming seems to be the best way...
                for coeffIdx = 1:length(fitCoeffNames)
                    fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx) = fitCoeffValue(coeffIdx);
                end
            end
%             
%             [bestOverallDist,fTmp]    = mode(bestDist);
%             disp(strcat(distName{bestOverallDist},' is the best fitting distribution in',{' '},num2str(fTmp),' of',{' '},num2str(length(fftData(:,1))),' channels'));

            
            % Set the x0 for remaining fits
            x0              = median(fittedStruct.x0(:));
            
            %% Set some limits for the next fitting
            % Gamma should not be wider than the Welch estimate, but may be
            % narrower.
            g               = median(abs(fittedStruct.g(:)));
            iqrG            = iqr(abs(fittedStruct.g(:)));
            maxG            = g + 3*iqrG;
            
            % b should be about right from the Welch estimate so constrain.
            b               = median(fittedStruct.b(:));
            iqrB            = iqr(fittedStruct.b(:));
            maxB            = b + 3*iqrB;
            minB            = b - 3*iqrB;
            
            % A may vary greatly from channel to channel
            
            %% Fit again with new limits and the original fft data
            for chanIdx = 1:length(fftData(:,1))
%                 for fitIdx = 1:length(peakWithSlope)
%                     % Guesses for fit (channel level)
%                     fftNeighbourChanData   = fftNeighbourData(chanIdx,:);
% 
%                     peakFitOptions              = fitoptions('method','Nonlinear', 'normal', 'off','weight',w);
%                     peakFitOptions.StartPoint   = [fittedStruct.A(chanIdx,bestOverallDist), x0, g, b, fittedStruct.c(chanIdx,bestOverallDist)];
%                     peakFitOptions.Lower        = [0, x0, 0, minB, -inf];
%                     peakFitOptions.Upper        = [inf, x0, maxG, maxB, inf];
%                     peakWithSlopeFitType        = fittype(peakWithSlope{fitIdx});
%                     [fittedModel,goodness]                 = fit(neighbourFreqIndices', fftNeighbourChanData', peakWithSlope{bestOverallDist},peakFitOptions);                
%                     
%                     rmse(fitIdx)                = goodness.rmse;
% 
%                     fitCoeffValue   = coeffvalues(fittedModel);
%                     fitCoeffNames   = coeffnames(fittedModel);
% 
%                     % Dynamic naming seems to be the best way...
%                     for coeffIdx = 1:length(fitCoeffNames)
%                         fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx,bestOverallDist) = fitCoeffValue(coeffIdx);
%                     end
%                 end
%                 % Decide which distribution to use
%                 [~,fitOrder] = sort(rmse,'ascend');
%             
%                 bestDist(chanIdx)    = fitOrder(1);

                % Guesses for fit (channel level)
                fftNeighbourChanData        = fftNeighbourData(chanIdx,:);

                peakFitOptions              = fitoptions('method','Nonlinear', 'normal', 'off','weight',w);
                peakFitOptions.StartPoint   = [0.5, A, x0, g, b, c];
                peakFitOptions.Lower        = [-inf,0, x0-g, 0, -inf, -inf];
                peakFitOptions.Upper        = [inf,inf, x0+g, inf, inf, inf];
                peakWithSlopeFitType        = fittype(voigtWithSlope);
                [fittedModel,goodness]      = fit(neighbourFreqIndices', fftNeighbourChanData', peakWithSlopeFitType, peakFitOptions);
%                     
                rmse(chanIdx)               = goodness.rmse;
                fitCoeffValue               = coeffvalues(fittedModel);
                fitCoeffNames               = coeffnames(fittedModel);

                % Dynamic naming seems to be the best way...
                for coeffIdx = 1:length(fitCoeffNames)
                    fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx) = fitCoeffValue(coeffIdx);
                end
                
            end
%             
%          	[bestOverallDist,fTmp]    = mode(bestDist);
%             disp(strcat(distName{bestOverallDist},' is the best fitting distribution in',{' '},num2str(fTmp),' of',{' '},num2str(length(fftData(:,1))),' channels'));

            
            
            maxG            = max(abs(fittedStruct.g(:)));
            
            % Get the width of cfg.strength gamma of fitted lorentzian.
            tmpStartIdx             = nearest(neighbourFreqIndices,(x0 - (cfg.strength * maxG)));
            tmpEndIdx               = nearest(neighbourFreqIndices,(x0 + (cfg.strength * maxG)));

            indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
            
            % Save the freq limits that have been modified. 
            filtFreq(peakIdx,1:2)   = [fftFreq(indicesToReplace(1)),fftFreq(indicesToReplace(end))];

            if length(indicesToReplace) > length(neighbourData(chanIdx,:))
                ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
            end
            
            for chanIdx = 1:length(fftData(:,1))
                switch cfg.method
                    case 'subtractPeak'
                        % Get the peak on a slope
                        peakOnly                = justVoigt(fittedStruct.m(chanIdx),fittedStruct.A(chanIdx),fittedStruct.x0(chanIdx),fittedStruct.g(chanIdx),indicesToReplace);
                        
                        % Remove it from the original data.
                        replacementData         = fftNeighbourData(chanIdx,tmpStartIdx:tmpEndIdx) - peakOnly;

                    case 'chopPeak'
                        % Get the peak on a slope
                        fittedModelData         = voigtWithSlope(fittedStruct.m(chanIdx),fittedStruct.A(chanIdx),fittedStruct.x0(chanIdx),fittedStruct.g(chanIdx),fittedStruct.b(chanIdx),fittedStruct.c(chanIdx),neighbourFreqIndices);
                        tmpStartVal             = fittedModelData(tmpStartIdx);
                        tmpEndVal               = fittedModelData(tmpEndIdx);
                        
                        % Linear interpolation across the peak.
                        replacementData         = interp1([neighbourFreqIndices(tmpStartIdx),neighbourFreqIndices(tmpEndIdx)],[tmpStartVal,tmpEndVal],indicesToReplace,'linear');
                    case 'leaveSlope'
                        slopeOnly         = justSlope(fittedStruct.b(chanIdx),fittedStruct.c(chanIdx),neighbourFreqIndices);
                        replacementData     = slopeOnly(tmpStartIdx:tmpEndIdx);
                
                end
                
                if strcmp(cfg.log,'yes')
                    replacementData = exp(replacementData);
                end
                
                % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
                fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), replacementData);

            end
        end
    case 'no'
        filtFreq        = zeros(1,2);

        % Check there are 2 or more peaks
        if length(peakFreq) < 2
            error('Only one peak found. Please set cfg.independentPeaks to yes');
        end
        
        % Find the bounds for the width
        peakLowerBound              = peakFreq - peakWidth;
        peakUpperBound              = peakFreq + peakWidth;
        widthFreqIndicesBound       = nearest(fftFreq,[min(peakLowerBound),max(peakUpperBound)]);
        
        % And for the specified neighbourhood.
        neighbourLowerBound         = peakFreq(1) - cfg.Neighwidth;
        neighbourUpperBound         = peakFreq(end) + cfg.Neighwidth;
        neighbourFreqIndicesBound   = nearest(fftFreq,[neighbourLowerBound,neighbourUpperBound]);
        neighbourFreqIndices        = neighbourFreqIndicesBound(1):neighbourFreqIndicesBound(end);
        
        if strcmp(cfg.log,'yes')
            % Get the chan avg data for those bounds.
            neighbourData               = log(scaledWelchPow(:,neighbourFreqIndices));
        elseif strcmp(cfg.log,'no')
            % Get the chan avg data for those bounds.
            neighbourData               = scaledWelchPow(:,neighbourFreqIndices);
        end
               
        % Guesses for fit (channel level)
        A           = zeros(size(peakFreq));
        x0          = zeros(size(peakFreq));
        g           = zeros(size(peakFreq));
        for peakIdx = 1:length(peakFreq)
            % Guesses for fit (peak level)
            % Amplitude
            % Floor
            A(peakIdx)              = peakProminence(peakIdx);
            
            % Centre
            x0(peakIdx)             = nearest(fftFreq,peakFreq(peakIdx));

            % Gamma. Sqrt to account for smoothing of estimate.
            g(peakIdx)              = nearest(fftFreq,peakWidth(peakIdx));
        end

        
        for chanIdx = 1:length(fftData(:,1))
            % Guesses for fit (channel level)
            neighbourChanData   = neighbourData(chanIdx,:);
            % Slope
            quarterLength   = round(length(neighbourChanData)/4);
            endVals         = mean(neighbourChanData((end - quarterLength):end));
            startVals       = mean(neighbourChanData(1:1+quarterLength));
            eighthLength    = round(quarterLength/2);
            b               = (endVals - startVals) / (neighbourFreqIndices(end - eighthLength) - neighbourFreqIndices(eighthLength));

            % Constant - What is c equal to for the slope go through the middle?
            middleValue     = mean([endVals,startVals]);
            c               = middleValue - (b * mean(x0));

            % Fit two peaks with slope to the channel data.
            fittedModel     = fit(neighbourFreqIndices', neighbourChanData', twoPeaksWithSlope,...
                                'StartPoint', [A(1), x0(1), g(1), A(2), x0(2), g(2), b, c],...
                                'Lower',[-inf, x0(1)-g(1), -inf, -inf, x0(2)-g(2), -inf, -inf, -inf],...
                                'Upper',[inf, x0(1)+g(1), inf, inf, x0(2)+g(2), inf, inf, inf]);

            fitCoeffValue   = coeffvalues(fittedModel);
            fitCoeffNames   = coeffnames(fittedModel);

            % Dynamic naming seems to be the best way...
            for coeffIdx = 1:length(fitCoeffNames)
                fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx) = fitCoeffValue(coeffIdx);
            end
        end
        
        % Better estimate of x0
        x0(1)           = median(fittedStruct.x01);
        x0(2)           = median(fittedStruct.x02);
        
        
        % It might be worth using the maxG as a limit for the fitted
        % model. 
        medianG(1)        = median(abs(fittedStruct.g1));
        iqrG(1)            = iqr(abs(fittedStruct.g1));
        maxG(1)            = medianG(1) + 3*iqrG(1);
        minG(1)            = abs(medianG(1) - 3*iqrG(1));
        medianG(2)        = median(abs(fittedStruct.g2));
        iqrG(2)            = iqr(abs(fittedStruct.g2));
        maxG(2)            = medianG(2) + 3*iqrG(2);
        minG(2)            = abs(medianG(2) - 3*iqrG(2));
        
        for chanIdx = 1:length(fftData(:,1))
            % Guesses for fit (channel level)
            neighbourChanData   = neighbourData(chanIdx,:);
            % Slope
            quarterLength   = round(length(neighbourChanData)/4);
            endVals         = mean(neighbourChanData((end - quarterLength):end));
            startVals       = mean(neighbourChanData(1:1+quarterLength));
            eighthLength    = round(quarterLength/2);
            b               = (endVals - startVals) / (neighbourFreqIndices(end - eighthLength) - neighbourFreqIndices(eighthLength));

            % Constant - What is c equal to for the slope go through the middle?
            middleValue     = mean([endVals,startVals]);
            c               = middleValue - (b * mean(x0));

            % Fit two peaks with slope to the channel data.
            fittedModel     = fit(neighbourFreqIndices', neighbourChanData', twoPeaksWithSlope,...
                                'StartPoint', [A(1), x0(1), medianG(1), A(2), x0(2), medianG(2), b, c],...
                                'Lower',[-inf, x0(1), minG(1), -inf, x0(2), minG(2), -inf, -inf],...
                                'Upper',[inf, x0(1), maxG(1), inf, x0(2), maxG(2), inf, inf]);

            fitCoeffValue   = coeffvalues(fittedModel);
            fitCoeffNames   = coeffnames(fittedModel);

            % Dynamic naming seems to be the best way...
            for coeffIdx = 1:length(fitCoeffNames)
                fittedStruct.(fitCoeffNames{coeffIdx})(chanIdx) = fitCoeffValue(coeffIdx);
            end
        end
        
        maxG(1)            = max(abs(fittedStruct.g1));
        maxG(2)             = max(abs(fittedStruct.g2));
        maxG                = max(maxG);
        
        % Get the width of fitted lorentzians.
        x0Min                   = min([fittedStruct.x01(1);fittedStruct.x02(1)]);
        x0Max                   = max([fittedStruct.x01(1);fittedStruct.x02(1)]);
        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * maxG));
        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * maxG));

        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
        
        
        filtFreq(1,1:2)   = [fftFreq(indicesToReplace(1)),fftFreq(indicesToReplace(end))];

        if length(indicesToReplace) > length(neighbourData)
            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
        end
        
        for chanIdx = 1:length(fftData(:,1))
            switch cfg.method
                case 'removePeak'
                    % Get the lorentzian component.
                    peaksOnly             = justTwoPeaks(fittedStruct.A1(chanIdx),fittedStruct.x01(chanIdx),fittedStruct.g1(chanIdx),fittedStruct.A2(chanIdx),fittedStruct.x02(chanIdx),fittedStruct.g2(chanIdx),indicesToReplace);

                    % Remove it from the original data.
                    replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - peaksOnly;


                case 'leaveSlope'
                    % Get the peak on a slope
                    fittedModelData         = twoPeaksWithSlope(fittedStruct.A1(chanIdx),fittedStruct.x01(chanIdx),fittedStruct.g1(chanIdx),fittedStruct.A2(chanIdx),fittedStruct.x02(chanIdx),fittedStruct.g2(chanIdx),fittedStruct.b(chanIdx),fittedStruct.c(chanIdx),neighbourFreqIndices);
                    tmpStartVal             = fittedModelData(tmpStartIdx);
                    tmpEndVal               = fittedModelData(tmpEndIdx);

                    % Linear interpolation across the peak.
                    replacementData  = interp1([neighbourFreqIndices(tmpStartIdx),neighbourFreqIndices(tmpEndIdx)],[tmpStartVal,tmpEndVal],indicesToReplace,'linear');

            end
            
            if strcmp(cfg.log,'yes')
                replacementData = exp(replacementData);
            end
            
            % Eulers formula: replace noise components with new amplitude combined with phase, that is retained from the original data
            fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), replacementData);
        end
end

% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
% to ensure a real valued signal after iFFT
filteredTimeseries = ifft(fftData,[],2,'symmetric');

% Put it back into original structure.
filt            = data;
filt.trial{1}   = filteredTimeseries;

