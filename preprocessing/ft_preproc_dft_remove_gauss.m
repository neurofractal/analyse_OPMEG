function [filt] = ft_preproc_dft_remove_gauss(cfg, data)
% Function to find peaks in spectrum, model their shape and remove them 
% before transforming back into time domain. 
% 
% Function is based on Search Results ft_preproc_dftfilter. In particular,
% the spectrum interpolation method (Leske & Dalal, 2019, NeuroImage 189,
% doi: 10.1016/j.neuroimage.2019.01.026).
% 
% Overview of function
% 1)   Data transformed into the frequency domain via a discrete Fourier 
%       transform (DFT). 
% 2)   Peaks in the specified frequency range are identified and some
%       characteristics extracted about them.
% 3)   The characteristics are used to model them onto data from each
%       channel.
% 4)   A specified method is used to remove the peaks.
% 4)   The signal is transformed back into the time domain via inverse DFT
%       (iDFT). 
%
% Use as
%   [filt] = ft_preproc_dft_remove_gause(cfg, data)
%   where
%   data                    = Fieldtrip data structure with continuous
%                               data. Usually raw data.
%   cfg                     = Structure containing the following:
%   cfg.foi                 = [start end], freq range to remove peaks 
%                               between e.g. [40 60];
%   cfg.Neighwidth          = width (in Hz) of peaks to evaluate, e.g. 2;
%   cfg.minPeakDistance     = minimum distance in Hz between peaks.
%   cfg.method              = 'leaveSlopeMinusGauss', 'removeGauss' 
%                               'removeLorentzian', or 
%                               'leaveSlopeMinusLorentzian'
%   cfg.independentPeaks    = 'yes' or 'no'. Whether to model peaks 
%                               independently or together.
%   cfg.log                 = 'yes' or 'no'. Whether to log the PSD before
%                               fitting.
%   cfg.strength            = g or s width to remove. ELABORATE
%   cfg.trialLength         = epoch length in seconds. Determines the
%                               smoothness of data for peak identifying.

%% Extract info from input structure
timeSeries          = data.trial{1};
nsamples            = length(timeSeries);
samplingFreq        = data.fsample;
fftFreq             = samplingFreq * linspace(0, 1, nsamples);

% Run the fft
fftData             = fft(timeSeries,nsamples,2);

% Get pow
fftPow              = abs(fftData);
avgFftPow           = median(fftPow,1);

% Tidy up
clear timeSeries samplingFreq nsamples 

% Select pow data for the foi.
foiIdx            	= ~(fftFreq < cfg.foi(1) | fftFreq > cfg.foi(2));
avgFftPowFOI        = avgFftPow(foiIdx);
if strcmp(cfg.log,'yes')
    avgFftPowFOI           = log(avgFftPow(foiIdx));
elseif strcmp(cfg.log,'no')
    avgFftPowFOI           = avgFftPow(foiIdx);
end
freqFOI             = fftFreq(foiIdx);

% The user will be asked if they are satisfied with peak identification.
firstRun            = true;
peakCount           = 0;
while firstRun || changePeaks
    % Get a staring point for peak prominence in the first run.
    if firstRun
        if strcmp(cfg.independentPeaks,'yes')
            [~,~,~,allPeakProm]   = findpeaks(avgFftPowFOI,freqFOI,'MinPeakDistance',cfg.minPeakDistance);
        elseif strcmp(cfg.independentPeaks,'no')
            [~,~,~,allPeakProm]   = findpeaks(avgFftPowFOI,freqFOI);
        end

        if isempty(allPeakProm)
            error('No peaks detected');
        else
            allPeakProm         = sort(allPeakProm,'descend');
            peakCount           = 1;
            minPeakProminence   = allPeakProm(peakCount);
        end
        
        % Prepare figure on first run
        figure;
    end
    
    % Identify peaks based on prominence. Extract halfheight width.
    if strcmp(cfg.independentPeaks,'yes')
        [~,peakFreq,peakWidth,peakProminence]  = findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'MinPeakDistance',cfg.minPeakDistance,'WidthReference','halfprom');
    elseif strcmp(cfg.independentPeaks,'no')
        [~,peakFreq,peakWidth,peakProminence]  = findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'WidthReference','halfprom');
    end

    % Show user the peaks
    if strcmp(cfg.independentPeaks,'yes')
        findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'Annotate','extents','MinPeakDistance',cfg.minPeakDistance,'WidthReference','halfprom');
    elseif strcmp(cfg.independentPeaks,'no')
        findpeaks(avgFftPowFOI,freqFOI,'MinPeakProminence',minPeakProminence,'Annotate','extents','WidthReference','halfprom');
    end
    
    % Ask if the user wants to change the number of peaks
    diffPeaks           = 3;
    while ~any(eq(diffPeaks,[0 1 2]))
        diffPeaks = input('To continue, enter 0\nFor less peaks, enter 1\nFor more peaks, enter 2:\n');
    end

    % Set minPeakProminence based on user input
    if diffPeaks == 1
        if minPeakProminence == max(allPeakProm)
            disp('Already at minimum peaks. Enter either 2 or 0');
            changePeaks     = true;
        else
            peakCount           = peakCount - 1;
            minPeakProminence   = allPeakProm(peakCount);
            changePeaks     = true;
        end
    elseif diffPeaks == 2
        if peakCount == length(allPeakProm)
            disp('No more peaks detected. Consider changing FOI')
            changePeaks     = true;
        else
            peakCount           = peakCount + 1;
            minPeakProminence   = allPeakProm(peakCount);
            changePeaks     = true;
        end
    else
        changePeaks = false;
        close(gcf);
        
        % Sort peaks by prominence
        [peakProminence,peakSortIdx]    = sort(peakProminence,'descend');
        peakFreq                        = peakFreq(peakSortIdx);
        peakWidth                       = peakWidth(peakSortIdx);
    end

    if firstRun
        firstRun = false;
    end
end


% Tidy up
clear firstRun changePeaks minPeakProminence diffPeaks noPeaks...
    initialPeakPromList freqFOI foiIdx smoothPowFOI...
    allPeakProm peakSortIdx peakCount avgFftPowFOI freqFOI

%% Remove Gaussian component from amplitude in spectrum
switch cfg.independentPeaks
    case 'yes'
        for peakIdx = 1:length(peakFreq)
            % Find the neighbourhood indices
            neighbourLowerBound         = peakFreq(peakIdx) - (cfg.Neighwidth);
            neighbourUpperBound         = peakFreq(peakIdx) + (cfg.Neighwidth);
            neighbourFreqIndicesBound   = nearest(fftFreq,[neighbourLowerBound,neighbourUpperBound]);
            neighbourFreqIndices        = neighbourFreqIndicesBound(1):neighbourFreqIndicesBound(end);
            
%             % And the frequencies
%             neighbourFreq               = fftFreq(neighbourFreqIndices);
            
            if strcmp(cfg.log,'yes')
                % Get the chan avg data for those bounds.
                neighbourData               = log(avgFftPow(neighbourFreqIndices));
            elseif strcmp(cfg.log,'no')
                % Get the chan avg data for those bounds.
                neighbourData               = avgFftPow(neighbourFreqIndices);
            end

            % Guesses for fit (peak level)
            % Amplitude
            A               = peakProminence(peakIdx);
            
            % Centre
            x0              = nearest(fftFreq,peakFreq(peakIdx));
            
            % Gamma (assumes fftFreq starts at 0)
            g               = nearest(fftFreq,peakWidth(peakIdx));
            
            % Guesses for fit (channel level)
            % Slope
            quarterLength   = round(length(neighbourData)/4);
            endVals         = mean(neighbourData((end - quarterLength):end));
            startVals       = mean(neighbourData(1:1+quarterLength));
            eighthLength    = round(quarterLength/2);
            b               = (endVals - startVals) / (neighbourFreqIndices(end - eighthLength) - neighbourFreqIndices(eighthLength));

            % Constant
            % What is c equal to for the slope go through the middle?
            middleValue     = mean([endVals,startVals]);
            c               = middleValue - (b * x0);
            
            % Define function for Gaussian/Lorentzian with and without
            % slopes.
            gaussWithSlope = @(A, x0, g, b, c, x)...
                (A*exp(-2*((x-x0)/g).^2)) + (b*x) + c;

            gauss = @(A, x0, g, x)...
                (A*exp(-2*((x-x0)/g).^2));

            slope = @(b, c, x)...
                (b.*x) + c;

            lorentzianWithSlope = @(A, x0, g, b, c, x)...
                ((A .* (g.^2))./((g.^2) + ((x - x0).^2))) + (b.*x) + c;

            lorentzian = @(A, x0, g, x)...
                ((A .* (g.^2))./((g.^2) + ((x - x0).^2)));
            
            switch cfg.method
                case 'leaveSlopeMinusGauss'
                    % Fit a gauss with slope to the data.
                    fittedAvgModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlope,...
                                                'StartPoint', [A, x0, g, b, c]);

                    bestGuess           = gaussWithSlope(A, x0, g, b, c, neighbourFreqIndices');
                    
                case 'removeGauss'
                    % Fit a gauss with slope to the data.
                    fittedAvgModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlope,...
                                                'StartPoint', [A, x0, g, b, c]);
                    
                    bestGuess           = gaussWithSlope(A, x0, g, b, c, neighbourFreqIndices');
                    
                case 'removeLorentzian'
                    % Fit a Lorentzian with slope
                    fittedAvgModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlope,...
                        'StartPoint', [A, x0, g, b, c]);
                    
                    bestGuess = lorentzianWithSlope(A, x0, g, b, c, neighbourFreqIndices');
                    
                case 'leaveSlopeMinusLorentzian'
                    % Fit a lortenzian with slope to the data.
                    fittedAvgModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlope,...
                        'StartPoint', [A, x0, g, b, c]);
                    
                    bestGuess = lorentzianWithSlope(A, x0, g, b, c, neighbourFreqIndices');
            end
            
% %             See how good the fit was.
%             slopeOnly = slope(b, c, neighbourFreqIndices');
%             hold on
%             plot(neighbourFreqIndices, neighbourData);
%             plot(neighbourFreqIndices,bestGuess);
%             plot(fittedAvgModel);
%             plot(neighbourFreqIndices, slopeOnly);
%             hold off
            
            % Redefine functions with constants for peak and width.
            gaussWithSlopeC = @(A, b, c, x)...
                (A.*exp(-2.*((x-fittedAvgModel.x0)./fittedAvgModel.g).^2)) + (b.*x) + c;

            gaussC = @(A, x)...
                (A.*exp(-2.*((x-fittedAvgModel.x0)./fittedAvgModel.g).^2));
            
            lorentzianWithSlopeC = @(A, b, c, x)...
                ((A .* (fittedAvgModel.g.^2))./((fittedAvgModel.g.^2) + ((x - fittedAvgModel.x0).^2))) + (b.*x) + c;

            lorentzianC = @(A, x)...
                ((A .* (fittedAvgModel.g.^2))./((fittedAvgModel.g.^2) + ((x - fittedAvgModel.x0).^2)));
            
            
            % Fit channels independently
            for chanIdx = 1:length(fftData(:,1))
                if strcmp(cfg.log,'yes')
                    % Get the data being fit.
                    neighbourData   = log(fftPow(chanIdx,neighbourFreqIndices));
                elseif strcmp(cfg.log,'no')
                    % Get the data being fit.
                    neighbourData   = fftPow(chanIdx,neighbourFreqIndices);
                end
                
                switch cfg.method
                    case 'leaveSlopeMinusGauss'
                        % Fit a gauss with slope to the data.
                        fittedModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlopeC,...
                                                    'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of the fitted Lorentz
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (cfg.strength * abs(fittedAvgModel.g)));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (cfg.strength * abs(fittedAvgModel.g)));
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx:tmpEndIdx);

                        fittedModelData         = gaussWithSlopeC(fittedModel.A,fittedModel.b,fittedModel.c,neighbourFreqIndices);
                        tmpStartVal             = fittedModelData(tmpStartIdx);
                        tmpEndVal               = fittedModelData(tmpEndIdx);
                        % Linear interpolation across the peak.
                        replacementData         = interp1([neighbourFreqIndices(tmpStartIdx),neighbourFreqIndices(tmpEndIdx)],[tmpStartVal,tmpEndVal],indicesToReplace,'linear');

                    case 'removeGauss'
                        % Fit a gauss with slope to the data.
                        fittedModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlopeC,...
                                                    'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of the fitted Lorentz
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (cfg.strength * abs(fittedAvgModel.g)));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (cfg.strength * abs(fittedAvgModel.g)));
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx:tmpEndIdx);

                        % Find the Gauss component.
                        onlyGauss               = gaussC(fittedModel.A, indicesToReplace);

                        % Remove the Gauss from the original data.
                        replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyGauss;

                    case 'removeLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlopeC,...
                            'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of cfg.strength gamma of fitted lorentzian.
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (cfg.strength * abs(fittedAvgModel.g)));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (cfg.strength * abs(fittedAvgModel.g)));

                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end

                        % Get the lorentzian component.
                        onlyLorentz             = lorentzianC(fittedModel.A, indicesToReplace);

                        % Remove it from the original data.
                        replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyLorentz;

                    case 'leaveSlopeMinusLorentzian'
                        % Fit a lortenzian with slope to the data.
                        fittedModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlopeC,...
                            'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of the fitted Lorentz
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (cfg.strength * abs(fittedAvgModel.g)));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (cfg.strength * abs(fittedAvgModel.g)));
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx:tmpEndIdx);

                        fittedModelData         = lorentzianWithSlopeC(fittedModel.A,fittedModel.b,fittedModel.c,neighbourFreqIndices);
                        tmpStartVal             = fittedModelData(tmpStartIdx);
                        tmpEndVal               = fittedModelData(tmpEndIdx);
                        % Linear interpolation across the peak.
                        replacementData         = interp1([neighbourFreqIndices(tmpStartIdx),neighbourFreqIndices(tmpEndIdx)],[tmpStartVal,tmpEndVal],indicesToReplace,'linear');

                        
                        
%                         % Replace under peak with just the slope
%                         replacementData         = slope(fittedModel.b,fittedModel.c, freqIndicesToReplace);
%                         
%                         
%                         % Get just the slope
%                         roughSlope                      = slope(fittedModel.b,fittedModel.c, freqToReplace);
%                         
%                         neighbourDataWithReplacement    = neighbourData;
%                         neighbourDataWithReplacement(tmpStartIdx:tmpEndIdx) = roughSlope;
%                         
%                         fittedSlope                     = fit(neighbourFreq',neighbourDataWithReplacement',slope,...
%                                                             'StartPoint',[fittedAvgModel.b, fittedAvgModel.c]);
%                         replacementData                 = fittedSlope(tmpStartIdx:tmpEndIdx);
%                         replacementData                 = replacementData';
%                         % The slope will always be slightly too low.
%                         hold on
%                         plot(neighbourFreqIndices,neighbourData)
%                         plot(freqIndicesToReplace,replacementData);
%                         plot(neighbourFreqIndices,fittedModelData);
%                         plot(fittedSlope,neighbourFreqIndices',tmp')
%                         plot(indicesToReplace,replacementData2);
                        
                end
%                 fittedData = lorentzianWithSlopeC(fittedModel.A,fittedModel.b,fittedModel.c, neighbourFreqIndices);
%                 
%                 % debug plots
%                 subplot(3,1,1);
%                 hold on
%                 plot(neighbourFreq, neighbourData)
%                 hold off
%                 subplot(3,1,2);
%                 hold on
%                 plot(freqToReplace,replacementData)
%                 plot(fittedModel)
%                 hold off
%                 subplot(3,1,3);
%                 hold on
%                 plot(neighbourFreq,fittedData)
%                 plot(neighbourFreq,bestGuess)
%                 hold off
%                 
                % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
                if strcmp(cfg.log,'yes')
                    fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), exp(replacementData));
                elseif strcmp(cfg.log,'no')
                    fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), replacementData);
                end
                
            end
        end
    case 'no'
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
            peakWidthData               = log(abs(avgData(widthFreqIndicesBound(1):widthFreqIndicesBound(end))));
            neighbourData               = log(abs(avgData(neighbourFreqIndices)));
        elseif strcmp(cfg.log,'no')
            % Get the chan avg data for those bounds.
            peakWidthData               = abs(avgData(widthFreqIndicesBound(1):widthFreqIndicesBound(end)));
            neighbourData               = abs(avgData(neighbourFreqIndices));
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
            g(peakIdx)              = sqrt(nearest(fftFreq,peakWidth(peakIdx)));
        end
        
        % Slope
        endVals             = mean(neighbourData((end - (round(g/4))):end));
        startVals           = mean(neighbourData(1:1+(round(g/4))));
        b                   = (endVals - startVals) / (neighbourFreqIndicesBound(end) - neighbourFreqIndicesBound(1));

        % Constant
        middleValue     = mean([endVals,startVals]);
        c               = middleValue - (b * median(neighbourFreqIndicesBound));
        
        % Define slope
        slope = @(b, c, x)...
                (b*x) + c;

        % Fit to the avg of channel data to get best guess.
        if length(peakFreq) == 2
            twoLorentzianWithSlope = @(A1, x01, g1, A2, x02, g2, b, c, x)...
                    ((A1 .* (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2))) + (b.*x) + c;

            twoLorentzian = @(A1, x01, g1, A2, x02, g2, x)...
                    ((A1 * (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2)));

            fittedAvgModel     = fit(neighbourFreqIndices', neighbourData', twoLorentzianWithSlope,...
                                'StartPoint', [A(1), x0(1), g(1), A(2), x0(2), g(2), b, c]);

            % Get the width of fitted lorentzians.
            x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02]);
            x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02]);
            gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2]));
            tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax)); % Magic number
            tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax)); % Magic number
            indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

            if length(indicesToReplace) > length(neighbourData)
                ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
            end

        elseif length(peakFreq) == 3
            threeLorentzianWithSlope = @(A1, x01, g1, A2, x02, g2, A3, x03, g3, b, c, x)...
                    ((A1 .* (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2))) + ((A3 * (g3.^2))./((g3.^2) + ((x - x03).^2))) + (b.*x) + c;

            threeLorentzian = @(A1, x01, g1, A2, x02, g2, A3, x03, g3, x)...
                    ((A1 .* (g1.^2))./((g1.^2) + ((x - x01).^2))) + ((A2 .* (g2.^2))./((g2.^2) + ((x - x02).^2))) + ((A3 * (g3.^2))./((g3.^2) + ((x - x03).^2)));

            fittedAvgModel     = fit(neighbourFreqIndices', neighbourData', threeLorentzianWithSlope,...
                                'StartPoint', [A(1), x0(1), g(1), A(2), x0(2), g(2), A(3), x0(3), g(3), b, c]);
            
            % Get the width of fitted lorentzians.
            x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
            x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
            gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2,fittedAvgModel.g3]));
            tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
            tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));
        end

        % Now apply per channel with some constants from the avg.  
        for chanIdx = 1:length(fftData(:,1))
            if strcmp(cfg.log,'yes')
                peakWidthData = log(abs(fftData(chanIdx,widthFreqIndicesBound(1):widthFreqIndicesBound(end))));
                neighbourData   = log(abs(fftData(chanIdx,neighbourFreqIndices)));
            elseif strcmp(cfg.log,'no')
                peakWidthData = abs(fftData(chanIdx,widthFreqIndicesBound(1):widthFreqIndicesBound(end)));
                neighbourData   = abs(fftData(chanIdx,neighbourFreqIndices));
            end
            

            if length(peakFreq) == 2
                 twoLorentzianWithSlopeC = @(A1, A2, b, c, x)...
                    ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2))) + (b*x) + c;

                twoLorentzianC = @(A1, A2, x)...
                    ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2)));

                switch cfg.method
                    case 'leaveSlopeMinusGauss'
                        
                    case 'removeGauss'
                        
                    case 'removeLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel     = fit(neighbourFreqIndices', neighbourData', twoLorentzianWithSlopeC,...
                                            'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of fitted lorentzians.
                        x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02]);
                        x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02]);
                        gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2]));
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));

                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end

                        % Get the lorentzian component.
                        onlyLorentz             = twoLorentzianC(fittedModel.A1, fittedModel.A2, indicesToReplace);

                        % Remove it from the original data.
                        replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyLorentz;

                        bestGuess       = twoLorentzianWithSlopeC(fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.b, fittedAvgModel.c,neighbourFreqIndices');

                    case 'leaveSlopeMinusLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel     = fit(neighbourFreqIndices', neighbourData', twoLorentzianWithSlopeC,...
                                            'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of fitted lorentzians.
                        x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02]);
                        x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02]);
                        gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2]));
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));

                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end
                        
                        % Get just the slope
                        replacementData         = slope(fittedModel.b,fittedModel.c, indicesToReplace);
                        
                        bestGuess       = twoLorentzianWithSlopeC(fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.b, fittedAvgModel.c,neighbourFreqIndices');
                end
                
            elseif length(peakFreq) == 3
                threeLorentzianWithSlopeC = @(A1, A2, A3, b, c, x)...
                        ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2))) + ((A3 * (fittedAvgModel.g3^2))./((fittedAvgModel.g3^2) + ((x - fittedAvgModel.x03).^2))) + (b*x) + c;
                
                threeLorentzianC = @(A1, A2, A3, x)...
                        ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2))) + ((A3 * (fittedAvgModel.g3^2))./((fittedAvgModel.g3^2) + ((x - fittedAvgModel.x03).^2)));

                
                
                switch cfg.method
                    case 'leaveSlopeMinusGauss'
                        
                    case 'removeGauss'
                        
                    case 'removeLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel     = fit(neighbourFreqIndices', neighbourData', threeLorentzianWithSlopeC,...
                                    'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.A3, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of fitted lorentzians.
                        x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                        x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                        gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2,fittedAvgModel.g3]));
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));
                        
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end

                        % Get the lorentzian component.
                        onlyLorentz             = threeLorentzianC(fittedModel.A1, fittedModel.A2, fittedAvgModel.A3, indicesToReplace);

                        % Remove it from the original data.
                        replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyLorentz;

                        bestGuess       = threeLorentzianWithSlopeC(fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.A3, fittedAvgModel.b, fittedAvgModel.c, neighbourFreqIndices');

                    case 'leaveSlopeMinusLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel     = fit(neighbourFreqIndices', neighbourData', threeLorentzianWithSlopeC,...
                                    'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.A3, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of fitted lorentzians.
                        x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                        x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                        gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2,fittedAvgModel.g3]));
                        tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (cfg.strength * gMax));
                        tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (cfg.strength * gMax));                        
                        
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);
                        
                        if length(indicesToReplace) > length(neighbourData)
                            ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                        end
                        
                        % Get just the slope
                        replacementData         = slope(fittedModel.b,fittedModel.c, indicesToReplace);
                        
                        bestGuess       = threeLorentzianWithSlopeC(fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.A3, fittedAvgModel.b, fittedAvgModel.c, neighbourFreqIndices');
                end
            end
%             % debug plots
%             subplot(3,1,1);
%             hold on
%             plot(neighbourFreqIndices, neighbourData)
%             plot(fittedModel)
%             hold off
%             subplot(3,1,2);
%             hold on
%             plot(indicesToReplace,replacementData)
%             plot(fittedModel)
%             hold off
%             subplot(3,1,3);
%             hold on
%             plot(fittedModel)
%             plot(bestGuess)
%             hold off

            % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
            if strcmp(cfg.log,'yes')
                fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), exp(replacementData));
            elseif strcmp(cfg.log,'no')
                fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), replacementData);
            end
        end
end

% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
% to ensure a real valued signal after iFFT
filteredTimeseries = ifft(fftData,[],2,'symmetric');

% Put it back into original structure.
filt            = data;
filt.trial{1}   = filteredTimeseries;

