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
%   cfg.independentPeaks    = 'yes' or 'no'. Whether to model peaks independently or together.

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
        if strcmp(cfg.independentPeaks,'yes')
            [peakVal,peakFreq,halfHeightWidth]  = findpeaks(powMedLog,freq,'MinPeakProminence',minPeakProminence,'WidthReference','halfheight','MinPeakDistance',cfg.minPeakDistance);
        elseif strcmp(cfg.independentPeaks,'no')
            [peakVal,peakFreq,halfHeightWidth]  = findpeaks(powMedLog,freq,'MinPeakProminence',minPeakProminence,'WidthReference','halfheight');
        end
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

%% Remove Gaussian component from amplitude in spectrum
% Run the fft
fftData = fft(timeSeries,nsamples,2); % calculate fft to obtain spectrum that will be interpolated
avgData = median(abs(fftData),1);


switch cfg.independentPeaks
    case 'yes'
        for peakIdx = 1:length(peakFftFreq)
            % Find constants for all channels from the channel avg first
            
            % Find the indices for the frequencies of the peak.
            threeStdDevFreqIndices  = nearest(fftFreq,threeStdDevFftFreqRange(1,peakIdx)):nearest(fftFreq,threeStdDevFftFreqRange(2,peakIdx));
            neighbourFreqIndices    = nearest(fftFreq,neighbourFftFreqRange(1,peakIdx)):nearest(fftFreq,neighbourFftFreqRange(2,peakIdx));

            % Guesses for fit (all channels)
            % Samples Width
            widthSamples = length(threeStdDevFreqIndices) / 6; % Magic number...

            
            % Define function for Gaussian/Lorentzian with and without
            % slopes.
            gaussWithSlope = @(A, x0, s, b, c, x)...
                (A*exp(-2*((x-x0)/s).^2)) + (b*x) + c;

            gauss = @(A, x0, s, x)...
                (A*exp(-2*((x-x0)/s).^2));

            slope = @(b, c, x)...
                (b*x) + c;

            lorentzianWithSlope = @(A, x0, g, b, c, x)...
                ((A * (g^2))./((g^2) + ((x - x0).^2))) + (b*x) + c;

            lorentzian = @(A, x0, g, x)...
                ((A * (g^2))./((g.^2) + ((x - x0).^2)));
            
            
            % Get the data being fit
            neighbourData   = log(abs(avgData(neighbourFreqIndices)));
            threeStdDevData = log(abs(avgData(threeStdDevFreqIndices)));

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
                    fittedAvgModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlope,...
                                                'StartPoint', [(peak - fitFloor), centre, widthSamples, slopeC, constant]);

                case 'removeGauss'
                    % Fit a gauss with slope to the data.
                    fittedAvgModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlope,...
                                                'StartPoint', [(peak - fitFloor), centre, widthSamples, slopeC, constant]);

                case 'removeLorentzian'
                    % Fit a Lorentzian with slope
                    fittedAvgModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlope,...
                        'StartPoint', [(peak - fitFloor), centre, widthSamples, slopeC, constant]);

                case 'leaveSlopeMinusLorentzian'
                    % Fit a lortenzian with slope to the data.
                    fittedAvgModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlope,...
                        'StartPoint', [(peak - fitFloor), centre, widthSamples, slopeC, constant]);
            end
            
            % Redefine functions with constants for peak and width.
            gaussWithSlopeC = @(A, b, c, x)...
                (A*exp(-2*((x-fittedAvgModel.x0)/fittedAvgModel.s).^2)) + (b*x) + c;

            gaussC = @(A, x)...
                (A*exp(-2*((x-fittedAvgModel.x0)/fittedAvgModel.s).^2));

            lorentzianWithSlopeC = @(A, b, c, x)...
                ((A * (fittedAvgModel.g^2))./((fittedAvgModel.g^2) + ((x - fittedAvgModel.x0).^2))) + (b*x) + c;

            lorentzianC = @(A, x)...
                ((A * (fittedAvgModel.g^2))./((fittedAvgModel.g.^2) + ((x - fittedAvgModel.x0).^2)));
            
            % Fit channels independently (may need to revise this)
            for chanIdx = 1:length(fftData(:,1))
                % Get the data being fit.
                neighbourData   = log(abs(fftData(chanIdx,neighbourFreqIndices)));
                threeStdDevData = log(abs(fftData(chanIdx,threeStdDevFreqIndices)));

                switch cfg.method
                    case 'leaveSlopeMinusGauss'
                        % Fit a gauss with slope to the data.
                        fittedModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlopeC,...
                                                    'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of 3 std.dev of the fitted Gauss
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (3 * fittedAvgModel.s));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (3 * fittedAvgModel.s));
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                        % Get just the slope
                        replacementData         = slope(fittedModel.b,fittedModel.c, indicesToReplace);

                    case 'removeGauss'
                        % Fit a gauss with slope to the data.
                        fittedModel    = fit(neighbourFreqIndices', neighbourData', gaussWithSlopeC,...
                                                    'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of 3 std.dev of the fitted Gauss
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (3 * fittedAvgModel.s));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (3 * fittedAvgModel.s));
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                        % Find the Gauss component.
                        onlyGauss               = gaussC(fittedModel.A, indicesToReplace);

                        % Remove the Gauss from the original data.
                        replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyGauss;

                    case 'removeLorentzian'
                        % Fit a Lorentzian with slope
                        fittedModel = fit(neighbourFreqIndices', neighbourData', lorentzianWithSlopeC,...
                            'StartPoint', [fittedAvgModel.A, fittedAvgModel.b, fittedAvgModel.c]);

                        % Get the width of 6 gamma of fitted lorentzian.
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (6 * abs(fittedAvgModel.g)));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (6 * abs(fittedAvgModel.g)));

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

                        % Get the width of 3 std.dev of the fitted Gauss
                        tmpStartIdx             = nearest(neighbourFreqIndices,fittedAvgModel.x0 - (6 * abs(fittedAvgModel.g)));
                        tmpEndIdx               = nearest(neighbourFreqIndices,fittedAvgModel.x0 + (6 * abs(fittedAvgModel.g)));
                        indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                        % Get just the slope
                        replacementData         = slope(fittedModel.b,fittedModel.c, indicesToReplace);
                end
                
                % debug plots
                subplot(2,1,1);
                hold on
                plot(neighbourFreqIndices, neighbourData)
                plot(fittedModel)
                hold off
                subplot(2,1,2);
                hold on
                plot(indicesToReplace,replacementData)
                plot(fittedModel)
                hold off
                % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
                fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), exp(replacementData));
            end
        end
    case 'no'
        % Find the total width indices for all peaks together
        maxThreeStdDevFreqIndices   = nearest(fftFreq,min(threeStdDevFftFreqRange,[],'all')):nearest(fftFreq,max(threeStdDevFftFreqRange,[],'all'));
        neighbourFreqIndices        = nearest(fftFreq,min(neighbourFftFreqRange,[],'all')):nearest(fftFreq,max(neighbourFftFreqRange,[],'all'));
        
        % Find constants from the avg data.
        maxThreeStdDevData      = log(abs(avgData(maxThreeStdDevFreqIndices)));
        neighbourData           = log(abs(avgData(neighbourFreqIndices)));

        % Guesses for fit (channel level)
        % Slope
        slopeC       = (neighbourData(end) - neighbourData(1)) / (maxThreeStdDevFreqIndices(end) - maxThreeStdDevFreqIndices(1));

        % Constant
        constant    = neighbourData(1) - (slopeC *  neighbourFreqIndices(1));

        A           = zeros(size(peakFreq));
        x0          = zeros(size(peakFreq));
        g           = zeros(size(peakFreq));
        b           = slopeC;
        c           = constant;
        for peakIdx = 1:length(peakFreq)
            % Find indices for the peak.
            threeStdDevFreqIndices  = nearest(fftFreq,threeStdDevFftFreqRange(1,peakIdx)):nearest(fftFreq,threeStdDevFftFreqRange(2,peakIdx));
%             threeStdDevData         = log(abs(avgData(threeStdDevFreqIndices)));
            threeStdDevData         = abs(avgData(threeStdDevFreqIndices));
            
            % Guesses for fit (peak level)
            % Floor
            fitFloor        = mean(threeStdDevData,2);

            % Peak
            peak            = threeStdDevData(round(length(threeStdDevData)/2));

            % Centre freq
            centre          = threeStdDevFreqIndices(round(length(threeStdDevFreqIndices)/2));

            % Guesses for fit (all channels)
            % Gamma
            widthSamples    = length(threeStdDevFreqIndices) / 6;

            % Guesses for fit (channel level)
            A(peakIdx)      = peak - fitFloor;
            x0(peakIdx)     = centre;
            g(peakIdx)      = widthSamples;
        end

        if length(peakFreq) == 2
            twoLorentzianWithSlope = @(A1, x01, g1, A2, x02, g2, b, c, x)...
                    ((A1 * (g1^2))./((g1^2) + ((x - x01).^2))) + ((A2 * (g2^2))./((g2^2) + ((x - x02).^2))) + (b*x) + c;

            twoLorentzian = @(A1, x01, g1, A2, x02, g2, x)...
                    ((A1 * (g1^2))./((g1^2) + ((x - x01).^2))) + ((A2 * (g2^2))./((g2^2) + ((x - x02).^2)));

            % Start values for gamma have been set to 1. The guess is not
            % well worked at the moment. 
            fittedAvgModel     = fit(neighbourFreqIndices', neighbourData', twoLorentzianWithSlope,...
                                'StartPoint', [A(1), x0(1), 1, A(2), x0(2), 1, b, c]);

            % Get the width of fitted lorentzians.
            x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02]);
            x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02]);
            gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2]));
            tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax));
            tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax));
            indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

            if length(indicesToReplace) > length(neighbourData)
                ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
            end

        elseif length(peakFreq) == 3
            threeLorentzianWithSlope = @(A1, x01, g1, A2, x02, g2, A3, x03, g3, b, c, x)...
                    ((A1 * (g1^2))./((g1^2) + ((x - x01).^2))) + ((A2 * (g2^2))./((g2^2) + ((x - x02).^2))) + ((A3 * (g3^2))./((g3^2) + ((x - x03).^2))) + (b*x) + c;

            threeLorentzian = @(A1, x01, g1, A2, x02, g2, A3, x03, g3, x)...
                    ((A1 * (g1^2))./((g1^2) + ((x - x01).^2))) + ((A2 * (g2^2))./((g2^2) + ((x - x02).^2))) + ((A3 * (g3^2))./((g3^2) + ((x - x03).^2)));

            fittedAvgModel     = fit(neighbourFreqIndices', neighbourData', threeLorentzianWithSlope,...
                                'StartPoint', [A(1), x0(1), 1, A(2), x0(2), 1, A(3), x0(3), 1, b, c]);
            
            % Get the width of fitted lorentzians.
            x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
            x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
            gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2,fittedAvgModel.g3]));
            tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax));
            tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax));
        end

        
        % Now apply per channel with some constants from the avg.
        for chanIdx = 1:length(fftData(:,1))
            maxThreeStdDevData = log(abs(fftData(chanIdx,maxThreeStdDevFreqIndices)));
            neighbourData   = log(abs(fftData(chanIdx,neighbourFreqIndices)));

            if length(peakFreq) == 2
                twoLorentzianWithSlopeC = @(A1, A2, b, c, x)...
                        ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2))) + (b*x) + c;

                twoLorentzianC = @(A1, A2, x)...
                        ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2)));

                fittedModel     = fit(neighbourFreqIndices', neighbourData', twoLorentzianWithSlopeC,...
                                    'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.b, fittedAvgModel.c]);
                                
                % Get the width of fitted lorentzians.
                x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02]);
                x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02]);
                gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2]));
                tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax));
                tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax));
                
                indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                if length(indicesToReplace) > length(neighbourData)
                    ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                end

                % Get the lorentzian component.
                onlyLorentz             = twoLorentzianC(fittedModel.A1, fittedModel.A2, indicesToReplace);

                % Remove it from the original data.
                replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyLorentz;

            elseif length(peakFreq) == 3
                threeLorentzianWithSlopeC = @(A1, A2, A3, b, c, x)...
                        ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2))) + ((A3 * (fittedAvgModel.g3^2))./((fittedAvgModel.g3^2) + ((x - fittedAvgModel.x03).^2))) + (b*x) + c;
                
                threeLorentzianC = @(A1, A2, A3, x)...
                        ((A1 * (fittedAvgModel.g1^2))./((fittedAvgModel.g1^2) + ((x - fittedAvgModel.x01).^2))) + ((A2 * (fittedAvgModel.g2^2))./((fittedAvgModel.g2^2) + ((x - fittedAvgModel.x02).^2))) + ((A3 * (fittedAvgModel.g3^2))./((fittedAvgModel.g3^2) + ((x - fittedAvgModel.x03).^2)));

                fittedModel     = fit(neighbourFreqIndices', neighbourData', threeLorentzianWithSlopeC,...
                                    'StartPoint', [fittedAvgModel.A1, fittedAvgModel.A2, fittedAvgModel.A3, fittedAvgModel.b, fittedAvgModel.c]);
                
                % Get the width of fitted lorentzians.
                x0Min                   = min([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                x0Max                   = max([fittedAvgModel.x01,fittedAvgModel.x02,fittedAvgModel.x03]);
                gMax                    = max(abs([fittedAvgModel.g1,fittedAvgModel.g2,fittedAvgModel.g3]));
                tmpStartIdx             = nearest(neighbourFreqIndices,x0Min - (6 * gMax));
                tmpEndIdx               = nearest(neighbourFreqIndices,x0Max + (6 * gMax));
                
                indicesToReplace        = neighbourFreqIndices(tmpStartIdx):neighbourFreqIndices(tmpEndIdx);

                if length(indicesToReplace) > length(neighbourData)
                    ft_error('Frequencies being replaced are wider than the specified width. Increase neighbour width');
                end

                % Get the lorentzian component.
                onlyLorentz             = threeLorentzianC(fittedModel.A1, fittedModel.A2, fittedModel.A3, indicesToReplace);

                % Remove it from the original data.
                replacementData         = neighbourData(tmpStartIdx:tmpEndIdx) - onlyLorentz;
            end
            % debug plots
            subplot(2,1,1);
            hold on
            plot(neighbourFreqIndices, neighbourData)
            plot(fittedModel)
            hold off
            subplot(2,1,2);
            hold on
            plot(indicesToReplace,replacementData)
            plot(fittedModel)
            hold off
            % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
            fftData(chanIdx,indicesToReplace) = bsxfun(@times, exp(bsxfun(@times,angle(fftData(chanIdx,indicesToReplace)),1i)), exp(replacementData));
        end
end

% complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
% to ensure a real valued signal after iFFT
filteredTimeseries = ifft(fftData,[],2,'symmetric');

% Put it back into original structure.
filt            = data;
filt.trial{1}   = filteredTimeseries;

