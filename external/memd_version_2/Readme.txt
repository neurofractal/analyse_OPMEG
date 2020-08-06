
MAIN Matlab CODE.
The Multivariate Empirical Mode Decomposition (MEMD) code - it is much computationally improved as compared with the old version, with up to an order of magnitude speedup.

Filename:  memd.m
	

SUPPORTING Matlab CODES:
If you would like to look into the filterbank property of MEMD or to plot the Hilbert-Huang spectrum, the following Matlab programs are available.


Filename:   g-noise-imfs.mat 
	Contains the variable 'imf_emd' which has IMFs obtained by applying EMD to separate 8 Gaussian noise channels
	of length N=5000;



Filename: filt_bank.m
	Plots the filter bank structure for multivariate IMFs
 	
	The input represents IMFs from 8 channel Gaussian noise of length L (from g-noise-imfs.mat)

 	The data set is then divided into 1000 element segments, which are then averaged together. 
 	The data set must therefore have be at least 1000 elements for this program to work. 

 	To save the computaional effort, we have already saved the imfs from
 	8-channel Gaussian noise obtained via MEMD (imf8_memd.mat) and via
 	standard EMD (imf_emd). Therefore, to plot the filter bank structure, you
 	can use the following code:

 	Example 1: (filter bank structure of MEMD) - load g-noise-imfs; filt_bank(imf8_memd);
 	Example 2: ( filter bank structure of standard EMD) - load g-noise-imfs; filt_bank(imf_emd);

Filename:    INST_FREQ_local.m, spectrogram_emd.m, disp_hhs.m:

	The IMF matrix (N X M : where N represents the number of IMF and M is the data length) can be fed into  `INST_FREQ_local.m’. 
	This function computes the Hilbert-Huang spectrum using the Hilbert transform, the instantaneous amplitude and the instantaneous 
	frequency for each IMF. The dimensions of the output matrices are same as that of the input IMF matrix. 

-	spectrogram_emd.m
	 transforms the 2D instantaneous amplitude and frequency into 3D spectrum matrix 
	 to plot the instantaneous amplitude contours on the time-frequency plane. 

-	disp_hhs.m then plots the 3D HHS spectrum of the spectrum matrix obtained as an output from spectrogram_emd.m.

	Here is an example:

	[instAmp,instFreq] = INST_FREQ_local(imf);           % calculating the instantaneous amplitude and the instantaneous frequency

	% Note that IMF matrix must be in the following format: (NxM), where N and M represent the number of 
	IMFs and the data length repectively. 

	spect = spectrogram_emd(instFreq,instAmp,1000);      % producing 3D Hilbert-Huang spectrum

	disp_hhs(spect);						     % Part of the free EMD toolbox provided by G. Rilling and P. Flandrin, and downloaded from
										 http://perso.ens-lyon.fr/patrick.flandrin/emd.html 
 
 
The following .mat files contain dataset which may serve as an input to
memd.m. Just load any of these files in matlab and send the resulting output 
vector as an input to the function nemd.

1) syn_12channel_inp.mat: contains synthetically generated 12 channel
data set (with combination of 5 tones (sinewaves) and noise added to few channels).

2) syn_16channel_inp.mat: contains synthetically generated 16 channel
data set (combined 6 tones (sinewaves) and noise added to some channels)

3) syn_hex_inp.mat: contains synthetically generated 6 channel data set
(combined 4 sinewaves and noise added to some channels - see the
Multivariate EMD paper and the Supplementary Material for more detaiil)

4) taichi_hex_inp.mat: contains hexavariate real world taichi dataset.
(two 3D recordings from intertial bodysensors combined into a single
hexavariate signal [left wrist and left ankle])