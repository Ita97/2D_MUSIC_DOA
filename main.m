clear;close;
load rx_signal.mat x fc SensorArray 

% load rx_signal:
% x - rx waveform [SxT] with S number of samples and R number of receiving antennas, 
% fc - carrier frequency, 
% Sensor Array - a phased.NRRectangularPanelArray of dimension 12x12 

sps = true; % spacial smoothing
fb = true; % forward and backward smoothing
plot = true; % plot spectrum

phased.NRRectangularPanelArray
K = 1; % num signals

az_range = -60:0.5:60;
el_range = -45:0.5:45;


%% Custom MUSIC
[est_aoa, P] = MUSIC_DOA_2D(x, SensorArray, fc, az_range, el_range, K,sps,fb,plot);
%est_aoa=wrapTo180(est_aoa+reachable_gNBs(idx).orientation')
real_aoa = [55.8299,   -7.3617]'
est_aoa
