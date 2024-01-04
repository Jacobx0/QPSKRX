clear
clc
close all
%% Parameters
M = 4; % For QPSK
symbolRate = 1e6; % Symbol rate: 1 Msps
fs = 2e6; % Sampling rate: 2 MSPS
numSymbols = 10000; % Number of QPSK symbols to generate
rolloffFactor = 0.5; % Roll-off factor for the RRC filter
filterSpan = 6; % Span of the RRC filter in symbol durations
sps = fs / symbolRate; % Samples per symbol
satelliteVelocity = 7500;
lightSpeed = 3e8;
carrierFreq = 12e9;
delayVal = 1.3; % Samples of timing error
phaseOffset = 25;

%% Transmitter
% Generate random symbols
dataSymbols = randi([0 M-1], numSymbols, 1);

% QPSK Modulation
qpskModulator = comm.QPSKModulator('BitInput', false);
modulatedSignal = qpskModulator(dataSymbols);

% Root Raised Cosine Filter
rrcFilterTx = comm.RaisedCosineTransmitFilter(...
    'RolloffFactor', rolloffFactor, ...
    'FilterSpanInSymbols', filterSpan, ...
    'OutputSamplesPerSymbol', sps);

% Apply RRC filter (pulse shaping)
txSignal = rrcFilterTx(modulatedSignal);
% Adding time delay
varDelay = dsp.VariableFractionalDelay;
delaySig = varDelay(txSignal,delayVal);


%% Channel
% Doppler Shift Calculation
dopplerShift = (satelliteVelocity / lightSpeed) * carrierFreq;

% Apply Doppler effect (Frequency shift)
t = (0:length(delaySig)-1)'/fs;
rxSignal = delaySig .* exp(1j * 2 * pi * dopplerShift * t);

freqOffset = comm.PhaseFrequencyOffset(PhaseOffset=phaseOffset, ...
                                       FrequencyOffset=dopplerShift, ...
                                       SampleRate=fs);
rxSignal = freqOffset(delaySig);
% Parameters for noise addition
SNR_dB = 20; 
% Add AWGN to the received signal
rxSignalWithNoise = awgn(rxSignal, SNR_dB, 'measured');

%% RECEIVER
% Raise the received signal to the fourth power
signalToFourth = rxSignalWithNoise.^4;

% Perform FFT on the fourth-power signal
N = length(signalToFourth);
fftSignal = fft(signalToFourth, N);
freqAxis = linspace(-fs/2, fs/2, N); % Frequency axis

% Shift the FFT and find the peak within the baud rate range
fftShifted = fftshift(fftSignal);

[~, peakIndex] = max(abs(fftShifted));
% Calculate the original frequency offset (divide by 4)
% Adjusting the peakIndex to align with the shifted frequency axis

estimatedOffset = fs/4 + freqAxis(peakIndex) / 4;
% Correct the frequency offset in the original signal
correctedSignal = rxSignalWithNoise .* exp(-1j * 2 * pi * estimatedOffset * (0:N-1).' / fs);

% Create a root-raised cosine filter
rrcFilterRx = rcosdesign(rolloffFactor, filterSpan, fs/symbolRate);

% Apply the filter to the corrected signal
mfSignal = conv(correctedSignal, rrcFilterRx, 'same');



symbolSyncE = comm.SymbolSynchronizer("TimingErrorDetector","Early-Late (non-data-aided)","SamplesPerSymbol",sps);
[rxSyncE,errE] = symbolSyncE(mfSignal);
% Fine Frequency/phase estimation/correction
fine = comm.CarrierSynchronizer( ...
    'DampingFactor',0.4, ...
    'NormalizedLoopBandwidth',0.001, ...
    'SamplesPerSymbol',1, ...
    'Modulation','QPSK');

fineSig = fine(rxSyncE);


constDiagram = comm.ConstellationDiagram( ...
    'ReferenceConstellation',qammod(0:M-1,M), ...
    'ChannelNames',{'Before convergence','After convergence'}, ...
    'ShowLegend',true);

constDiagram([fineSig(1:1000) fineSig(9001:end)])

scope = spectrumAnalyzer(SampleRate=fs);
scope(txSignal)
