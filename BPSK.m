% MATLAB implementation of BPSK modulation, constellation diagrams,
% time-domain plots, and frequency-domain spectra.
% Translated from Python code.

clear; close all; clc;

%% Configuration
fs = 44100;             % sampling rate [Hz]
baud = 900;             % symbol rate [symbols/s]
Nbits = 4000;           % number of bits
f0 = 1800;              % carrier Frequency [Hz]
Ns = floor(fs/baud);    % number of Samples per Symbol (Python int() truncates)
N = Nbits * Ns;         % Total Number of Samples
t = (0:N-1)/fs;         % time vector [s] (row vector)

% Limit for representation of time domain signals for better visibility.
symbolsToShow = 20;
timeDomainVisibleLimit = min(Nbits/baud, symbolsToShow/baud);

% Limit for representation of frequency domain signals for better visibility.
% (Note: Python f = r_[0:N/2.0]/N*fs implies N/2 points for spectrum)
% MATLAB f_axis for plotting FFT typically uses N/2 or N/2+1 points.
% Frequency axis for plotting FFT results (single-sided spectrum)
f_axis_full = linspace(0, fs/2, floor(N/2)+1); % For N_fft/2+1 points

sideLobesToShow = 9;
sideLobeWidthSpectrum = baud;
lowerLimit = max(0, f0 - sideLobeWidthSpectrum*(1 + sideLobesToShow));
upperLimit = f0 + sideLobeWidthSpectrum*(1 + sideLobesToShow);

carrier1 = cos(2*pi*f0*t); % carrier1 will be a row vector

%% Helper functions (equivalent to Python's GetBpskSymbol and BpskSymbolMapper)
% For BPSK:
% Bit 0 (False) -> Symbol 0 -> Maps to +Amplitude (cos(0))
% Bit 1 (True)  -> Symbol 1 -> Maps to -Amplitude (cos(pi))

% No explicit GetBpskSymbol needed if we map bits directly.
% BpskSymbolMapper logic will be incorporated directly.

%% BPSK Modulation
% Modulator Input
inputBits = randn(Nbits,1) > 0; % Logical column vector (true/false)

% Digital-to-Analog Conversion (Baseband Signal)
% Python: inputSignal = (np.tile(inputBits*2-1,(1,Ns))).ravel()
% inputBits*2-1: True -> 1, False -> -1
symbols_pm1 = 2*double(inputBits) - 1; % Column vector of +1 and -1
inputSignal = repelem(symbols_pm1, Ns); % Column vector, (Nbits*Ns) x 1

% Multiplicator / mixer
% Ensure correct dimensions for element-wise multiplication.
% inputSignal is column, carrier1 is row. MATLAB's .* handles this by broadcasting.
BPSK_signal = inputSignal .* carrier1'; % Make carrier1 column or inputSignal row.
                                       % Or rely on broadcasting like: inputSignal .* carrier1
                                       % Let's make carrier explicit column for clarity here.
BPSK_signal = inputSignal .* cos(2*pi*f0*t'); % t' makes time column, so carrier is column

%% Prepare BPSK Constellation Diagram
amplitude = 1;

% Generate noise.
noiseStandardDeviation = 0.12;
% Python noise1/noise2 are for I and Q components for Rx_symbols
noiseI = noiseStandardDeviation * randn(Nbits, 1);
noiseQ = noiseStandardDeviation * randn(Nbits, 1);

% Transmitted symbols (ideal constellation points)
% Based on Python logic: True bit (1) -> -amplitude; False bit (0) -> +amplitude
Tx_symbols = zeros(Nbits, 1);
Tx_symbols(inputBits == true) = -amplitude; % Bit 1 (True) maps to -amp
Tx_symbols(inputBits == false) = amplitude; % Bit 0 (False) maps to +amp

% Received symbols (add noise to Tx_symbols)
% Rx_symbols = mapped_symbol + (noiseI + 1j*noiseQ)
Rx_symbols = Tx_symbols + (noiseI + 1i*noiseQ);

%% Plot of BPSK Constellation Diagram
figure('Name', 'Constellation Diagram BPSK');
% sgtitle('Constellation Diagram BPSK'); % For newer MATLAB versions

subplot(2,2,1);
plot(real(Tx_symbols), imag(Tx_symbols), '.');
title('Tx Symbols');
xlabel('Inphase [V]');
ylabel('Quadrature [V]');
xlim([-2 2]); xticks([-1.5 0 1.5]);
ylim([-2 2]); yticks([-1.5 0 1.5]); % Python sharey='row' implies consistent Y
grid on;

subplot(2,2,2);
plot(real(Tx_symbols), imag(Tx_symbols),'-o'); % Plot with trajectory
title('Tx with Trajectory');
xlabel('Inphase [V]');
% ylabel('Quadrature [V]'); % Not needed if sharey
xlim([-2 2]); xticks([-1.5 0 1.5]);
ylim([-2 2]); yticks([-1.5 0 1.5]);
grid on;

subplot(2,2,3);
plot(real(Rx_symbols), imag(Rx_symbols), '.');
title('Rx Symbols');
xlabel('Inphase [V]');
ylabel('Quadrature [V]');
xlim([-2 2]); xticks([-1.5 0 1.5]);
ylim([-2 2]); yticks([-1.5 0 1.5]);
grid on;

subplot(2,2,4);
plot(real(Rx_symbols), imag(Rx_symbols),'-o'); % Plot with trajectory
title('Rx with Trajectory');
xlabel('Inphase [V]');
% ylabel('Quadrature [V]'); % Not needed if sharey
xlim([-2 2]); xticks([-1.5 0 1.5]);
ylim([-2 2]); yticks([-1.5 0 1.5]);
grid on;

%% Plot of BPSK Time-Domain Signals
figure('Name', 'BPSK Modulation');
% sgtitle('BPSK Modulation');

subplot(3,1,1);
plot(t, inputSignal); % t is row, inputSignal is column. MATLAB plot handles this.
title('Baseband Signal (Upsampled Symbols: +/-1)');
xlabel('Time [s]');
xlim([0, timeDomainVisibleLimit]);
ylabel('Amplitude [V]');
grid on; % Python 'dotted' linestyle, MATLAB default is solid for grid

subplot(3,1,2);
plot(t, carrier1); % t is row, carrier1 (original def) is row.
title('Carrier Signal');
xlabel('Time [s]');
xlim([0, timeDomainVisibleLimit]);
ylabel('Amplitude [V]');
grid on;

subplot(3,1,3);
plot(t, BPSK_signal); % t is row, BPSK_signal is column.
title('BPSK Modulated Signal');
xlabel('Time [s]');
xlim([0, timeDomainVisibleLimit]);
ylabel('Amplitude [V]');
grid on;

% Adjust spacing if needed, e.g., by manually setting subplot positions
% or using tightlayout in newer MATLAB if sgtitle causes overlap.
% For older MATLAB, can use:
% suptitle_text = 'BPSK Modulation';
% h_suptitle = suptitle(suptitle_text); % May require downloading suptitle for older versions
% set(h_suptitle, 'FontSize', 12);

%% Plot of Modulated Signal and Spectrum
figure('Name', 'BPSK Spectra');
% sgtitle('BPSK Spectra');

% --- Magnitude Spectrum ---
N_sig = length(BPSK_signal); % Should be N
X_fft = fft(BPSK_signal);
P2_abs = abs(X_fft/N_sig);
% Single-sided spectrum for real signals
if mod(N_sig,2) == 0 % N_sig is Even
    P1_mag = P2_abs(1:N_sig/2+1);
    P1_mag(2:end-1) = 2*P1_mag(2:end-1);
else % N_sig is Odd
    P1_mag = P2_abs(1:(N_sig+1)/2);
    P1_mag(2:end) = 2*P1_mag(2:end);
end
f_plot_axis = linspace(0, fs/2, length(P1_mag));

subplot(4,1,1);
plot(f_plot_axis, P1_mag);
title('Magnitude Spectrum');
xlabel('Frequency [Hz]');
ylabel('|X(f)|');
xlim([lowerLimit, upperLimit]);
grid on;

% --- Log. Magnitude Spectrum (dB) ---
subplot(4,1,2);
plot(f_plot_axis, 20*log10(P1_mag));
title('Log. Magnitude Spectrum (dB)');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xlim([lowerLimit, upperLimit]);
grid on;

% --- Power Spectrum Density (PSD) ---
subplot(4,1,3);
% Python ax.psd(BPSK_signal,NFFT=len(t),Fs=fs)
% NFFT=len(t) which is N_sig. Use a rectangular window and no overlap
% to be closer to a direct FFT-based PSD if that was implied by NFFT=len(t).
% Or use MATLAB's pwelch defaults for a Welch estimate.
% pwelch(BPSK_signal, rectwin(N_sig), 0, N_sig, fs);
pwelch(BPSK_signal, [], [], [], fs); % Using pwelch defaults for window, overlap, NFFT
title('Power Spectral Density (PSD)');
xlabel('Frequency [Hz]');
ylabel('Power/Frequency [dB/Hz]'); % pwelch default is dB/Hz
xlim([lowerLimit, upperLimit]);
% grid on; % pwelch typically adds its own grid

% --- Carrier Spectrum ---
N_carrier = length(carrier1); % Should be N
C_fft = fft(carrier1);
P2_carrier_abs = abs(C_fft/N_carrier);
if mod(N_carrier,2) == 0 % N_carrier is Even
    P1_carrier_mag = P2_carrier_abs(1:N_carrier/2+1);
    P1_carrier_mag(2:end-1) = 2*P1_carrier_mag(2:end-1);
else % N_carrier is Odd
    P1_carrier_mag = P2_carrier_abs(1:(N_carrier+1)/2);
    P1_carrier_mag(2:end) = 2*P1_carrier_mag(2:end);
end
f_plot_carrier_axis = linspace(0, fs/2, length(P1_carrier_mag));

subplot(4,1,4);
plot(f_plot_carrier_axis, P1_carrier_mag);
title('Carrier Spectrum');
xlabel('Frequency [Hz]');
ylabel('|C(f)|');
xlim([lowerLimit, upperLimit]);
grid on;

% plt.subplots_adjust(hspace=0.5) % In MATLAB, adjust subplot positions manually if needed
% For example:
% all_axes = findall(gcf, 'type', 'axes');
% for i = 1:length(all_axes)
%     pos = get(all_axes(i), 'Position');
%     set(all_axes(i), 'Position', [pos(1) pos(2) pos(3) pos(4)*0.85]); % Reduce height slightly
% end

disp('MATLAB script finished.');