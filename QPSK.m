% Main script
fs = 44100;                  % Tần số lấy mẫu (sampling rate)
baud = 900;                  % Tốc độ ký hiệu (symbol rate)
Nbits = 4000;                % Số bit
f0 = 1800;                   % Tần số sóng mang (carrier frequency)
Ns = floor(fs / baud);       % Số mẫu trên mỗi ký hiệu
N = Nbits * Ns;              % Tổng số mẫu
t = (0:N-1) / fs;            % Vector thời gian
f = (0:floor(N/2)) * fs / N; % Vector tần số cho phổ

symbolsToShow = 20;
timeDomainVisibleLimit = min(Nbits / baud, symbolsToShow / baud);

% Kiểm tra số bit phải là bội của 2 (vì QPSK dùng 2 bit trên mỗi ký hiệu)
if mod(Nbits, 2) == 0
    % Tạo dãy bit ngẫu nhiên
    inputBits = rand(Nbits,1) > 0.5;  % Mảng logic (0 hoặc 1)
    inputSignal = kron(2*double(inputBits)' - 1, ones(1,Ns));  % Tạo tín hiệu đầu vào (±1)

    % Tạo sóng mang
    carrier1 = cos(2*pi*f0*t);        % Sóng mang I (cos)
    carrier2 = cos(2*pi*f0*t + pi/2); % Sóng mang Q (sin)

    % Tách bit thành kênh I và Q (serial-to-parallel)
    I_bits = inputBits(1:2:end);  % Bit lẻ cho kênh I
    Q_bits = inputBits(2:2:end);  % Bit chẵn cho kênh Q

    % Tạo tín hiệu I và Q (mỗi bit lặp lại 2*Ns mẫu)
    I_signal = kron(2*double(I_bits)' - 1, ones(1,2*Ns));
    Q_signal = kron(2*double(Q_bits)' - 1, ones(1,2*Ns));

    % Điều chế tín hiệu
    I_signal_modulated = I_signal .* carrier1;
    Q_signal_modulated = Q_signal .* carrier2;
    QPSK_signal = I_signal_modulated + Q_signal_modulated;

    % Chuẩn bị cho biểu đồ chòm sao
    dataSymbols = zeros(length(I_bits),1);
    for x = 1:length(I_bits)
        dataSymbols(x) = GetQpskSymbol(I_bits(x), Q_bits(x));
    end

    amplitude_I_signal = 1;
    amplitude_Q_signal = 1;
    noiseStandardDeviation = 0.065;
    noise1 = noiseStandardDeviation * randn(size(dataSymbols));
    noise2 = noiseStandardDeviation * randn(size(dataSymbols));

    Tx_symbols = zeros(size(dataSymbols), 'like', 1i);
    Rx_symbols = zeros(size(dataSymbols), 'like', 1i);
    for x = 1:length(dataSymbols)
        Tx_symbols(x) = QpskSymbolMapper(dataSymbols(x), amplitude_I_signal, amplitude_Q_signal);
        Rx_symbols(x) = QpskSymbolMapper(dataSymbols(x), amplitude_I_signal, amplitude_Q_signal, noise1(x), noise2(x));
    end

    % Vẽ biểu đồ chòm sao
    figure;
    subplot(2,2,1);
    plot(real(Tx_symbols), imag(Tx_symbols), '.');
    title('Tx');
    xlabel('Inphase [V]');
    ylabel('Quadrature [V]');
    xlim([-2,2]);
    ylim([-2,2]);

    subplot(2,2,2);
    plot(real(Tx_symbols), imag(Tx_symbols), '-', real(Tx_symbols), imag(Tx_symbols), '.');
    title('Tx with Trajectory');
    xlabel('Inphase [V]');
    ylabel('Quadrature [V]');
    xlim([-2,2]);
    ylim([-2,2]);

    subplot(2,2,3);
    plot(real(Rx_symbols), imag(Rx_symbols), '.');
    title('Rx');
    xlabel('Inphase [V]');
    ylabel('Quadrature [V]');
    xlim([-2,2]);
    ylim([-2,2]);

    subplot(2,2,4);
    plot(real(Rx_symbols), imag(Rx_symbols), '-', real(Rx_symbols), imag(Rx_symbols), '.');
    title('Rx with Trajectory');
    xlabel('Inphase [V]');
    ylabel('Quadrature [V]');
    xlim([-2,2]);
    ylim([-2,2]);

    % Vẽ tín hiệu trong miền thời gian
    figure;
    subplot(6,1,1);
    plot(t, inputSignal);
    title('Digital Data Signal');
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    xlim([0, timeDomainVisibleLimit]);
    grid on;

    subplot(6,1,2);
    plot(t, I_signal);
    title('Digital I-Signal');
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    xlim([0, timeDomainVisibleLimit]);
    grid on;

    subplot(6,1,3);
    plot(t, I_signal_modulated);
    title('Modulated I-Signal');
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    xlim([0, timeDomainVisibleLimit]);
    grid on;

    subplot(6,1,4);
    plot(t, Q_signal);
    title('Digital Q-Signal');
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    xlim([0, timeDomainVisibleLimit]);
    grid on;

    subplot(6,1,5);
    plot(t, Q_signal_modulated);
    title('Modulated Q-Signal');
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    xlim([0, timeDomainVisibleLimit]);
    grid on;

    subplot(6,1,6);
    plot(t, QPSK_signal);
    title('QPSK Signal Modulated');
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    xlim([0, timeDomainVisibleLimit]);
    grid on;

    % Vẽ tín hiệu trong miền tần số
    figure;
    subplot(3,1,1);
    spectrum = abs(fft(QPSK_signal));
    plot(f, spectrum(1:length(f)));
    title('Magnitude Spectrum');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude');
    xlim([0,6000]);
    grid on;

    subplot(3,1,2);
    log_spectrum = 20*log10(spectrum);
    plot(f, log_spectrum(1:length(f)));
    title('Log. Magnitude Spectrum');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [dB]');
    xlim([0,6000]);
    grid on;

    subplot(3,1,3);
    [pxx, f_pxx] = pwelch(QPSK_signal, [], [], [], fs);
    plot(f_pxx, 10*log10(pxx));
    title('Power Spectral Density (PSD)');
    xlabel('Frequency [Hz]');
    ylabel('Power/Frequency [dB/Hz]');
    xlim([0,6000]);
    grid on;
else
    disp('Error! Number of bits has to be a multiple of 2.');
end

% Hàm ánh xạ bit thành ký hiệu QPSK (0, 1, 2, 3)
function symbol = GetQpskSymbol(bit1, bit2)
    if ~bit1 && ~bit2
        symbol = 0;
    elseif ~bit1 && bit2
        symbol = 1;
    elseif bit1 && ~bit2
        symbol = 2;
    elseif bit1 && bit2
        symbol = 3;
    else
        symbol = -1;
    end
end

% Hàm ánh xạ ký hiệu QPSK thành số phức
function signal = QpskSymbolMapper(symbol, amplitude_I, amplitude_Q, noise1, noise2, phaseOffset1, phaseOffset2)
    % Xử lý các tham số tùy chọn
    if nargin < 4, noise1 = 0; end
    if nargin < 5, noise2 = 0; end
    if nargin < 6, phaseOffset1 = 0; end
    if nargin < 7, phaseOffset2 = 0; end

    % Tính biên độ tổng
    amplitude = sqrt(amplitude_I^2 + amplitude_Q^2);

    % Ánh xạ ký hiệu thành góc pha
    if symbol == 0
        angle = 45 * pi / 180;
    elseif symbol == 1
        angle = 135 * pi / 180;
    elseif symbol == 2
        angle = 225 * pi / 180;
    elseif symbol == 3
        angle = 315 * pi / 180;
    else
        signal = 0 + 0i;
        return;
    end

    % Tạo tín hiệu số phức
    signal = amplitude * (cos(angle + phaseOffset1) + 1i * sin(angle + phaseOffset2)) + (noise1 + 1i*noise2);
end