% Clear environment and close figures
clear; close all; clc;

%-------------------------------------#
%---------- Configuration ------------#
%-------------------------------------#
fs = 44100;                  % sampling rate
baud = 900;                  % symbol rate (QAM symbols per second)
Nbits = 4000;                % number of bits (ensure it's a multiple of 4)
f0 = 1800;                   % carrier Frequency

if Nbits < 0, Nbits = 0; end 

if Nbits > 0 && baud > 0 && fs > 0
    num_QAM_symbols = Nbits / 4; 
    if num_QAM_symbols == 0 % Handle Nbits < 4
        N = 0;
        t = [];
        Ns_for_repmat = 0;
    else
        Ns_per_QAM_symbol_ideal = fs / baud; 
        N_ideal = num_QAM_symbols * Ns_per_QAM_symbol_ideal;
        % Ensure N is an integer for total samples
        % Ensure Ns_for_repmat is an integer and N is a multiple of num_QAM_symbols
        if Ns_per_QAM_symbol_ideal >= 1
            Ns_for_repmat = floor(Ns_per_QAM_symbol_ideal);
            N = num_QAM_symbols * Ns_for_repmat;
            if N > 0
                t = (0:N-1)'/fs;      
            else
                t = []; % if Ns_for_repmat results in N=0
            end
        else 
            % Case: fs/baud < 1, meaning less than 1 sample per QAM symbol ideal
            % This is problematic. Set N to 0 or handle error.
            N = 0;
            t = [];
            Ns_for_repmat = 0;
            warning('fs/baud results in < 1 sample per QAM symbol. Simulation may not be meaningful.');
        end
    end
else
    N = 0;
    t = [];
    num_QAM_symbols = 0;
    Ns_for_repmat = 0;
end


symbolsToShow = 20; 
if baud > 0 && num_QAM_symbols > 0
    duration_symbols_to_show = symbolsToShow / baud;
    total_qam_duration = num_QAM_symbols / baud;
    timeDomainVisibleLimit = min(total_qam_duration, duration_symbols_to_show);
else
    timeDomainVisibleLimit = 0; 
end


pyC0_blue   = [0.1216, 0.4667, 0.7059]; 
pyC1_orange = [1.0000, 0.4980, 0.0549]; 
pyC2_green  = [0.1725, 0.6275, 0.1725]; 
pyC3_red    = [0.8392, 0.1529, 0.1569]; 
pyC4_purple = [0.5804, 0.4039, 0.7412]; 

%----------------------------------------#
%---------- QAM16 Modulation ------------#
%----------------------------------------#

if Nbits > 0
    inputBits = randn(Nbits,1) > 0; 
else
    inputBits = [];
end

inputSignal_fordisplay = []; 
if N > 0 && ~isempty(inputBits) && Ns_for_repmat > 0 && mod(Nbits, 4) == 0 
    Ns_per_bit_display = Ns_for_repmat / 4; 
    if Ns_per_bit_display >= 1
        inputSignal_fordisplay_temp = repmat(inputBits', floor(Ns_per_bit_display), 1);
        inputSignal_fordisplay = inputSignal_fordisplay_temp(:);
        if length(inputSignal_fordisplay) > N
            inputSignal_fordisplay = inputSignal_fordisplay(1:N);
        elseif length(inputSignal_fordisplay) < N && N > length(inputSignal_fordisplay)
            inputSignal_fordisplay = [inputSignal_fordisplay; zeros(N-length(inputSignal_fordisplay),1)];
        end
    else
        inputSignal_fordisplay = zeros(N,1); 
        % warning('Ns_per_bit_display < 1. "Digital Data Signal" plot might be misleading.');
    end
else
    if N > 0, inputSignal_fordisplay = zeros(N,1); end % fill with zeros if N exists but no bits to form signal
end


if(mod(Nbits, 4) == 0) && N > 0 && num_QAM_symbols > 0
    
    amplitude1_I_signal = 0.25;
    amplitude2_I_signal = 0.75;
    amplitude1_Q_signal = 0.25;
    amplitude2_Q_signal = 0.75;

    carrier1 = cos(2*pi*f0*t);
    carrier2 = cos(2*pi*f0*t + pi/2); 

    I1_bits = inputBits(1:4:end); 
    I2_bits = inputBits(2:4:end);
    Q1_bits = inputBits(3:4:end);
    Q2_bits = inputBits(4:4:end);
    
    I_symbols_levels = zeros(num_QAM_symbols,1);
    Q_symbols_levels = zeros(num_QAM_symbols,1);
    for x = 1:num_QAM_symbols
        I_symbols_levels(x) = TwoBitToAmplitudeMapper(I1_bits(x), I2_bits(x), amplitude1_I_signal, amplitude2_I_signal);
        Q_symbols_levels(x) = TwoBitToAmplitudeMapper(Q1_bits(x), Q2_bits(x), amplitude1_Q_signal, amplitude2_Q_signal);
    end
    
    I_signal_baseband_temp = repmat(I_symbols_levels', Ns_for_repmat, 1);
    I_signal_baseband = I_signal_baseband_temp(:);
    
    Q_signal_baseband_temp = repmat(Q_symbols_levels', Ns_for_repmat, 1);
    Q_signal_baseband = Q_signal_baseband_temp(:);
    
    I_signal_modulated = I_signal_baseband .* carrier1;
    Q_signal_modulated = Q_signal_baseband .* carrier2;
    QAM16_signal = I_signal_modulated + Q_signal_modulated;

    dataSymbols_int = zeros(num_QAM_symbols,1);
    for x = 1:num_QAM_symbols
        dataSymbols_int(x) = GetQam16Symbol(I1_bits(x), I2_bits(x), Q1_bits(x), Q2_bits(x));
    end

    noiseStandardDeviation = 0.055;
    noise1_src = noiseStandardDeviation * randn(num_QAM_symbols,1);
    noise2_src = noiseStandardDeviation * randn(num_QAM_symbols,1);
    noise3_src = noiseStandardDeviation * randn(num_QAM_symbols,1);
    noise4_src = noiseStandardDeviation * randn(num_QAM_symbols,1);

    Tx_symbols = complex(zeros(num_QAM_symbols,1));
    Rx_symbols = complex(zeros(num_QAM_symbols,1));

    for x = 1:num_QAM_symbols
        Tx_symbols(x) = Qam16SymbolMapper(dataSymbols_int(x), ...
                                          amplitude1_I_signal, amplitude1_Q_signal, ...
                                          amplitude2_I_signal, amplitude2_Q_signal);
                                          
        Rx_symbols(x) = Qam16SymbolMapper(dataSymbols_int(x), ...
                                          amplitude1_I_signal, amplitude1_Q_signal, ...
                                          amplitude2_I_signal, amplitude2_Q_signal, ...
                                          noise1_src(x), noise2_src(x), ...
                                          noise3_src(x), noise4_src(x));
    end
    
    %-------------------------------------------#
    %---------- Data Representation ------------#
    %-------------------------------------------#
    
    fig1 = figure('Name', 'Time Domain Signals');
    
    t_plot = []; % Initialize t_plot
    samples_to_show_in_plot = 0; % Initialize
    if N > 0 && timeDomainVisibleLimit > 0
        samples_to_show_in_plot_calc = floor(timeDomainVisibleLimit * fs);
        if samples_to_show_in_plot_calc == 0 && N > 0 
             samples_to_show_in_plot_calc = min(N, floor(symbolsToShow * Ns_for_repmat));
             if samples_to_show_in_plot_calc == 0 && N > 0 
                samples_to_show_in_plot_calc = min(N,1);
             end
        end
        if samples_to_show_in_plot_calc > N, samples_to_show_in_plot_calc = N; end
        samples_to_show_in_plot = samples_to_show_in_plot_calc; % Assign to final variable
        
        if samples_to_show_in_plot > 0 && ~isempty(t)
            t_plot = t(1:samples_to_show_in_plot);
        end
    end

    plot_info_signals = {
        {'Digital Data Signal (Source Code/ Block Diagram: "inputBits")', inputSignal_fordisplay, pyC1_orange};
        {'Digital I-Signal (Source Code/ Block Diagram: "I\_signal")', I_signal_baseband, pyC2_green};
        {'Modulated I-Signal (Source Code/ Block Diagram: "I\_signal\_modulated")', I_signal_modulated, pyC2_green};
        {'Digital Q-Signal (Source Code/ Block Diagram: "Q\_signal")', Q_signal_baseband, pyC3_red};
        {'Modulated Q-Signal (Source Code/ Block Diagram: "Q\_signal\_modulated")', Q_signal_modulated, pyC3_red};
        {'QAM16 Signal Modulated (Source Code/ Block Diagram: "QAM16\_signal")', QAM16_signal, pyC4_purple}
    };

    for i = 1:length(plot_info_signals)
        subplot(6,1,i);
        current_plot_data_vec = plot_info_signals{i}{2};
        if ~isempty(t_plot) && ~isempty(current_plot_data_vec) && length(current_plot_data_vec) >= samples_to_show_in_plot && samples_to_show_in_plot > 0
            plot(t_plot, current_plot_data_vec(1:samples_to_show_in_plot), 'Color', plot_info_signals{i}{3}); 
        else
            plot([],[]); 
        end
        title(plot_info_signals{i}{1});
        ylabel('Amplitude [V]');
        if ~isempty(t_plot) && t_plot(end) > 0 
            xlim([0, t_plot(end)]); 
        elseif N > 0 && timeDomainVisibleLimit > 0 
            xlim([0, timeDomainVisibleLimit]);
        else
            xlim([0,1]); 
        end
        grid on;
        if i == length(plot_info_signals), xlabel('Time [s]'); end 
    end

    allax = findall(gcf, 'Type', 'axes');
    for k=1:length(allax)
         allax(k).TitleHorizontalAlignment = 'center'; 
    end

    % Figure 2: Spectrum Plots 
    if N > 0 && ~isempty(QAM16_signal)
        fig2 = figure('Name', 'Spectrum Analysis');
        L_sig = length(QAM16_signal); 
        
        P1 = []; f_axis = []; % Initialize
        if L_sig > 0
            Y = fft(QAM16_signal);
            P2_calc = abs(Y/L_sig);
            P1 = P2_calc(1:floor(L_sig/2)+1);
            if length(P1) > 1 && L_sig > 1 
                 P1(2:end-1) = 2*P1(2:end-1); 
            end
            f_axis_calc = fs*(0:(floor(L_sig/2)))/L_sig;
            if isempty(f_axis_calc) && L_sig == 1 % Case L_sig=1 (DC only)
                f_axis = 0;
                if isempty(P1), P1 = P2_calc; end % P1 should be P2_calc directly
            elseif L_sig > 0
                f_axis = f_axis_calc;
            end
        end

        subplot(3,1,1); 
        if ~isempty(f_axis) && ~isempty(P1), plot(f_axis, P1, 'Color', pyC1_orange); else plot([],[]); end
        title('Magnitude Spectrum (Source Code/ Block Diagram: "QAM16\_signal")');
        ylabel('Magnitude (energy)'); 
        if fs > 0, xlim([0, min(6000, fs/2)]); else xlim([0, 6000]); end
        grid on;

        subplot(3,1,2); 
        noise_floor_val_log = 1e-9; % Adjusted for log spectrum "floor"
        if ~isempty(f_axis) && ~isempty(P1), plot(f_axis, 20*log10(P1 + noise_floor_val_log), 'Color', pyC1_orange); else plot([],[]); end
        title('Log. Magnitude Spectrum (Source Code/ Block Diagram: "QAM16\_signal")');
        ylabel('Magnitude (dB)'); 
        if fs > 0, xlim([0, min(6000, fs/2)]); else xlim([0, 6000]); end
        % ylim([-150, -30]); % Example: May need to adjust based on `noise_floor_val_log`
        grid on;

        subplot(3,1,3); 
        if L_sig >= 1 % pwelch or periodogram needs at least 1 sample
            % Use pwelch configured to be more like Python's psd(NFFT=len(signal))
            [Pxx, F_psd] = pwelch(QAM16_signal, rectwin(L_sig), 0, L_sig, fs);
            if ~isempty(F_psd) && ~isempty(Pxx), plot(F_psd, 10*log10(Pxx + eps), 'Color', pyC0_blue); else plot([],[]); end
        else
            plot([],[]);
        end
        title('Power Spectrum Density (PSD) (Source Code/ Block Diagram: "QAM16\_signal")');
        ylabel('Power/Frequency (dB/Hz)');
        xlabel('Frequency (Hz)'); 
        if fs > 0, xlim([0, min(6000, fs/2)]); else xlim([0, 6000]); end
        % ylim([-180, -10]); % Example: Adjust based on results
        grid on;
        
        allax2 = findall(gcf, 'Type', 'axes'); 
        for k=1:length(allax2)
             allax2(k).TitleHorizontalAlignment = 'center'; 
        end
    elseif N > 0 % N > 0 but QAM16_signal is empty (should not happen if N>0 and QAM syms > 0)
        fig2 = figure('Name', 'Spectrum Analysis (No QAM16_signal data)');
        for i=1:3, subplot(3,1,i); title('No QAM16\_signal data for spectrum'); plot([],[]); end
    end

    % Figure 3: Constellation Diagrams
    fig3 = figure('Name', 'Constellation Diagrams');
    
    commonLim = 1.5; 
    axTicks_ref = [-1, -0.5, 0, 0.5, 1];
    scale_val = 1; 
    
    if num_QAM_symbols > 0 && ~isempty(Tx_symbols) % Check num_QAM_symbols as well
        all_symbols_real_vec = [real(Tx_symbols(:)); real(Rx_symbols(:))];
        all_symbols_imag_vec = [imag(Tx_symbols(:)); imag(Rx_symbols(:))];
        
        if ~isempty(all_symbols_real_vec) && ~isempty(all_symbols_imag_vec)
            max_abs_val_symbols = max(abs([all_symbols_real_vec; all_symbols_imag_vec]));
            if isempty(max_abs_val_symbols) || isnan(max_abs_val_symbols) || max_abs_val_symbols == 0 
                max_abs_val_symbols = 1; 
            end
            commonLim_candidate_val = 1.1 * max_abs_val_symbols;
            if ~(isempty(commonLim_candidate_val) || commonLim_candidate_val == 0 || isnan(commonLim_candidate_val) || isinf(commonLim_candidate_val))
                commonLim = commonLim_candidate_val;
            end
        end
    end
    
    values_for_max_tick_vec = abs(axTicks_ref(axTicks_ref ~= 0));
    if isempty(values_for_max_tick_vec)
        max_abs_tick_val_calc = 0; 
    else
        max_abs_tick_val_calc = max(values_for_max_tick_vec); 
    end
    
    denominator_base_val_calc = max(max_abs_tick_val_calc, 1); 
    denominator_val_final_calc = 1.1 * denominator_base_val_calc; 
    
    if (denominator_val_final_calc == 0) || isnan(denominator_val_final_calc) || isinf(denominator_val_final_calc)
        scale_val_candidate_calc = 1; 
    else
        scale_val_candidate_calc = commonLim / denominator_val_final_calc;
    end
    
    if ~(isnan(scale_val_candidate_calc) || isinf(scale_val_candidate_calc))
        scale_val = scale_val_candidate_calc;
    else
        scale_val = 1; 
    end

    axTicks = axTicks_ref * scale_val;

    subplot(2,2,1); 
    if num_QAM_symbols > 0 && ~isempty(Tx_symbols), plot(real(Tx_symbols), imag(Tx_symbols),'.', 'Color',pyC1_orange, 'MarkerSize',10); else plot([],[]); end
    title('Tx (Source Code/ Block Diagram: "Tx\_symbols")');
    xlabel('Inphase [V]'); ylabel('Quadrature [V]');
    xlim([-commonLim,commonLim]); ylim([-commonLim,commonLim]);
    xticks(axTicks); yticks(axTicks);
    grid on; axis square;

    subplot(2,2,2); 
    if num_QAM_symbols > 0 && ~isempty(Tx_symbols)
        plot(real(Tx_symbols), imag(Tx_symbols),'-', 'Color', pyC0_blue); hold on;
        plot(real(Tx_symbols), imag(Tx_symbols),'.', 'Color', pyC1_orange, 'MarkerSize',10); hold off;
    else plot([],[]); end
    title('Tx with Trajectory (Source Code/ Block Diagram: "Tx\_symbols")');
    xlabel('Inphase [V]'); ylabel('Quadrature [V]');
    xlim([-commonLim,commonLim]); ylim([-commonLim,commonLim]);
    xticks(axTicks); yticks(axTicks);
    grid on; axis square;

    subplot(2,2,3); 
    if num_QAM_symbols > 0 && ~isempty(Rx_symbols), plot(real(Rx_symbols), imag(Rx_symbols),'.', 'Color',pyC1_orange, 'MarkerSize',10); else plot([],[]); end
    title('Rx (Source Code/ Block Diagram: "Rx\_symbols")');
    xlabel('Inphase [V]'); ylabel('Quadrature [V]');
    xlim([-commonLim,commonLim]); ylim([-commonLim,commonLim]);
    xticks(axTicks); yticks(axTicks);
    grid on; axis square;

    subplot(2,2,4); 
    if num_QAM_symbols > 0 && ~isempty(Rx_symbols)
        plot(real(Rx_symbols), imag(Rx_symbols),'-', 'Color', pyC0_blue); hold on;
        plot(real(Rx_symbols), imag(Rx_symbols),'.', 'Color', pyC1_orange, 'MarkerSize',10); hold off;
    else plot([],[]); end
    title('Rx with Trajectory (Source Code/ Block Diagram: "Rx\_symbols")');
    xlabel('Inphase [V]'); ylabel('Quadrature [V]');
    xlim([-commonLim,commonLim]); ylim([-commonLim,commonLim]);
    xticks(axTicks); yticks(axTicks);
    grid on; axis square;
    
    allax3 = findall(gcf, 'Type', 'axes'); 
    for k=1:length(allax3)
         allax3(k).TitleHorizontalAlignment = 'center'; 
    end

else % Conditions for main simulation not met
    if N == 0 && Nbits > 0 && (baud <= 0 || fs <=0) 
        disp("Error! Nbits > 0 but baud or fs is not positive, resulting in 0 total samples (N). Adjust parameters.");
    elseif N == 0 && Nbits > 0 && fs/baud < 1/num_QAM_symbols && num_QAM_symbols > 0 % Ns_for_repmat becomes 0
        disp("Error! Nbits > 0 but fs/baud is too small, resulting in Ns_for_repmat = 0 and N = 0. Adjust parameters.");
    elseif Nbits == 0
        disp("Nbits is 0. No modulation to perform or plot.");
    elseif mod(Nbits,4) ~=0 && Nbits > 0
        disp("Error! Number of bits has to be a multiple of 4. Number of Bits entered: "+ num2str(Nbits)+".");
    else
        disp("An unspecified configuration issue prevented the simulation. Check Nbits, baud, fs.");
    end
    
    % Create empty figures to prevent script from appearing to do nothing if it errors out early
    % (and if no other figures were created)
    if isempty(findall(0,'type','figure')) % Check if any figures were already created
        figure('Name', 'Time Domain Signals (No Data)'); subplot(1,1,1); title('No valid data for time domain plots');
        figure('Name', 'Spectrum Analysis (No Data)'); subplot(1,1,1); title('No valid data for spectrum plots');
        figure('Name', 'Constellation Diagrams (No Data)'); subplot(1,1,1); title('No valid data for constellation plots');
    end
end

% ------------------------------------- %
% ---------- LOCAL FUNCTIONS ---------- %
% (Copy and paste the three local functions here as before)
% ------------------------------------- %

% Maps a two bit input to a certain amplitude provided
function amplitude = TwoBitToAmplitudeMapper(bit1, bit2, amplitude1, amplitude2)
    if (~bit1 && ~bit2)
        amplitude = -amplitude1;
    elseif (~bit1 && bit2)
        amplitude = -amplitude2;
    elseif (bit1 && ~bit2)
        amplitude = amplitude1;
    elseif (bit1 && bit2)
        amplitude = amplitude2;
    else
        amplitude = 0; % Should not happen with boolean inputs
    end
end

% Used for symbol creation. Returns a decimal number from a 4 bit input
function symbol = GetQam16Symbol(bit1, bit2, bit3, bit4)
    if(~bit1 && ~bit2 && ~bit3 && ~bit4) %0000
        symbol = 0;
    elseif(~bit1 && ~bit2 && ~bit3 && bit4) %0001
        symbol = 1;
    elseif(~bit1 && ~bit2 && bit3 && ~bit4) %0010
        symbol = 2;
    elseif(~bit1 && ~bit2 && bit3 && bit4) %0011
        symbol = 3;
    elseif(~bit1 && bit2 && ~bit3 && ~bit4) %0100
        symbol = 4;
    elseif(~bit1 && bit2 && ~bit3 && bit4) %0101
        symbol = 5;
    elseif(~bit1 && bit2 && bit3 && ~bit4) %0110
        symbol = 6;
    elseif(~bit1 && bit2 && bit3 && bit4) %0111
        symbol = 7;
    elseif(bit1 && ~bit2 && ~bit3 && ~bit4) %1000
        symbol = 8;
    elseif(bit1 && ~bit2 && ~bit3 && bit4) %1001
        symbol = 9;
    elseif(bit1 && ~bit2 && bit3 && ~bit4) %1010
        symbol = 10;
    elseif(bit1 && ~bit2 && bit3 && bit4) %1011
        symbol = 11;
    elseif(bit1 && bit2 && ~bit3 && ~bit4) %1100
        symbol = 12;
    elseif(bit1 && bit2 && ~bit3 && bit4) %1101
        symbol = 13;
    elseif(bit1 && bit2 && bit3 && ~bit4) %1110
        symbol = 14;
    elseif(bit1 && bit2 && bit3 && bit4) %1111
        symbol = 15;
    else
        symbol = -1; % Error case
    end
end

% Maps a given symbol to a complex signal. Optionally, noise and phase offset can be added.
function complex_signal = Qam16SymbolMapper(symbol_val, amplitude_I1, amplitude_Q1, amplitude_I2, amplitude_Q2, noise1, noise2, noise3, noise4, phaseOffset1, phaseOffset2)
    if nargin < 6, noise1 = 0; end
    if nargin < 7, noise2 = 0; end
    if nargin < 8, noise3 = 0; end
    if nargin < 9, noise4 = 0; end
    if nargin < 10, phaseOffset1 = 0; end
    if nargin < 11, phaseOffset2 = 0; end

    complex_signal = 0 + 0i; % Initialize

    if(symbol_val == 0) %0000
        complex_signal = sqrt(amplitude_I1^2 + amplitude_Q1^2)*(cos(deg2rad(225) + phaseOffset1)+ 1i *sin(deg2rad(225) + phaseOffset2)) + (noise1 + 1i*noise2);
    elseif(symbol_val == 1) %0001
        complex_signal = sqrt(amplitude_I2^2 + amplitude_Q1^2 )*(cos(deg2rad(198) + phaseOffset1)+ 1i *sin(deg2rad(198) + phaseOffset2)) + (noise3 + 1i*noise2);
    elseif(symbol_val == 2) %0010
        complex_signal = sqrt(amplitude_I1^2 + amplitude_Q2^2)*(cos(deg2rad(251) + phaseOffset1)+ 1i *sin(deg2rad(251) + phaseOffset2)) + (noise1 + 1i*noise4);
    elseif(symbol_val == 3) %0011
        complex_signal = sqrt(amplitude_I2^2 + amplitude_Q2^2)*(cos(deg2rad(225) + phaseOffset1)+ 1i *sin(deg2rad(225) + phaseOffset2)) + (noise3 + 1i*noise4);
    elseif(symbol_val == 4) %0100
        complex_signal = sqrt(amplitude_I1^2 + amplitude_Q1^2)*(cos(deg2rad(315) + phaseOffset1)+ 1i *sin(deg2rad(315) + phaseOffset2)) + (noise1 + 1i*noise2);
    elseif(symbol_val == 5) %0101
        complex_signal = sqrt(amplitude_I1^2 + amplitude_Q2^2)*(cos(deg2rad(288) + phaseOffset1)+ 1i *sin(deg2rad(288) + phaseOffset2)) + (noise1 + 1i*noise4);
    elseif(symbol_val == 6) %0110
        complex_signal = sqrt(amplitude_I2^2 + amplitude_Q1^2)*(cos(deg2rad(342) + phaseOffset1)+ 1i *sin(deg2rad(342) + phaseOffset2)) + (noise2 + 1i*noise3);
    elseif(symbol_val == 7) %0111
        complex_signal = sqrt(amplitude_I2^2 + amplitude_Q2^2)*(cos(deg2rad(315) + phaseOffset1)+ 1i *sin(deg2rad(315) + phaseOffset2)) + (noise3 + 1i*noise4);
    elseif(symbol_val == 8) %1000
        complex_signal = sqrt(amplitude_I1^2 + amplitude_Q1^2)*(cos(deg2rad(135) + phaseOffset1)+ 1i *sin(deg2rad(135) + phaseOffset2)) + (noise1 + 1i*noise2);
    elseif(symbol_val == 9) %1001
        complex_signal = sqrt(amplitude_I1^2 + amplitude_Q2^2)*(cos(deg2rad(108) + phaseOffset1)+ 1i *sin(deg2rad(108) + phaseOffset2)) + (noise1 + 1i*noise4);
    elseif(symbol_val == 10) %1010
        complex_signal = sqrt(amplitude_I2^2 + amplitude_Q1^2)*(cos(deg2rad(162) + phaseOffset1)+ 1i *sin(deg2rad(162) + phaseOffset2)) + (noise3 + 1i*noise2);
    elseif(symbol_val == 11) %1011
        complex_signal = sqrt(amplitude_I2^2 + amplitude_Q2^2)*(cos(deg2rad(135) + phaseOffset1)+ 1i *sin(deg2rad(135) + phaseOffset2)) + (noise3 + 1i*noise4);
    elseif(symbol_val == 12) %1100
        complex_signal =  sqrt(amplitude_I1^2 + amplitude_Q1^2)*(cos(deg2rad(45) + phaseOffset1)+ 1i *sin(deg2rad(45) + phaseOffset2)) + (noise1 + 1i*noise2);
    elseif(symbol_val == 13) %1101
        complex_signal = sqrt(amplitude_I2^2 + amplitude_Q1^2)*(cos(deg2rad(18) + phaseOffset1)+ 1i *sin(deg2rad(15) + phaseOffset2)) + (noise2 + 1i*noise3);
    elseif(symbol_val == 14) %1110
        complex_signal = sqrt(amplitude_I1^2 + amplitude_Q2^2)*(cos(deg2rad(71) + phaseOffset1)+ 1i *sin(deg2rad(75) + phaseOffset2)) + (noise1 + 1i*noise4);
    elseif(symbol_val == 15) %1111
        complex_signal = sqrt(amplitude_I2^2 + amplitude_Q2^2)*(cos(deg2rad(45) + phaseOffset1)+ 1i *sin(deg2rad(45) + phaseOffset2)) + (noise3 + 1i*noise4);
    else
        complex_signal = complex(0,0); 
    end
end