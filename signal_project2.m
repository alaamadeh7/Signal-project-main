%% SIGNALS PROCESSING AND ANALYSIS
%% --------------------------------

%% 1. SIGNAL GENERATION (Ramp, Unit Step)
n = -10:10; % Discrete time range

% Ramp Signal
ramp_signal = max(0, n);

% Unit Step Signal
unit_step_signal = double(n >= 0);

% Plot signals
figure;
subplot(2,1,1);
stem(n, ramp_signal, 'LineWidth', 1.5, 'Color', 'b');
grid on;
title('Ramp Signal');
xlabel('n'); ylabel('Amplitude');

subplot(2,1,2);
stem(n, unit_step_signal, 'LineWidth', 1.5, 'Color', 'r');
grid on;
title('Unit Step Signal');
xlabel('n'); ylabel('Amplitude');

%% 2. OPERATIONS ON SIGNALS
% Addition
addition = ramp_signal + unit_step_signal;

% Subtraction
subtraction = ramp_signal - unit_step_signal;

% Multiplication
multiplication = ramp_signal .* unit_step_signal;

% Division (with small epsilon to avoid division by zero)
division = zeros(size(ramp_signal));
division(unit_step_signal ~= 0) = ramp_signal(unit_step_signal ~= 0) ./ unit_step_signal(unit_step_signal ~= 0);

% Plot results
figure;
subplot(4,1,1);
stem(n, addition, 'LineWidth', 1.5, 'Color', 'g');
grid on;
title('Addition of Signals');

subplot(4,1,2);
stem(n, subtraction, 'LineWidth', 1.5, 'Color', 'm');
grid on;
title('Subtraction of Signals');

subplot(4,1,3);
stem(n, multiplication, 'LineWidth', 1.5, 'Color', 'k');
grid on;
title('Multiplication of Signals');

subplot(4,1,4);
stem(n, division, 'LineWidth', 1.5, 'Color', 'c');
grid on;
title('Division of Signals');

%% 3. SAMPLING OF SIGNALS
t = 0:0.01:2*pi; % Continuous time
x_cont = sin(t);

% Sampled Signal
T_s = 0.2; % Sampling interval
t_sampled = 0:T_s:max(t);
x_sampled = sin(t_sampled);

% Plot continuous and sampled signals
figure;
plot(t, x_cont, 'LineWidth', 1.5, 'Color', 'b'); hold on;
stem(t_sampled, x_sampled, 'r', 'LineWidth', 1.5);
grid on;
title('Sampling of a Signal');
legend('Continuous Signal', 'Sampled Signal');
xlabel('Time'); ylabel('Amplitude');

%% 4. CONVOLUTION OF SIGNALS - Manual and MATLAB Convolution
% Define range for the convolution result
n_conv = 0:length(ramp_signal) + length(unit_step_signal) - 2;

% Manual Convolution
conv_manual = zeros(1, length(n_conv));
for n_idx = 1:length(n_conv)
    sum_result = 0;
    for k = 1:length(ramp_signal)
        if (n_idx - k + 1 > 0) && (n_idx - k + 1 <= length(unit_step_signal))
            sum_result = sum_result + ramp_signal(k) * unit_step_signal(n_idx - k + 1);
        end
    end
    conv_manual(n_idx) = sum_result;
end

% MATLAB Convolution
conv_builtin = conv(ramp_signal, unit_step_signal);

% Plot comparison
figure;
subplot(2,1,1);
stem(n_conv, conv_manual, 'LineWidth', 1.5, 'Color', 'b');
grid on;
title('Manual Convolution');
xlabel('n');
ylabel('Amplitude');

subplot(2,1,2);
stem(n_conv, conv_builtin, 'LineWidth', 1.5, 'Color', 'r');
grid on;
title('MATLAB Convolution');
xlabel('n');
ylabel('Amplitude');

%% 5. FILTERING
% Noisy signal
noise = rand(1, length(ramp_signal)); % Random noise
noisy_signal = ramp_signal + noise;  % Add noise to ramp signal

% Averaging Filter
M = 3; % Filter size
h = ones(1, M) / M; % Filter coefficients
filtered_signal = filter(h, 1, noisy_signal); % Apply the filter

% Plot noisy and filtered signals
figure;
plot(n, noisy_signal, 'LineWidth', 1.5, 'Color', 'r'); hold on;
plot(n, filtered_signal, 'LineWidth', 1.5, 'Color', 'b');
grid on;
title('Filtering: Noisy vs. Filtered Signal');
legend('Noisy Signal', 'Filtered Signal');
xlabel('n (Sample Index)');
ylabel('Amplitude');


%% 6. FREQUENCY RESPONSE OF FILTER
[H, w] = freqz(h, 1, 512);
figure;
plot(w/pi, abs(H), 'LineWidth', 1.5, 'Color', 'b');
grid on;
title('Frequency Response of the Averaging Filter');
xlabel('Normalized Frequency (\times \pi rad/sample)');
ylabel('Magnitude');

%% 7. DISCRETE-TIME FOURIER TRANSFORM (DTFT)
omega = linspace(-pi, pi, 512);
X = fftshift(fft(ramp_signal, 512));

figure;
plot(omega, abs(X), 'LineWidth', 1.5, 'Color', 'b');
grid on;
title('Discrete-Time Fourier Transform (DTFT)');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude');

%% 8. DISCRETE FOURIER TRANSFORM (DFT)
N = 8; % Number of points
dft_signal = fft(ramp_signal, N);

figure;
stem(0:N-1, abs(dft_signal), 'LineWidth', 1.5, 'Color', 'g');
grid on;
title('Discrete Fourier Transform (DFT)');
xlabel('Frequency Index');
ylabel('Magnitude');

%% 9. Z-TRANSFORM
syms z n_sym;
n_sym = -10:10; % Discrete range
z_transform = sum(ramp_signal .* z.^(-n_sym));

disp('Z-Transform:');
disp(z_transform);