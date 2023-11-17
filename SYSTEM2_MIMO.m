clear all, close all, clc

% DIGITAL COMMUNICATION SYSTEM M-QAM
%% BINARY SOURCE

% Read the image
img=imread("C:\Users\ASUS\OneDrive\Documentos\MATLAB\tigre.jpg");
% Convert the original image to grayscale
img_Gray=im2gray(img);

% Binarized Image
level=graythresh(img_Gray);% Calculate binarization threshold
imagen_binaria=im2bw(img_Gray,level); 
dataIn = imagen_binaria(:)';

% ------------PLOTS--------------
figure(1);
subplot(3,1,1);
imshow(img);
title('imagen original ');
subplot(3,1,2)
imshow(img_Gray) ,axis image, colormap gray, colorbar;
title('imagen en escala de grises');
subplot(3,1,3)
imshow(imagen_binaria), axis image, colormap gray;
title('imagen binaria');

% Convert the binary image into a bit vector

vectorbits = reshape(dataIn, 1, []);
num_total_bits= length(vectorbits); % Total number of bits
duracion_senal= 1; % 1 second
tasa_de_bits= num_total_bits/duracion_senal; % Bits transmission rate
disp(['Tasa de bits: ' num2str(tasa_de_bits)  'bps']);
ancho_banda= 2*tasa_de_bits;

% Calculate the Fourier Transform of the original signal (before modulation)
espectro_original = fft(vectorbits);

% Waveform parameters
Rs = 1; % Number of symbols transmitted per second (baud rate)
sps = 7; % Samples per symbol (oversampling factor)
fs = sps * Rs; % Sampling frequency
Ts = 1 / fs; % Sampling period

% Calculate the frequency axis
f = (-length(vectorbits)/2:length(vectorbits)/2-1) * fs / length(vectorbits);

% Calculate the magnitude of the spectrum
magnitud_espectro_original = abs(espectro_original);

% Visualize the spectrum
figure(2),
plot(f, fftshift(magnitud_espectro_original));
title('Espectro de la señal original (antes de la modulación)');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;

%% BASEBAND MODULATION M-QAM

modulation = 'QAM';  % Type of modulation to use.
M =64; 
num_bits_por_simbolo = log2(M); 
Mpam=sqrt(M);

Residuo2 = mod(length(vectorbits),log2(M));
if (Residuo2 ~= 0)  
        Relleno = zeros(1,log2(M)-Residuo2);
        vectorbitspar = [vectorbits,Relleno];
        num_simbolos = length(vectorbitspar)/num_bits_por_simbolo; %Número de símbolos.
        bits=vectorbitspar;
else 
        num_simbolos = length(vectorbits)/num_bits_por_simbolo; %Número de símbolos.
        bits=vectorbits;
end

simbolos_complejos = zeros(1, num_simbolos);
for i = 1:num_simbolos
    indice = bi2de(bits((i-1)*log2(M)+1:i*log2(M)), 'left-msb');
    simbolos_complejos(i) = qammod(indice, M);
end

% % Division into Real and Imaginary parts of Filtered Signal 
qam_Real = real(simbolos_complejos);
qam_Imag = imag(simbolos_complejos);
%
qam_Real = qam_Real.';
qam_Imag = qam_Imag.'; 


% Modulated signal spectrum
fs = 8; 
N_Qam = length(simbolos_complejos); 
frequencies_Qam = linspace(-fs / 2, fs / 2, N_Qam);% Frequency vector

% Calculate the FFT of the modulated signal
spectrum = fftshift(fft(simbolos_complejos, N_Qam));

% ------------PLOTS--------------

% Scatter plot of the transmitted symbol constellation
figure(3);

subplot(2,1,1);
scatter(real(simbolos_complejos), imag(simbolos_complejos), 'filled');
title(['Constelación de símbolos transmitidos (', num2str(M), '-QAM)']);
xlim([-10 10]), ylim([-10 10])
grid on;

subplot(2, 1, 2);
plot(frequencies_Qam    , abs(spectrum));
title(['Espectro de la señal (',num2str(M), '-QAM)']);
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
grid on;

%% SHAPING FILTER

d = 1; % Defines the separation between points
U = 7;
fs = U * Rs; % Sampling frequency
ts = 1 / fs;
fc = 2 * Rs; % Carrier frequency
c = 8; % Number of zeros taken by the filters p(t) and q(t)
span = 8; % Number of zeros taken by the filter

% Preparation of sequences for filtering (Signal oversampling)
Sreal = [real(simbolos_complejos) zeros(1, 2 * span)];
Simg = [imag(simbolos_complejos) zeros(1, 2 * span)];

% Transmitter filter design
rollOff=0.5; % Roll-off factor
F=rcosdesign(rollOff,span,sps,'sqrt'); 

% Symbol filterin 
X1=filter(F,1,upsample(Sreal,sps)); 
X2=filter(F,1,upsample(Simg,sps)); 
tx_filtrada = upsample(simbolos_complejos,span); 
conformed_pulses=filter(F,1,tx_filtrada);

% Modulated signal spectrum
N_sc = length(conformed_pulses);
frequencies_Nc = linspace(-fs / 2, fs / 2, N_sc); 

% Calculate the FFT of the shaped signal
spectrum2 = fftshift(fft(conformed_pulses, N_sc));

% ------------PLOTS--------------
figure(4);
subplot(2,1,1);
stem(F(1:50),'filled');
title('Pulso Conformador');

subplot(2, 1, 2);
stem(tx_filtrada(1:40),'filled');
title('Sobremuestro');

figure(5);
subplot(2, 1, 1);
stem(conformed_pulses)
xlim([10 50]) 
title('Señal Conformada');

subplot(2, 1, 2);
plot(frequencies_Nc, abs(spectrum2));
title('Espectro de la señal conformada');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
grid on;

fvtool(F, 'Analysis','Impulse');
%% MIMO TRANSMITTER

% Split symbols into two parts
num_simbolos_divididos = floor(num_simbolos / 2);
simbolos_complejos_antena1 = simbolos_complejos(1:num_simbolos_divididos);
simbolos_complejos_antena2 = simbolos_complejos(num_simbolos_divididos+1:end);

X1_antena1 = filter(F, 1, upsample(real(simbolos_complejos_antena1), sps));
X2_antena2 = filter(F, 1, upsample(real(simbolos_complejos_antena2), sps));

% Create a time axis for the conforming signal
t = (0:length(X1_antena1)-1) * Ts;

% Spectrum of modulated signals for each antenna
spectrum_antena1 = fftshift(fft(X1_antena1, N_sc));
spectrum_antena2 = fftshift(fft(X2_antena2, N_sc));

%--------PLOTS-------------

% Plots for each antenna
figure(6);
subplot(2, 1, 1);
plot(t, real(X1_antena1));
hold on;
plot(t, real(X2_antena2), 'r');
legend('Antena 1', 'Antena 2');
title('Señal Conformada para cada antena');

subplot(2, 1, 2);
plot(frequencies_Nc, abs(spectrum_antena1));
hold on;
plot(frequencies_Nc, abs(spectrum_antena2), 'r');
legend('Antena 1', 'Antena 2');
title('Espectro de la señal conformada para cada antena');
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
grid on;

% Create the column vector X
X = [X1_antena1; X2_antena2];

%% MIMO FADING CHANNEL

% Create fading matrix H (simulated)
H = (randn(2) + 1i * randn(2))/sqrt(2); % Square matrix of fading coefficients

% Transmit the signal through the MIMO channel
Y_antena1 = H(1,1) * X1_antena1 + H(1,2) * X2_antena2;
Y_antena2 = H(2,1) * X1_antena1 + H(2,2) * X2_antena2;


ebno=100000; 
Es = sum(abs(simbolos_complejos).^2) / length(simbolos_complejos);
sigma = sqrt(Es / (2 * log2(M) * ebno));
% Generate Gaussian noise
ruido_antena1 = sigma * (randn(size(Y_antena1)) + 1i * randn(size(Y_antena1)));
ruido_antena2 = sigma * (randn(size(Y_antena2)) + 1i * randn(size(Y_antena2)));

% Add noise to the received signals
Y_antena1_con_ruido = Y_antena1 + ruido_antena1;
Y_antena2_con_ruido = Y_antena2 + ruido_antena2;

%% MIMO RECEIVER

Frx=F;  % Filter for reception

Y11=filter(Frx,1,Y_antena1_con_ruido);
Y21=filter(Frx,1,Y_antena2_con_ruido);

% Downsample the filtered signal at the receiver
Y12=downsample(Y11,sps);
Y22=downsample(Y21,sps);

% Estimation of received symbols
Y= [Y12;Y22];

% Check the invertibility of matrix H
if det(H) == 0
    error('La matriz H no es invertible.');
end

% Calculate the inverse of matrix H
H_invertida = inv(H);

% Estimate transmitted symbols
X_estimado = H_invertida * Y;

% Concatenate estimated symbols to reconstruct the image
simbolos_estimados = [X_estimado(1, :), X_estimado(2, :)];

% Received symbols 
figure;
scatter(Y12,Y22,'filled'), title('Constelación de símbolos recibidos')
set(gcf,'color','w');
xlim([-20 20]), ylim([-20 20])
grid on

%% BASE-BAND DEMODULATION

demodulated_symbols = qamdemod(simbolos_estimados, M);
rx_bits = de2bi(demodulated_symbols,num_bits_por_simbolo)';
dataOut = reshape(rx_bits(1:length(dataIn)), size(dataIn));

% Calculate the Fourier transform of the signal

N_dem = length(simbolos_estimados); 
fs6 =1/Ts ; 
frequencies = linspace(-fs6 / 2, fs6 / 2, N_dem);

% Calculate the FFT of the demodulated signal
spectrum6 = fftshift(fft(simbolos_estimados, N_dem));

figure;

plot(frequencies, abs(spectrum6)); 
title('Espectro de la Señal demodulada');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
grid on;


%% RECOVER IMAGE

[row, col] = size(imagen_binaria); 
ReconstructedImage = reshape(dataOut, row, col);
ReconstructedImage = logical(ReconstructedImage); 

% ------------PLOTS--------------

figure;
imwrite(ReconstructedImage, 'imagen_reconstruida_mimo.jpg'); 
FINAL1= imread('imagen_reconstruida.jpg');
image(FINAL1);
axis image , colormap gray
title('Imagen Reconstruida MIMO');

% Calculate the Bit Error Rate (BER)
num_errors = sum(sum(dataIn ~= dataOut));
ber = num_errors / numel(dataIn);
disp(['Bit Error Rate (BER): ', num2str(ber)]);
