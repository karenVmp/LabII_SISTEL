clear all, close all, clc

% DIGITAL COMMUNICATION SYSTEM M-QAM
%% BINARY SOURCE

% Read the image
img=imread("ubicación: \tigre.jpg");
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

%% PASS-BAND MODULATION

% Generation of modulated signals in phase and quadrature
t=0:Ts:(num_simbolos+2*span)*(1/Rs)-Ts; % Time vector
fc=2*Rs; 

X11= sqrt(2)*X1.*cos(2*pi*t*fc); % In-phase modulated signal
X21=-sqrt(2)*X2.*sin(2*pi*t*fc); % Quadrature modulated signal
X=X11+X21;% Signal to transmit

% Calculate the Fourier transform of the modulated signal
N = length(X); 
fs3 =5 ; 
frequencies = linspace(-fs3 / 2, fs3 / 2, N);

% Calculate the FFT of the modulated signal
spectrum3 = fftshift(fft(X, N));

% ------------PLOTS-------------

figure;
subplot(2, 1, 1);
plot(t,X);
title('Señal Modulada Pasa Banda');
xlabel('tiempo (t)')
ylabel('Amplitud')
xlim([0 50]) 

subplot(2, 1, 2);
plot(frequencies, abs(spectrum3)); % Escala en dB
title('Espectro de la Señal Modulada PB');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
grid on;

%% AWGN CHANNEL

% symbol energy
ebno=100000; 
num0=zeros(1,2*c);

Es = sum(abs(simbolos_complejos).^2) / length(simbolos_complejos);

sigma=sqrt(Es/(2*log2(M)*ebno));% Determine the noise variance
Z=sigma*randn(1,length(X));
Y=X+Z;

%% MULTIPATH CHANNEL
k=0.4;%Define a multipath coefficient k representing attenuation due to multipropagation.

t1=80;% Define a delay value t1 indicating the delay of the reflected ray with respect to the direct ray

Ymtray=k*circshift(Y,[0 t1]);

% Multipath
YMT= Y+Ymtray; 

% Calculate the Fourier transform of the signal
N = length(YMT);
fs4 =5 ; 
frequencies = linspace(-fs4 / 2, fs4 / 2, N);


% Calculate the FFT of the modulated signal
spectrum4 = fftshift(fft(YMT, N));

% ------------PLOTS--------------
figure;
subplot(3, 1, 1);
plot(t,Y);
title('Señal con ruido AWGN - Rayo Directa');
xlabel('tiempo (t)')
ylabel('Amplitud')
xlim([0 100]) 

subplot(3, 1, 2);
plot(t,Ymtray);
title('Señal con ruido AWGN - Rayo desviado');
xlabel('tiempo (t)')
ylabel('Amplitud')
xlim([0 100]) 

subplot(3, 1, 3);
plot(frequencies, abs(spectrum4)); % Escala en dB
title('Espectro de la Señal en canal Multitrayecto');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
grid on;

%% EQUALIZER

YMT1=circshift(YMT,[0 t1]);
YMTk1= k*YMT1;
YMTneg=-YMTk1;
Ysal1=YMTneg+YMT;

YMT2=circshift(YMT1,[0 t1]);
YMTk2=(k^2)*YMT2;
Ysal2=Ysal1+YMTk2;

% Calculate the Fourier transform of the signal
N = length(Ysal2);
fs5 =5 ; 
frequencies = linspace(-fs5 / 2, fs5 / 2, N);

% Calculate the FFT of the equalized signal
spectrum5 = fftshift(fft(Ysal2, N));

% ------------PLOTS--------------

figure;
subplot(3, 1, 1);
plot(t, YMT); 
title('Señal con multitrayecto YMT');
xlabel('tiempo (t)');
ylabel('Amplitud');
xlim([0 100]);

subplot(3, 1, 2);
plot(t, Ysal2);
title('Señal ecualizada');
xlabel('tiempo (t)');
ylabel('Amplitud');
xlim([0 100]);

subplot(3, 1, 3);
plot(frequencies, abs(spectrum5)); 
title('Espectro de la Señal ecualizada');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
grid on;
%% PASS-BAND DEMODULATION

% Initialize variables to contain the baseband signals.
Y1=sqrt(2)*Ysal2.*cos(2*pi*fc.*t); % In-phase demodulation
Y2=-sqrt(2)*Ysal2.*sin(2*pi*fc.*t); % Quadrature demodulation

%% COUPLING FILTER

Frx=F; 

% Filtering process of the signal at the channel output
Y11=filter(Frx,1,Y1);
Y21=filter(Frx,1,Y2);

% Sampling the filtered signal at the receiverr
Y12=downsample(Y11,sps);
Y22=downsample(Y21,sps);


%% BASE-BAND DEMODULATION

% Received symbols 
figure;
scatter(Y12,Y22,'filled'), title('Constelación de símbolos recibidos')
set(gcf,'color','w');
xlim([-5 5]), ylim([-5 5])
grid on

% Demapping - minimum distance criterion
Mpam=sqrt(M);

A=zeros(1,Mpam);
for i=1:Mpam
A (1,i)=(2*i-1-Mpam)*d/2;
end

T=Mpam-1; 
Tes=zeros(1,T);
for i=1:T
Tes(1,1)=(A(1,1)+A(1,1+1))/2;
end
Se1=Y12; Se2=Y22;

if M>2
    for i=2:(Mpam-2)
        Se1(Se1>Tes(1,i-1) & Se1<Tes(1,i))=A(1,i);
        Se1(Se1>Tes(1,i) & Se1 <Tes(1,i+1))=A(1,i+1);
        Se2(Se2>Tes(1,i-1) & Se2<Tes(1,i))=A(1,i);
        Se2(Se2>Tes(1,i) & Se2<Tes(1,i+1))=A(1,i+1);
    end
end

Se1=Se1(2*span+1:length(Se1));
Se2=Se2(2*span+1:length(Se2));

Sreal=Sreal(1:length(Sreal)-2*span); 
Simg=Simg(1:length(Simg)-2*span);

rx_qamsymbo = Se1 + 1j*Se2; 

%%QAM demodulation
demodulated_symbols = qamdemod(rx_qamsymbo, M);
rx_bits = de2bi(demodulated_symbols,num_bits_por_simbolo)';
dataOut = reshape(rx_bits(1:length(dataIn)), size(dataIn));

% Calculate the Fourier transform of the signal
N_dem = length(rx_qamsymbo); 
fs6 =1/Ts ; 
frequencies = linspace(-fs6 / 2, fs6 / 2, N_dem);

% Calculate the FFT of the demodulated signal
spectrum6 = fftshift(fft(rx_qamsymbo, N_dem));

figure;
subplot(2, 1, 1);
scatter(real(rx_qamsymbo), imag(rx_qamsymbo), 'filled');
title('Constelación de símbolos en recepción');
xlim([-5 5]), ylim([-5 5]);
grid on;

subplot(2, 1, 2);
plot(frequencies, abs(spectrum6)); 
title('Espectro de la Señal demodulada');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
grid on;


%% EYE DIAGRAM
eyediagram(rx_qamsymbo, span*2, span*2); 
title('Diagrama de Ojo después de la Demodulación');
%% RECOVER IMAGE

[row, col] = size(imagen_binaria); 
ReconstructedImage = reshape(dataOut, row, col);
ReconstructedImage = logical(ReconstructedImage); 

% ------------PLOTS--------------

figure;
imwrite(ReconstructedImage, 'imagen_reconstruida.jpg'); 
FINAL= imread('imagen_reconstruida.jpg');
image(FINAL);
axis image , colormap gray
title('Imagen Reconstruida');

% Calculate the Bit Error Rate (BER)
num_errors = sum(sum(dataIn ~= dataOut));
ber = num_errors / numel(dataIn);
disp(['Bit Error Rate (BER): ', num2str(ber)]);

