close all;
clear all;
clc;

%% generate the NB signals
signal_1 = uplinkwave()';
signal_2 = uplinkwave()';
signal_3 = uplinkwave()';
signal_4 = uplinkwave()';
signal_5 = uplinkwave()';
signal_6 = uplinkwave()';
signal_7 = uplinkwave()';
signal_8 = uplinkwave()';


fs_old = 1.92e06; % the orignal freq
% this is the first f axis for the signals before the upsampling
f = linspace(-fs_old/2 , fs_old/2 , length(signal_1));

%___________this is the plot for the orignal signal_______________________
figure;
plot(f/1000 , fftshift(abs(fft(signal_1))) );
title('original generated signal');
xlabel('freq');
ylabel('amplitude');


%% upsampling by 2


% _________________upsampling the signals by 2____________________________
up_s_signal_1 = upsample(signal_1 , 2);
up_s_signal_2 = upsample(signal_2 , 2);
up_s_signal_3 = upsample(signal_3 , 2);
up_s_signal_4 = upsample(signal_4 , 2);
up_s_signal_5 = upsample(signal_5 , 2);
up_s_signal_6 = upsample(signal_6 , 2);
up_s_signal_7 = upsample(signal_7 , 2);
up_s_signal_8 = upsample(signal_8 , 2);
%************************************************************************

%_________________freq axis for the upsampled signal_____________________
fs_new = 2 * fs_old;
f_axis = linspace(-fs_new/2 , fs_new/2 , length(up_s_signal_1));
%************************************************************************

figure;
plot(f_axis/1000 , fftshift(abs(fft(up_s_signal_1))));
title('the up sampled by 2');
xlabel('freq');
ylabel('amplitude');

%% Filtering the Signals after upsampling by 2

% __________load the filer coff for the w lease squer filter_____________
load H_Coff_Fliped;
H_Coff_Fliped = H_Coff_Fliped';
%************************************************************************



%________________filter the signal after upsampling by 2___________________
signal_1_after_upandfilter = conv(up_s_signal_1 , H_Coff_Fliped);
signal_2_after_upandfilter = conv(up_s_signal_2 , H_Coff_Fliped);
signal_3_after_upandfilter = conv(up_s_signal_3 , H_Coff_Fliped);
signal_4_after_upandfilter = conv(up_s_signal_4 , H_Coff_Fliped);
signal_5_after_upandfilter = conv(up_s_signal_5 , H_Coff_Fliped);
signal_6_after_upandfilter = conv(up_s_signal_6 , H_Coff_Fliped);
signal_7_after_upandfilter = conv(up_s_signal_7 , H_Coff_Fliped);
signal_8_after_upandfilter = conv(up_s_signal_8 , H_Coff_Fliped);
%**************************************************************************
f_axis_2 = linspace(-fs_new/2 , fs_new/2 , length(signal_1_after_upandfilter));

figure;
plot(f_axis_2/1000 , fftshift(abs(fft(signal_1_after_upandfilter))));
title('filtered signal 1');
xlabel('freq');
ylabel('amplitude');

%*******************************************************************************
%% applying the carriers IF 

% _________________________parameters for carriers_________________________
fs_for_carrier = 3.84e6;           % Sampling frequency (3.84 MHz)
t = (0:length(signal_1_after_upandfilter)-1) / fs_for_carrier;  % Time vector in seconds
%**************************************************************************

%______________________________carriers___________________________________
f_carrier_1 = 200e3;                 % Carrier frequency (200 kHz)
carrier_1 = exp(1j*2*pi*f_carrier_1*t);  % Complex exponential carrier

f_carrier_2 = 400e3;                 % Carrier frequency (200 kHz)
carrier_2 = exp(1j*2*pi*f_carrier_2*t);  % Complex exponential carrier


f_carrier_3 = 600e3;                 % Carrier frequency (200 kHz)
carrier_3 = exp(1j*2*pi*f_carrier_3*t);  % Complex exponential carrier

f_carrier_4 = 800e3;                 % Carrier frequency (200 kHz)
carrier_4 = exp(1j*2*pi*f_carrier_4*t);  % Complex exponential carrier

f_carrier_5 = 1000e3;                 % Carrier frequency (200 kHz)
carrier_5 = exp(1j*2*pi*f_carrier_5*t);  % Complex exponential carrier

f_carrier_6 = 1200e3;                 % Carrier frequency (200 kHz)
carrier_6 = exp(1j*2*pi*f_carrier_6*t);  % Complex exponential carrier

f_carrier_7 = 1400e3;                 % Carrier frequency (200 kHz)
carrier_7 = exp(1j*2*pi*f_carrier_7*t);  % Complex exponential carrier

f_carrier_8 = 1600e3;                 % Carrier frequency (200 kHz)
carrier_8 = exp(1j*2*pi*f_carrier_8*t);  % Complex exponential carrier

%*************************************************************************

%___________________signals after applying carrier________________________
modulated_signal_1 = signal_1_after_upandfilter .* carrier_1;
modulated_signal_2 = signal_2_after_upandfilter .* carrier_2;
modulated_signal_3 = signal_3_after_upandfilter .* carrier_3;
modulated_signal_4 = signal_4_after_upandfilter .* carrier_4;
modulated_signal_5 = signal_5_after_upandfilter .* carrier_5;
modulated_signal_6 = signal_6_after_upandfilter .* carrier_6;
modulated_signal_7 = signal_7_after_upandfilter .* carrier_7;
modulated_signal_8 = signal_8_after_upandfilter .* carrier_8;
%**************************************************************************

%___________________the last signal will be sent__________________________
Signal_to_be_sent = modulated_signal_1 + modulated_signal_2 + ...
    modulated_signal_3 + modulated_signal_4 + modulated_signal_5 ...
   + modulated_signal_6 + modulated_signal_7 + modulated_signal_8;
%**************************************************************************

figure;
plot(f_axis_2/1000 , fftshift(abs(fft(Signal_to_be_sent))));
title('Last Signal will be sent');
xlabel('freq');
ylabel('amplitude');

%**************************************************************************

%% in the receiver
%__________________________carriers_for__received__signal__________________
carrier_r_1 = exp(-1j*2*pi*f_carrier_1*t);  % Complex exponential carrier

carrier_r_2 = exp(-1j*2*pi*f_carrier_2*t);  % Complex exponential carrier

carrier_r_3 = exp(-1j*2*pi*f_carrier_3*t);  % Complex exponential carrier

carrier_r_4 = exp(-1j*2*pi*f_carrier_4*t);  % Complex exponential carrier

carrier_r_5 = exp(-1j*2*pi*f_carrier_5*t);  % Complex exponential carrier

carrier_r_6 = exp(-1j*2*pi*f_carrier_6*t);  % Complex exponential carrier

carrier_r_7 = exp(-1j*2*pi*f_carrier_7*t);  % Complex exponential carrier

carrier_r_8 = exp(-1j*2*pi*f_carrier_8*t);  % Complex exponential carrier
%*************************************************************************

%_______________________received_signals___________________________________
received_signal_1 = Signal_to_be_sent .* carrier_r_1;
received_signal_1 = conv(received_signal_1 , H_Coff_Fliped);
received_signal_1 = downsample(received_signal_1 , 2);


received_signal_2 = Signal_to_be_sent .* carrier_r_2;
received_signal_2 = conv(received_signal_2 , H_Coff_Fliped);
received_signal_2 = downsample(received_signal_2 , 2);

received_signal_3 = Signal_to_be_sent .* carrier_r_3;
received_signal_3 = conv(received_signal_3 , H_Coff_Fliped);
received_signal_3 = downsample(received_signal_3 , 2);

received_signal_4 = Signal_to_be_sent .* carrier_r_4;
received_signal_4 = conv(received_signal_4 , H_Coff_Fliped);
received_signal_4 = downsample(received_signal_4 , 2);

received_signal_5 = Signal_to_be_sent .* carrier_r_5;
received_signal_5 = conv(received_signal_5 , H_Coff_Fliped);
received_signal_5 = downsample(received_signal_5 , 2);

received_signal_6 = Signal_to_be_sent .* carrier_r_6;
received_signal_6 = conv(received_signal_6 , H_Coff_Fliped);
received_signal_6 = downsample(received_signal_6 , 2);

received_signal_7 = Signal_to_be_sent .* carrier_r_7;
received_signal_7 = conv(received_signal_7 , H_Coff_Fliped);
received_signal_7 = downsample(received_signal_7 , 2);

received_signal_8 = Signal_to_be_sent .* carrier_r_8;
received_signal_8 = conv(received_signal_8 , H_Coff_Fliped);
received_signal_8 = downsample(received_signal_8 , 2);
%*************************************************************************


%% ploting the received_signals _signal 1
f_axis_recived = linspace(-fs_old/2 , fs_old/2 , length(received_signal_1));
figure; 
plot(f_axis_recived/1000 , fftshift(abs(fft(received_signal_1))));
title('received signal 1');
xlabel('freq');
ylabel('amplitude');

%% Signal 2
figure; 
plot(f_axis_recived/1000 , fftshift(abs(fft(received_signal_2))));
title('received signal 2');
xlabel('freq');
ylabel('amplitude');

%% Signal 3
figure; 
plot(f_axis_recived/1000 , fftshift(abs(fft(received_signal_3))));
title('received signal 3');
xlabel('freq');
ylabel('amplitude');

%% Signal 4
figure; 
plot(f_axis_recived/1000 , fftshift(abs(fft(received_signal_4))));
title('received signal 4');
xlabel('freq');
ylabel('amplitude');

%% Signal 5
figure; 
plot(f_axis_recived/1000 , fftshift(abs(fft(received_signal_5))));
title('received signal 5');
xlabel('freq');
ylabel('amplitude');

%% Signal 6
figure; 
plot(f_axis_recived/1000 , fftshift(abs(fft(received_signal_6))));
title('received signal 6');
xlabel('freq');
ylabel('amplitude');

%% Signal 7
figure; 
plot(f_axis_recived/1000 , fftshift(abs(fft(received_signal_7))));
title('received signal 7');
xlabel('freq');
ylabel('amplitude');

%% Signal 8
figure; 
plot(f_axis_recived/1000 , fftshift(abs(fft(received_signal_8))));
title('received signal 8');
xlabel('freq');
ylabel('amplitude');