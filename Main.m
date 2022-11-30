clear all;
close all;
clc; 

% -----------------------------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------------------------
noOfSamples = 2048 ;				                                       % total number of samples
noSymb = ceil(noOfSamples/16) ;	                                           % total number of symbols (16 samples per symbol)
symbolRate = 25 * 1e9 ;	
samplingRate = 16 * symbolRate ; 	                                       % 16 samples per symbol
T = 1/samplingRate;				                                           % sampling period

z = input("Enter the distance of propogation in kilometers")*1000;         % distance of propagation (meters)
D = 1e-5	;						                                       % dispersion coefficient of the fiber
lamda = 1550e-9;					                                       % wavelength of the light
c = 3e8;    						                                       % speed of light(m/s)

[F,I]=modf(abs(D)*lamda*lamda*z/(2*c*T*T));
N = 2*I+1;							                                       % total number of taps for each case of fiber length.


% -------------------------------------------------------------------------------------------
% ALGORITHM IMPLIMENTATION
% -------------------------------------------------------------------------------------------
bitsSeq = randi([0 1],1,noSymb)	;	                                       %generates a random sequence of bits
f = linspace(-samplingRate/2, samplingRate/2, noOfSamples);

% Computation of input NRZ train pulses
win = [0,0,1,1,1,1,1,0,0];                                                 % Window function 	
squareNRZ = NRZG(bitsSeq,noOfSamples,floor(samplingRate/symbolRate));
INRZP=conv(squareNRZ, win)/sum(win);
inputNRZ = conv(ZPAD(squareNRZ,4,4), win)/sum(win);

% FFT of NRZ 
S=fft(inputNRZ);

% Required arrays
k = zeros(1,floor(N));
b = zeros(1,floor(N), 'like',i);
H = zeros(1,noOfSamples, 'like',i);
ss = zeros(1,noOfSamples, 'like',i);
SS = zeros(1,noOfSamples, 'like',i);
Fss = zeros(1,noOfSamples, 'like',i);

%Computation of the FIR coefficients.
[F,I]= modf(N/2);
k = sort(I, I+1);
b = sqrt(1j*c*T*T/(D*lamda*lamda*z)) * exp((-1j*3.14*c*T*T*k*k/(D*lamda*lamda*z)));

%Computation of the chromatic dispersion
H = ifftshift(exp(i*D*lamda*lamda*z*((2*pi*f).^2)/(4*pi*c)));

%Computation of the NRZ pulses with chromatic dispersion
SS = conv(H,S);
ss = ifft(SS);

%Counterbalance the chromatic dispersion with FIR filter.
Fss = conv(ss, b);

% -------------------------------------------------------------------------------------------
% Plotting the signals
% -------------------------------------------------------------------------------------------
% plotting the bit sequence
figure(1)
plot(linspace(0,noOfSamples*T,length(bitsSeq)),bitsSeq)
title("Input bit-sequence")

% plotting the windowed NRZ
figure(2)
plot(linspace(0,noOfSamples*T,length(inputNRZ)),inputNRZ)
title("Input NRZ with windowing")

% plotting NRZ with dispersion
figure(3)
plot(linspace(0,noOfSamples*T,length(ss)),ss)
title("input signal with dispersion")

% plotting the Compensated result 
figure(4)
plot(linspace(0,noOfSamples*T,length(Fss)),Fss)
title("signal after compensation")

% -------------------------------------------------------------------------------------------
% Useful Functions
% -------------------------------------------------------------------------------------------

% Extracting integer and fractional parts of the number
function [F,I]=modf(k)
    I=floor(k);
    F=k-I;
end

% NRZ generator
function OP=NRZG(BS,L,SPS)
    OP=[];
    for i=1:length(BS)
        for j=1:SPS
            OP(i)=BS(i);
        end
    end
end

% Zero padding
function O=ZPAD(A,x,y)
            O=[zeros(1,x) A zeros(1,y)];
end
