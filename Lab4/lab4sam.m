% SAM 
% LAB4
% SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lab4sam()
% I am calling my functions here in this folder
% Exercice 1 – DFT
sine();
cos4();
squar();
noise();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercice 2 – Sampling
addtwo();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1
% The DFT of a 5 Hz sin wave sampled with the sampling of fs = 50 Hz over 1000 (N = 1000)
% samples is computed
function []= sine()

f= 5; 
fs = 50;
t = 0: 1/fs : 1;
xn = sin(2 * f * t);
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(1);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');


end

% 1.2
% Compute the DFT of a cosine wave, how that differs from DFT of a sine wave

function [] =cos4()

f= 5; 
fs = 50;
t = 0: 1/fs : 1;
xn = cos(2 * f * t);
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(2);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

% cosines and sines differ only by phase shift. The difference in their performance arises from their boundary behavior.
% On an interval [0,?], the sine system {sin2?n?} satisfies the Dirichlet boundary condition, attaining zero value at 0,?.
% end

end
% 
% 1.3
% Use square wave using the same frequency and sampling frequency as the sin and cosine wave.
% (“signal.square” and “square” in python and matlab, respectively)


function [] =squar()

f= 5; 
fs = 50;
t = 0: 1/fs : 1;
xn = square(2 *pi * f * t);
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(3);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

end

% 1.4
% Use Gaussian noise with 10000, samples.
%   I am not using the gaussian wave with 10000 samples because my 
%   computer start behaving very badly.
function [] =noise()

fs = 50;
t = 0: 1/fs : 1;
xn4 = randn (100);
N = length(xn4);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn4));
figure (4)
title('Gaussian Noise')
% subplot(221); plot(t, xn4); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.1 Generate and display the following signals of 1 sec duration

% Frequencies : 10, 20, 25, 40, 50, 100, 150

function [] =addtwo()

f1= 5;
f2=20;
fs = 10;
t = 0: 1/fs : 1;
x1 = sin(2 *pi * f1 * t);
x2= cos(2 *pi * f2 * t);
xn= 4*x1 +3 *x2;
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(5);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

f1= 5;
f2=20;
fs = 20;
t = 0: 1/fs : 1;
x1 = sin(2 *pi * f1 * t);
x2= cos(2 *pi * f2 * t);
xn= 4*x1 +3 *x2;
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(6);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

f1= 5;
f2=20;
fs = 25;
t = 0: 1/fs : 1;
x1 = sin(2 *pi * f1 * t);
x2= cos(2 *pi * f2 * t);
xn= 4*x1 +3 *x2;
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(7);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

f1= 5;
f2=20;
fs = 40;
t = 0: 1/fs : 1;
x1 = sin(2 *pi * f1 * t);
x2= cos(2 *pi * f2 * t);
xn= 4*x1 +3 *x2;
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(8);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

f1= 5;
f2=20;
fs = 50;
t = 0: 1/fs : 1;
x1 = sin(2 *pi * f1 * t);
x2= cos(2 *pi * f2 * t);
xn= 4*x1 +3 *x2;
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(9);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

f1= 5;
f2=20;
fs = 100;
t = 0: 1/fs : 1;
x1 = sin(2 *pi * f1 * t);
x2= cos(2 *pi * f2 * t);
xn= 4*x1 +3 *x2;
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(10);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

f1= 5;
f2=20;
fs = 150;
t = 0: 1/fs : 1;
x1 = sin(2 *pi * f1 * t);
x2= cos(2 *pi * f2 * t);
xn= 4*x1 +3 *x2;
N = length(xn);
fr = (-N/2 : N/2-1)* fs/N;
xf = fftshift(fft(xn));
figure(11);
subplot(221); plot(t, xn); title('Signal'); xlabel('Time(sec)'); ylabel('Amplitude');
subplot(222); plot(fr, abs(xf)); title('Magnitude'); xlabel('Frequency'); ylabel('|X(f)|');
subplot(223); plot(fr, real(xf)); title('Real'); xlabel('Frequency'); ylabel('Re(X(f))');
subplot(224); plot(fr, imag(xf)); title('Imaginary'); xlabel('Frequency'); ylabel('Im(X(f))');

% Aliaisng
%  aliasing is an effect that causes different signals to become indistinguishable (or aliases of one another) 
%  when sampled. It also refers to the distortion or artifact that results when the signal reconstructed from 
%  samples is different from the original continuous signal.
 

% The maximum frequency of x[n] is 20 herth. so, when we sample the signal
 % the sampling frequency must be   GREATER THAN two times the maximum
 % frequency of the the signal to be sampled. as it is depicted from the the
 % from the DFt of the sampled signals there is no ALIASING when the fs>
 % 2*20=40 herz also called Nyquist Sampling frequency, when the sampling 
 %frequency is small the the lagrge frequency signal can't be recovered and
 %thus it is not displayed in the output
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [] =barcode()
% 
% img=
% 
% 
% 
% 
% end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%