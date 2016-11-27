% SAM 
% LAB4
% SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lab4samupdated()
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
% Exercice 3 
barcode();
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

 function [] =barcode()
% 
% Question-3
% I forgot to upload Earlier
% Completed by SAM AND SAED

folder_Im='C:\Users\Usama\Desktop\DSP-TP-1617\Lab4\images\1D-DFT';

cd (folder_Im)
files=dir;
Images_number=length(files)-2; 
 for i = 1:Images_number
     cell{i} = imread(files(i+2).name);
 end
%3.2
[y,x] = size(cell{1});  
 for idx = 2:Images_number
     [a,b]=size(cell{idx}) 
     if a*b<x*y
     y=a;
     x=b;
     end
 end
 p=1;q=1;
 for i=1:Images_number
    I = imread(files(i+2).name);
    I = double(I);
    [a,b,c]= size(I);
        if c==4 
            I=I(:,:,1);
        end
    I=I/(max(I(:)));  
    Resized_image = imresize(I, [y x]);   
% 3.3   This part of the code is prepared by Saed and Sam 
    z=round(y/2); 
    Profile_1 =Resized_image(z,:);
    N=x;
    fr = (-N/2 : N/2-1);
    x1 = ifftshift(fft(Profile_1)); 
    Profile_1D2 = Resized_image(z+1,:);
    Profile_1D3 = Resized_image(z-1,:); 
    Profile_1D4 = Resized_image(z+2,:); 
    Profile_1D5 = Resized_image(z-2,:); 
    x2 = ifftshift(fft(Profile_1D2));
    x3 = ifftshift(fft(Profile_1D3));
    x4 = ifftshift(fft(Profile_1D4));
    x5 = ifftshift(fft(Profile_1D5)); 
% 3.4
    v1=abs(x1);
    v2=abs(x2);
    v3=abs(x3);
    v4=abs(x4); 
    v5=abs(x5); 
    threshold1=max(v1);
     threshold2=abs(4*max(v1)-max(v2)-max(v3)-max(v4)-max(v5));
    threshold(i) = threshold1*threshold2;
    if (threshold(i)< 40)   
        Barcode_image(p)=i ;
        p=p+1; 
    else
        NonBarcode_image(q)=i; 
        q=q+1;
    end
 end
 figure (22);
 subplot(221);stem(threshold); title('Threshold');
 disp('Barcode');
 subplot(222); stem(Barcode_image); title('Non Barcode Image');
 disp('Non barcode');
 subplot(223);stem(NonBarcode_image); title('Barcode Image');
 Barcode=[1,2,6,44:54]; 
 NonBarcode = [3,4,5,7:43];
 x=length(Barcode);
 u=0;v=0;
 for i=1:length(Barcode)
   for j=1:length(Barcode_image)
     if Barcode_image(j)==Barcode(i)
         u=u+1;
     end
   end
 end
 
 for i=1:length(NonBarcode)
   for j=1:length(Barcode_image)
     if Barcode_image(j)==NonBarcode(i)
         v=v+1;
     end
   end
 end
 v=v+x-u;
 percentage_of_accuracy= ((Images_number-v)/Images_number)*100;
 sprintf('Percentage of correct distinction between the Barcode and NonBarcode images is: %d percent', round(percentage_of_accuracy)) % 85 % accurcy
 end
 end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%