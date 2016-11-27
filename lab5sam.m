% Lab 5
% MAIA
% USAMA PERVAIZ aka SAM

% % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
function lab5sam()

% Call All the Functions here in the main Function
%  Exercice 1  2D - DFT

% 1.1
syn();

%1.2
translated();

%1.3
rotated();

%1.4
syn2();
syn3();

%1.5
synthetic2();

%1.6
lena();

%1.7
sobel();
end
% % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function []= syn()
% 1.1
% We here compute the FFTof the Synthetic Image
img = zeros(301,301);
img(100:200, 140:160) = 255;

imgFreq = fftshift(fft2(img));
figure(1);
% We will show the Input Image , the magnitude and the phase of the FFT
subplot(131); imshow(img); title('Input Signal')
subplot(132); imagesc(abs(imgFreq)); colormap('gray'); title('Magnitude')    
subplot(133); imagesc(angle(imgFreq)/pi*180); colormap('gray'); title('Phase')

end

% % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [] = translated()
% 1.2
% We here compute the FFTof the Synthetic Image
imgTrans = zeros(301,301);
imgTrans(150:250, 160:180) = 255;

imgFreq = fftshift(fft2(imgTrans));

figure(2);
subplot(131); imshow(imgTrans); title('Input Signal')
subplot(132); imagesc(abs(imgFreq)); colormap('gray'); title('Magnitude')
subplot(133); imagesc(angle(imgFreq)/pi*180); colormap('gray'); title('Phase')

end

% % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [] = rotated()

% 1.3
% We here compute the FFTof the Synthetic Image
img = zeros(301,301);
img(100:200, 140:160) = 255;
imgRot = imrotate(img, 45);

imgFreq = fftshift(fft2(imgRot));
% We will show the Input Image , the magnitude and the phase of the FFT
figure(3);
subplot(131); imshow(imgRot); title('Input Signal')
subplot(132); imagesc(abs(imgFreq)); colormap('gray'); title('Magnitude')
subplot(133); imagesc(angle(imgFreq)/pi*180); colormap('gray'); title('Phase')

end

% Observations
% 
% Translate Image
% It has no change in magnitude but has change in phase spectrum. 
%  
% Rotated imaged
% It has change in both magnitude and phase spectrum.
 
% % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% 1.4
function [] = syn2()

img2 = zeros(301,301);
img2(20:120, 140:160) = 255;
img2(180:280, 140:160) = 255;
% We here compute the FFTof the Synthetic Image
imgFreq = fftshift(fft2(img2));

figure(4);
% We will show the Input Image , the magnitude and the phase of the FFT
subplot(131); imshow(img2); title('Input Signal')
subplot(132); imagesc(abs(imgFreq)); colormap('gray'); title('Magnitude')
subplot(133); imagesc(angle(imgFreq)/pi*180); colormap('gray'); title('Phase')

end

function [] =syn3()

img3 = zeros(301,301);
img3(100:200, 145:155) = 255 ;
% We here compute the FFTof the Synthetic Image
imgFreq = fftshift(fft2(img3));
figure(5);
% We will show the Input Image , the magnitude and the phase of the FFT
subplot(131); imshow(img3); title('Input Signal')
subplot(132); imagesc(abs(imgFreq)); colormap('gray'); title('Magnitude')
subplot(133); imagesc(angle(imgFreq)/pi*180); colormap('gray'); title('Phase')

end

% 1.5
% % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [] =synthetic2()

% We will generate the synthetic Image here
Im=0;
N=64;
T=1;
Ts=T/N;
Fs=1/Ts;
df=Fs/N;
Im(N/8:N/4,N/4+1:N/2)=1;
Im(1:N/4,N/2+1:N)=Im;
Im(N/4+1:N/2,:) = Im;
Im(N/2+1:3*N/4,:) = Im(1:N/4,:);
Im(3*N/4+1:N,:) = Im(1:N/4,:);

% We here compute the FFTof the Synthetic Image
imgFreq = fftshift(fft2(Im));
% Here is the Normalize Centered Frequency of the Image
Centered= imgFreq(N/2+1, N/2+1)/64^2
Centered2= mean (Im(:))
% We calculate the centered frequency by both means and they are equal 
figure(6);
% We will show the Input Image , the magnitude and the phase of the FFT
subplot(131); imshow(Im); title('Input Signal')
subplot(132); imagesc(abs(imgFreq)); colormap('gray'); title('Magnitude')
subplot(133); imagesc(angle(imgFreq)/pi*180); colormap('gray'); title('Phase')
%

% Plot If (u; 0) and If (0; v) with the correct frequency range. 
% Discuss your observaton.
% I will select the centered row and the centered column
 col = imgFreq(:,N/2+1); 
 row = imgFreq(N/2+1,:); 
 fr = (-N/2 : N/2-1);
 
%  I will display the magnitude of If (u; 0) and If (0; v)
 figure(7); 
 subplot(121); plot(fr,abs(row));title('Plot - If(u,0)')
 subplot(122); plot(fr,abs(col));title('Plot - If(0,v)');
  

end

% 1.6
% % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [] = lena()

% Load the \lena" image, and show the phase and magnitude, 
% then reconstruct the image using
% either frequency or phase.

Samlena= imread('/lena-grey.bmp');
% We here compute the FFTof the Lena Image
imgFreq = fftshift(fft2(Samlena));
figure(8);
% We will show the Input Image , the magnitude and the phase of the FFT
subplot(131); imshow(Samlena); title('Input Signal')
subplot(132); imagesc(abs(imgFreq)); colormap('gray'); title('Magnitude')
subplot(133); imagesc(angle(imgFreq)/pi*180); colormap('gray'); title('Phase')

% Reconstruction of the Input Image

mag1=abs(imgFreq);
s=log(1+fftshift(imgFreq));
phase1=angle(imgFreq);
r1=ifftshift(ifft2(mag1));
r2=ifft2(exp(1i*phase1));
figure(9); 

% Both by Magnitude and Phase
subplot(121);imshow(uint8(r1)); title(' Reconstruct by Magnitude') ;
subplot(122);imshow(r2,[]); title(' Reconstruct by Phase') ;
end

% 1.7
% % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [] =sobel()

% Apply the sobel filter only in vertical direction to `lena` 
% image in the frequency domain.
% sam
%  The sobel Filter
 sobelfilter=[-1 0 1; -2 0 2; -1 0 1];
 Samlena= imread('/lena-grey.bmp');
 
%  The size of lena and the sobel filter
 dimfilter = size(sobelfilter);
 diminput = size(Samlena);
%  To calculate the size for the Zer-Padding
 ZP = diminput + dimfilter - 1;

%  Compute the DFT of the image with additional shift
 dftInput = fft2(double(Samlena), ZP(1), ZP(2));
 dftFilter = fft2(double(sobelfilter), ZP(1), ZP(2));
 
%  Apply the multiplication in the Fourier space
 dftOutput = dftInput.*dftFilter;
%  
%  Compute the inverse Fourier transform
 InverseOutput = ifft2(dftOutput);
 
%  Crop the image at its original size
Output = InverseOutput(2:size(Samlena,1)+1, 2:size(Samlena,2)+1);
figure(10);
% Dispaly the Results
subplot(121);imshow(Samlena); title('Original Sam Lena') ;
subplot(122); imshow(Output,[]); title(' Sam Lena After Filtering') ;

% Observations
% We detect the edges of the lena Image
end
