function SAM()
% Question#1
% SUBMITTED BY : SAM ( USAMA PERVAIZ)
% MAIA,LAB1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.1
% MATLAB DIRAC FUNCTION
sam1= dirac(10,20);
figure(1); stem(sam1);
xlabel('n')
ylabel('$(k)')
title(' SAM1-PLOT DIRAC')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.2
% STEP FUNCTION
sam2=step(10,20);
figure(2);stem(sam2);
xlabel('n')
ylabel('H(K)=1 ; k>=0')
title(' SAM2-PLOT STEP')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.3
% RAMP FUNCTION
sam3=ramp(2,10,20);
figure(3); stem(sam3);
xlabel('n')
ylabel(' P(K)=K IF K>=0')
title('SAM3-PLOT RAMP')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.4
% GEOMETRIC FUNCTION 
sam4= geometric(2,10,20);
figure(4);stem(sam4);
xlabel('n')
ylabel(' G(K)=a^K IF K>=0')
title('SAM4-PLOT GEOMETRIC')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.5
% BOX FUNCTION
sam5= box(3,10,20);
figure(5);stem(sam5);
xlabel('n')
ylabel(' B(K)=1 FOR -a<=k<=a')
title('SAM5-PLOT BOX FUNCTION')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.6
% SINE FUNCTION
sam6=sinewave(10,1/100);
figure(6);stem(sam6);
xlabel('f=10 , fs=100')
ylabel(' sin(2*pi*f*n*t')
title('SAM6-Sin Function')

sam7=sinewave2(10,1/1000);
figure(7);stem(sam7);
xlabel('f = 10Hz and Ts = 1000 and Periods = 2')
ylabel(' sin(2*pi*f*n*t')
title('SAM7-Sin Function')

sam8=sinewave3(10,1/30);
figure(8);stem(sam8);
xlabel('f = 10, Ts = 30 over 2 periods.')
ylabel(' sin(2*pi*f*n*t')
title('SAM8-Sin Function')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1
% Histogram of the gaussian distribution
sam9 = rand_a(1000); 
% observation xn 1000
figure(13); histfit(sam9)
title('Fit of the Histogram Distribution for n=1000')

sam10 = rand_a2(10000);
figure(14); histfit(sam10)
title('Fit of the Histogram Distribution for n=10000')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2
% Same question as of 2.1 with the uniform law of the random process U and an observation xu.

sam11=uniform(1000);
figure(15); histfit(sam11)
title('Fit of the Uniform Distribution for n=1000')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.3
% Compute the autocorrelation of the two observations and plot them. Are these noises "white" ?

sam12=cogaussian();
figure(16); plot(sam12)
title('Autocorelation of Gausain Noise')

sam13=couniform();
figure(17); plot(sam13)
title('Autocorelation of Uniform Noise')
% Conclusion : The gaussian noise is the white noise as you can see from the graph that it is not dependent 
% on the frequency and the energy level is the same whereas in the case of uniform noise , it's dependent on
% frequency and it's energy distribution keep changing 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4
[sam14,sam15,sam16,sam17]=binary();
figure(18); plot(sam14);
title('Plot of the Signal S=s1+s2+s3')

figure(19); plot(sam15);
title('Corelation of the S & S1')

figure(20); plot(sam16);
title('Corelation of the S & S2')

figure(21); plot(sam17);
title('Corelation of the S & S3')


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTION-1
% 1.1
% By default n and N should be equal to 0 and 20, respectively. Plot the signal for n = 10.
function sam1= dirac(n,N)

if(n>N-1)
        Disp(' There is Error because n cant be greater then N-1');
        
else
    a=zeros(N,1);
    a(n)=1;
    sam1=a;
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2
function sam2= step(n,N)

if(n>N-1)
    
        Disp(' There is Error because n cant be greater then N-1');
        
else
    a=zeros(N,1);
    for i=n:N
        a(i)=1;
    
    end
    sam2=a;
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.3
function sam3= ramp(a,n,N)

if(n>N-1)
    
        Disp(' There is Error because n cant be greater then N-1');
        
else
            b = zeros(1,N);  
            
            for j = n:N
            b(j) = a*(j-n);
            end
            sam3 = b;
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.4
function sam4= geometric(a,n,N)

if(n>N-1)
    
        Disp(' There is Error because n cant be greater then N-1');
        
else
            b = zeros(1,N); 
            
            for j = n:N
            b(j) = a^(j-n);
            end
            sam4 = b;
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.5
function sam5= box(a,n,N)

if n >(N-a) || n< (1+a)
            Disp('There is Error because n should be inferior then N-a or greater than 1+a'); 
    else
            
            f = zeros(1,N);  %
            
            for j = n-a:n+a
            f(j) = 1;
            end
            
                
            sam5 = f;
            
end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.6
function sam6 =  sinewave(f,Ts)  
                    
%  The f=10, so time period for single cycle is 0.1 second 
   StopTime = .10;
   t = (0:Ts:StopTime-Ts);                   
   sam6 = sin(2*pi*f*t);

end

function sam7 =  sinewave2(f,Ts)  
                    
%  The f=10, so time period for single cycle is 0.1 second 
%  We want to plot two cycles of the sinewaveform
   StopTime = .20;
   t = (0:Ts:StopTime-Ts);                   
   sam7 = sin(2*pi*f*t);

end

function sam8 =  sinewave3(f,Ts)  
                    
%  The f=10, so time period for single cycle is 0.1 second 
%  We want to plot two cycles of the sinewaveform
   StopTime = .20;
   t = (0:Ts:StopTime-Ts);                   
   sam8 = sin(2*pi*f*t);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1
% 'Fit of the Histogram Distribution for n=1000'
function sam9 =  rand_a(n) %Function Defination 

sam9 = randn(n,1);
y=sam9;
[h,xh]=hist(y);
figure(9); plot(y);
title('Plot the Random Function')

figure(10); hist(y);
title('Histogram of the Random Function')

mu= mean (y);
s= std (y);

%HISTFIT FUNCTION DO ALL THIS WORK SO THERE IS NO NEED TO
%DO THAT AGAIN BY GAUSSIAN FUNCTION , SO YOU CAN COMMENT THE
%CODE AND JUST USE THE HISTFIT FUNCTION
p1 = -.5 * ((xh - mu)/s) .^ 2;
p2 = (s * sqrt(2*pi));
samf = exp(p1) ./ p2;

figure(11); plot(samf);
title('Plot of the Gaussian Function')

figure(12); hist(samf);
title('Histogram of the Gaussian Function')


end

function sam10 =  rand_a2(n)
% Fit of the Histogram Distribution for n=10000
sam10 = randn(n,1);

% Increase the number of n to 10000
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2
% Fit of the Uniform Distribution for n=1000
function sam11 =  uniform(n)

sam11 = rand(n,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.3
% Autocorelation of Gausain Noise
function sam12 =  cogaussian()

sam12 = randn(1000,1);
sam12=xcorr(sam12,'biased');
% Use of the xcorr function for calculating the coorelation
end
% 'Autocorelation of Uniform Noise'
function sam13 =  couniform()

sam13 = rand(1000,1);
sam13=xcorr(sam13,'biased');
% Conclusion : The gaussian noise is the white noise as you can see from the graph that it is not dependent 
% on the frequency and the energy level is the same whereas in the case of uniform noise , it's dependent on
% frequency and it's energy distribution keep changing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4
% Generate a whole signal s containing these signals at dirrerent shifts
function [sam14,sam15,sam16,sam17] =  binary()

s1=round(rand(1,50));
s2=round(rand(1,50));
s3=round(rand(1,50));
a=zeros(300,1);
    for i= 1:50
        a(i)=s1(i);   
    end
    for i=51 :300
        s1(i)=0;
    end
    
    for i=1:50
        a(i+100)=s2(i);
    end
    for i=51 :300
        s2(i)=0;
    end
    
    for i=1:50
        a(i+200)=s3(i);
    end
    for i=51 :300
        s3(i)=0;
    end
 
sam14=a;

length(sam14)
length(s1)
% Compute the cross-correlation
% between the whole signals and s1; s2; s3. Comments the results.
sam15=corrcoef(sam14,s1); 
sam16=corrcoef(sam14,s2);
sam17=corrcoef(sam14,s3);

% EXPLAIN??
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

