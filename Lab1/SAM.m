function SAM()
sam1= dirac(10,20);
figure(1); stem(sam1);
xlabel('n')
ylabel('$(k)')
title(' SAM1-PLOT DIRAC')

sam2=step(10,20);
figure(2);stem(sam2);
xlabel('n')
ylabel('H(K)=1 ; k>=0')
title(' SAM2-PLOT STEP')

sam3=ramp(2,10,20);
figure(3); stem(sam3);
xlabel('n')
ylabel(' P(K)=K IF K>=0')
title('SAM3-PLOT RAMP')

sam4= geometric(2,10,20);
figure(4);stem(sam4);
xlabel('n')
ylabel(' G(K)=a^K IF K>=0')
title('SAM4-PLOT GEOMETRIC')

sam5= box(3,10,20);
figure(5);stem(sam5);
xlabel('n')
ylabel(' B(K)=1 FOR -a<=k<=a')
title('SAM5-PLOT BOX FUNCTION')


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

sam9 = rand_a(1000); 
% observation xn 1000
% 4 is the mean of the data and 400 is the standard deviation of data
% I am plotting the normal distribution here
figure(13); histfit(sam9)
title('Fit of the Histogram Distribution for n=1000')

sam10 = rand_a2(10000);
figure(14); histfit(sam10)
title('Fit of the Histogram Distribution for n=10000')

sam11=uniform(1000);
figure(15); histfit(sam11)
title('Fit of the Uniform Distribution for n=1000')

sam12=cogaussian();
figure(16); plot(sam12)
title('Autocorelation of Gausain Noise')

sam13=couniform();
figure(17); plot(sam13)
title('Autocorelation of Uniform Noise')
end

% QUESTION-1
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

function sam3= ramp(a,n,N)

if(n>N-1)
    
        Disp(' There is Error because n cant be greater then N-1');
        
else
   b = zeros(1,N);  % 
            
            for j = n:N
            b(j) = a*(j-n);
            end
            sam3 = b;
    
end
end

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

p1 = -.5 * ((xh - mu)/s) .^ 2;
p2 = (s * sqrt(2*pi));
samf = exp(p1) ./ p2;

figure(11); plot(samf);
title('Gaussian Function')

figure(12); hist(samf);
title('Histogram of the Gaussian Function')


end

function sam10 =  rand_a2(n)

sam10 = randn(n,1);

% Increase the number of n to 10000
end

function sam11 =  uniform(n)

sam11 = rand(n,1);

end

function sam12 =  cogaussian()

sam12 = randn(1000,1);
sam12=xcorr(sam12,'biased');

end
function sam13 =  couniform()

sam13 = rand(1000,1);
sam13=xcorr(sam13,'biased');

end