function SAM()
% Question#1
% SUBMITTED BY : SAM ( USAMA PERVAIZ)
% MAIA,LAB2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%REMINDER
% STEP FUNCTION
sam1=step(4,20);
figure(1);stem(sam1);
xlabel('n')
ylabel('H(K)=1 ; k>=0')
title(' SAM1-PLOT STEP')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REMINDER
% SINE FUNCTION
sam2=sinewave(1,20);
figure(2);plot(sam2);
xlabel('f=1 , fs=20')
ylabel(' sin(2*pi*f*n*t')
title('SAM2-Sin Function')

sam3=sinewave2(1,20);
figure(3);plot(sam3);
xlabel('f=1 , fs=20')
ylabel(' sin(2*pi*(f/fs)*n')
title('SAM3-Sin Function')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXERCISE 1
% CAUSALITY
% 1.1
for k= 1:19
    sam4(k) = sam1(k)/2+ (sam1(k+1))/2;
end
figure(4); stem(sam4)
xlabel('K=4, K=20 ')
ylabel('Y(K)')
title('SAM4-NON CAUSAL SYSTEM')
%Comments:
%The system contains k+1 in the its system equation which makes it 
%non-causal because it contain future values and If a system whose present 
%response depends on future values of the inputs is called as a non-causal system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2
for k= 2:20
    sam5(k) = sam1(k)/2+ (sam1(k-1))/2;
end
figure(5); stem(sam5)
xlabel('K=4, K=20')
ylabel('Y(K)')
title('SAM5-CaAUSAL SYSTEM')
%Comments:
%A system is said to be causal system if its output depends on present and
%past inputs only not on future inputs, The above system bacame
%causal by replacing k+1 by k-1; because the system does not depend on
%future inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXERCISE 2
% STABILITY
% 2.1
sam6 = 0; 
for i = 2:1:20
    sam6(i) = sam6(i-1)+sam1(i);
end
figure(6)
stem(sam6)
xlabel('K=4, K=20')
ylabel('Y(K)')
title('SAM6-ACCUMULATION')
%Comments :
%Unstable system because it keep adding previous values till our range and
%it doesn't approach to zero, The signal keep increasing with constant
%addition of 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.2A
% MATLAB DIRAC FUNCTION
sam7= dirac(4,20);
figure(7); stem(sam7);
xlabel('n')
ylabel('$(k)')
title(' SAM7-PLOT DIRAC')
%Comments :
%As we are using a dirac function so it is 1 til all its range which makes
%it stable and constant for all range of values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.2B
%Same question as of 1.1 with dirac as an input.
d= sam7;
sam8 = 0;
for i = 2:1:20
    sam8(i) = sam8(i-1)+d(i);
end
figure(8)
stem(sam8) 
xlabel('n')
ylabel('Y(K)')
title('SAM8-STABLE SYSTEM')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.3
d= sam7;
sam9 = 0;
for i = 2:1:20
    sam9(i) = d(i)+ 2*(sam9(i-1));
end

figure(9)
stem(sam9);
xlabel('n ')
ylabel('Y(K)')
title('SAM9-Unstable System')
%Comments :
%Signal Exponential increasing so its not a stable system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2.4
d= sam7;
sam10 = 0;
for i = 2:1:20
      sam10(i) = d(i)+ (sam10(i-1)/3);
end
figure(10)
stem(sam10) 
xlabel('n')
ylabel('Y(K)')
title('SAM10-STABLE SYSTEM')
%Comments:
%Signal is increasing in negative direction approaching
%to zero which makes it a stable sytem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXERCISE 3
% INVARIANCE AND LINEARITY
% 3.1
samxa=[0 0 0 0 1 2 3 4 5 0 0 0 0 0 0 0 0 0 0];
samxb=[0 0 0 0 0 0 0 0 0 4 3 2 1 0 0 0 0 0 0];
sama(1)=0;
samb(1)=0;

for i=2:1:19-1
    sama(i)=3*samxa(i-1)-2*samxa(i)+samxa(i+1);
end
figure(11);
stem(sama) 
xlabel('x(a) ')
ylabel('SAM-f(a)')
title('System SAM-Xa')

for i=2:1:19-1
    samb(i)=3*samxb(i-1)-2*samxb(i)+samxb(i+1);
end
figure(12)
stem(samb) 
xlabel('x(b) ')
ylabel('SAM-f(b)')
title('System SAM-Xb')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3.3a
h=[1,-4,8];
f=samxa+samxb; 
f1=conv(f,h);

figure(13)
stem(f1);
xlabel('X = conv(samxa) + conv(samxb) ')
ylabel('f)')
title('Convolution of two signal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3.3b
f2=conv(samxa,h)+conv(samxb,h); 
figure(14)
stem(f1);
xlabel('X = conv(samxa + samxb) ')
ylabel('f)')
title('Convolution of two signal')
%Comments:
%LTI system can be characterized entirely by a single function called the 
%system's impulse response , So we take impulse response h and check
%whether the sytem is LTI System or not, So From results we can clearly see
%that there is no change in the sytem,SO sytem in LTI system

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
% 1.6
function sam6 =  sinewave(f,fs)  
                    
  t = [0:1/fs:20];                  
   sam6 = sin(2*pi*f*t);

end

function sam7 =  sinewave2(f,fs)

   t = [0:1/fs:20];
   sam7 = sin(2*pi*(f/fs)*t);

end


