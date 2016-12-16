
%%%%  Sam Pervaiz %%%%
%%%%%    LAB6&7   %%%%%

function SamLab6_7()


%  Question 1

[x,y] = butter(3,0.5,'low');
[H, W]= freqz(x,y);
figure(1);
subplot(2,2,1);
plot(W/pi, abs(H)); 
title('Low Pass Butterworth');

[x,y] = butter(3,0.5, 'high');
[H, W]= freqz(x,y);
subplot(2,2,2);
plot(W/pi, abs(H)); 
title('High Pass Butterworth');

[x,y] = butter(3,[0.6 0.9]);
[H, W]= freqz(x,y);
subplot(2,2,3);
plot(W/pi, abs(H));
title('Band Pass Butterworth');

[x,y] = butter(3,[0.6 0.9] , 'stop');
[H, W]= freqz(x,y);
subplot(2,2,4);
plot(W/pi, abs(H));
title('Band Stop Butterworth');

[x,y] = cheby1(3,1,0.5,'low');  
[H, W]= freqz(x,y);
figure(2);
subplot(2,2,1);
plot(W/pi, abs(H)); 
title('Low Pass Chebychev-I');

[x,y] = cheby1(3,1, 0.5, 'high');
[H, W]= freqz(x,y);
subplot(2,2,2);
plot(W/pi, abs(H));
title('High Pass Chebychev-I');

[x,y] = cheby1(3,1, [0.5 0.7]);
[H, W]= freqz(x,y);
subplot(2,2,3);
plot(W/pi, abs(H));
title('Band Pass Chebychev-I');

[x,y] = cheby1(3,1, [0.5 0.7] , 'stop');
[H, W]= freqz(x,y);
subplot(2,2,4);
plot(W/pi, abs(H)); 
title('Bamd Stop Chebychev-I');
figure(3);

[x,y] = cheby1(3,1,0.5,'low');
[H, W]= freqz(x,y);
subplot(2,2,1);
plot(W/pi, abs(H));
title('Low Pass Chebyshev-I, with an order of 3');

[x,y] = cheby1(5,1,0.5,'low');
[H, W]= freqz(x,y);
subplot(2,2,2);
plot(W/pi, abs(H)); 
title('Low Pass Chebyshev-I, with an order of 5');

[x,y] = cheby1(10,1,0.5,'low');
[H, W]= freqz(x,y);
subplot(2,2,3);
plot(W/pi, abs(H));
title('Low Pass Chebyshev-I, with an order of 10');

[x,y] = cheby1(20,1,0.5,'low');
[H, W]= freqz(x,y);
subplot(2,2,4);
plot(W/pi, abs(H));
title('Low Pass Chebyshev-I, with an Order of 20');

%Observations : The gradient of the Low Pass filter gets more ideal, also
% the ripple of the Low Pass filter increases by increasing the order 

% Question 2 

SamDirac=dirac(40, 20);
figure(4);   
plot(SamDirac); 
title('Dirac Function'); xlabel('k'); ylabel('x(k)');
Scale = 0.5; 
Ts = 1; 
Alpha = Scale*Ts;  
AlphaSignal = exp(-Alpha) ;
SamDiracAnticausalSmooth = zeros(length(SamDirac),1);
DiracLength=length(SamDirac)-2:-1:1;

for i =  DiracLength 
 SamDiracAnticausalSmooth(i) = Scale*Alpha*AlphaSignal*dirac(i+1)+(2*AlphaSignal)*SamDiracAnticausalSmooth(i+1)-(AlphaSignal^2)*SamDiracAnticausalSmooth(i+2);
end

figure(5)
stem (SamDiracAnticausalSmooth) ;
title('Anticausal Smoothing'); 
SamDiracCausalSmooth = zeros(length(SamDirac),1);

for i = 3:length(SamDirac);
 SamDiracCausalSmooth(i) = -Scale*Alpha*AlphaSignal*dirac(i-1)+(2*AlphaSignal)*SamDiracAnticausalSmooth(i-1)-(AlphaSignal^2)*SamDiracAnticausalSmooth(i-2);
end

figure(6)
stem (SamDiracCausalSmooth) ;
title('Causal Smmothing'); 
step10= step(40,10)
step30= step(40,30)
StepFinal=step10-step30;
figure(7)
stem(StepFinal)
StepCausal = zeros(length (StepFinal),1) ;

for i = 3 : length(StepFinal)
 StepCausal(i) = StepFinal(i)+AlphaSignal*(Alpha-1)*StepFinal(i-1)+(2*AlphaSignal)*StepCausal(i-1)-(AlphaSignal^2)*StepCausal(i-2) ;
end
figure(8)
stem (StepCausal) ;
title('Causal Deravative');
StepAnticausal = zeros(length (StepFinal),1);
StepLength = length(StepFinal)-2 : -1 : 1 ; 
for i = StepLength  
 StepAnticausal(i) = AlphaSignal*(Alpha+1)*step(i+1)-(AlphaSignal^2)*step(i+2)+(2*AlphaSignal)*StepAnticausal(i+1)-(AlphaSignal^2)*StepAnticausal(i+2) ;  
end

figure(9)
stem (StepAnticausal);
title('Anticausal Deravative');


% Question 3

SamBarbara = imread('C:\Users\Usama Perviz\Desktop\DSP-TP-1617\Lab6-7\images\barbara.gif');
figure(10);
imshow(SamBarbara);  
Case1st = zeros(size(SamBarbara)); 
Case2nd = zeros(size(SamBarbara)); 

for i = 1:size(SamBarbara, 2)
    SamImage1 = SamBarbara(:,i);
    CasualResponse = zeros(length (SamImage1),1) ;
    for i = 3 : length(SamImage1)
     CasualResponse(i) = SamImage1(i)+AlphaSignal*(Alpha-1)*SamImage1(i - 1)+(2*AlphaSignal)*CasualResponse(i-1)-(AlphaSignal^2)*CasualResponse(i-2) ;
    end
    AntiCausalResponse = zeros(length (SamImage1 ),1) ;
    BarbaraLength = length(SamImage1)-2 : -1 : 1 ;
    for i =  BarbaraLength
     AntiCausalResponse(i) = AlphaSignal*(Alpha+1)*SamImage1(i+1)-(AlphaSignal^2)*SamImage1(i+2)+(2*AlphaSignal)*AntiCausalResponse(i+1)-(AlphaSignal^2)*AntiCausalResponse(i+2) ;
    end
    TotalResponse1 = CasualResponse + AntiCausalResponse;
    Case1st(:,i) = TotalResponse1;            
end

figure(11);
imshow (Case1st, []);

for i = 1:size(SamBarbara, 2)
    image_2 = SamBarbara(i,:);
    CasualResponse = zeros(length (image_2),1) ;
    for i = 3 : length(image_2)
     CasualResponse(i) = image_2(i)+AlphaSignal*(Alpha-1)*image_2(i - 1)+(2*AlphaSignal)*CasualResponse(i-1)-(AlphaSignal^2)*CasualResponse(i-2) ;
    end
    AntiCausalResponse = zeros(length (image_2 ),1) ;
    BarbaraLength = length(image_2)-2 : -1 : 1 ;
    for i =  BarbaraLength
     AntiCausalResponse(i) = AlphaSignal*(Alpha+1)*image_2(i+1)-(AlphaSignal^2)*image_2(i+2)+(2*AlphaSignal)*AntiCausalResponse(i+1)-(AlphaSignal^2)*AntiCausalResponse(i+2) ;
    end
    TotalResponse2 = CasualResponse + AntiCausalResponse;
    Case2nd(i,:) = TotalResponse2;   
end

figure(12);
imshow (Case2nd, []);



function F1 =  dirac(n,N)       
    if ((n<1)||(n>N))
        disp('n > N-1');  
        F1= 0;
    else
        s = zeros(1,N);
        s(n) = 1 ;
        F1 = s;   
    end
end


function F1 =  step(n,N) 
    if ((n<1)||(n>N))
        disp('n > N-1'); 
        F1= 0;    
    else   
        s = zeros(N,1);   
        for i = n+1:N
            s(i) = 1 ;
        end
        F1 = s;
        figure(2)
        subplot(3,1,1)
        stem(F1) ;                        
        xlabel('X'); ylabel('Y')   
    end
end 
end