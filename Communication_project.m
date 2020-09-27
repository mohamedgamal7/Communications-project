 
clear;
clc;

%% Control flags
%value of the desired amplitide 
A=4;
%number of waveforms 
w=500;
%number of bits
b=100;
%sampling frequency is 100 hz  
fs=100;



%% Creating 500 waveforms each with 100 random bits 
distr=makedist('Binomial');
Data = cell(1,w);
for K = 1 : w
  Data{K} = random(distr, 1,b);
end


%% Unipolar Realizations


for K1 = 1 : w
 UNI{K1} = Data{K1}*A;
  UNI_repeated{K1}=repmat( UNI{K1},7,1);
  UNI_ready{K1}=reshape( UNI_repeated{K1},size( UNI_repeated{K1},1)* size( UNI_repeated{K1},2),1);
end

figure ;
plot(UNI_ready{1},'r');
ylim([-2 7]);
title('An example of a waveform of UNipolar line code');
xlabel('number of bits') ;
ylabel('amplitude') ;

%% Polar non return to zero Realizations

for K1 = 1 : w
 NRZ{K1} = ((2*Data{K1})-1)*A;
 NRZ_repeated{K1}=repmat(NRZ{K1},7,1);
 NRZ_ready{K1}=reshape(NRZ_repeated{K1},size(NRZ_repeated{K1},1)* size(NRZ_repeated{K1},2),1);
end

figure ;
plot(NRZ_ready{1},'k');
ylim([-9 9]);
title('An example of a waveform of polar NRZ line code');
xlabel('number of bits') ;
ylabel('amplitude') ;


%% Polar Return to zero Realizations
for K1 = 1 : w
 RZ{K1}=((2*Data{K1})-1)*A ;
 RZ_repeated{K1}=repmat (RZ{K1},4,1);
RZ_repeated{K1}=cat(1,RZ_repeated{K1},zeros(3,100));
  RZ_ready{K1}=reshape(RZ_repeated{K1},size(RZ_repeated{K1},1)* size(RZ_repeated{K1},2),1);
end
figure ;
plot(RZ_ready{1},'b');
ylim([-9 9]);
title('An example of a waveform of polar RZ line code');
xlabel('number of bits') ;
ylabel('amplitude') ;
  
%% Applying time shifts to create different initial time shifts
for K1 = 1:w 
    random_shift = randi([19 30], 1, 1);
  UNI_ready{1,K1} = circshift( UNI_ready{1,K1},random_shift);
   NRZ_ready{K1}=circshift( NRZ_ready{K1},random_shift);
  RZ_ready{K1}=circshift( RZ_ready{K1},random_shift);
end

%% Creating cell matrices where each coloumn represents a time instance of each realization
UNI=num2cell(transpose(cell2mat(UNI_ready)),1);
NRZ=num2cell(transpose(cell2mat(NRZ_ready)),1);
RZ=num2cell(transpose(cell2mat(RZ_ready)),1);



%% Q1 :compute the statistical mean

%% calculating the Statistical mean of realizations  

for i=1:700
      meanUNI{i}=mean(UNI{i});
end
%mean for RZ
for i=1:700
      meanRZ{i}=mean(RZ{i});   
end
%mean for NRZ

for i=1:700
    meanNRZ{i}=mean(NRZ{i});
    
end

%% ploting the statistical mean of each line code
figure;
subplot(3,1,1);
plot(cell2mat(meanUNI),'r');
ylim([0 5]);
title('Statistical mean of Unipolar Realizations ');
xlabel('time instance') ;
ylabel('Mean value') ;

subplot(3,1,2);
plot(cell2mat(meanNRZ),'k');
ylim([-5 5]);
title('Statistical mean of polar NRZ Realizations');
xlabel('time instance') ;
ylabel('Mean value') ;

subplot(3,1,3);
plot(cell2mat(meanNRZ),'b');
ylim([-5 5]);
title('Statistical mean of polar  RZ Realizations ');
xlabel('time instance') ;
ylabel('Mean value') ;

%% Q3 : Determine the ensemble autocorrelation function Rx(?)

%% calculating the autocorelation for the line code ensemble
 
%autocorrelation of UNIpolar
for tao=0:1:699

autoUNI{(tao)+1}=mean(UNI{1}.*UNI{1+tao});
%autoUNI(tao+1)=mean(UNI{700}.*UNI{700-tao});
end
%autocorrelation of polarNRZ
for tao=0:1:699
  autoNRZ{tao+1}=mean(NRZ{1}.*NRZ{1+tao});
  %autoNRZ(tao+1)=mean(NRZ{700}.*NRZ{700-tao});
end
%autocorrelation of polarRZ
for tao=0:1:699   
    autoRZ{tao+1}=mean(RZ{1}.*RZ{1+tao});  
    %autoRZ(tao+1)=mean(RZ{700}.*RZ{700-tao});
end
%% Ploting the statistical autocorrelation of each line code
figure;
subplot(3,1,1);
t=0:699 ;
plot(-t,cell2mat(autoUNI),'r',t,cell2mat(autoUNI),'r');
xlim([-690 690]);
title('Statistical autocorrelation of UNI');
xlabel('tao') ;
ylabel('Rx(tao)') ;

subplot(3,1,2);
t=0:699 ;
plot(-t,cell2mat(autoNRZ),'b',t,cell2mat(autoNRZ),'b');
xlim([-690 690]);
title('Statistical autocorrelation of NRZ');
xlabel('tao') ;
ylabel('Rx(tao)') ;

subplot(3,1,3);
t=0:699 ;
plot(-t,cell2mat(autoRZ),'g',t,cell2mat(autoRZ),'g');
xlim([-690 690]);
title('Statistical autocorrelation of RZ');
xlabel('tao') ;
ylabel('Rx(tao)') ;
%% Q2 : Is the process stationary ?
%since :
%The mean is constant across ensemble (1)
%The autocorelation is a function of tao only not time (2)
%therefore:
%The procees is stationary 

%% Q4 :Compute the time mean and auto correlation of one wave form 

%getting the time mean of a random wave form from any of the 500 wavefroms
RandomWaveNum=randi([1 500]);

onewaveformUNI=circshift(UNI_ready{RandomWaveNum},[1 90]);
onewaveformNRZ=circshift(NRZ_ready{RandomWaveNum},[1 90]);
onewaveformRZ=circshift(RZ_ready{RandomWaveNum},[1 90]);

tmeanuni=zeros(500,1);
tmeannrz=zeros(500,1);
tmeanrz=zeros(500,1);

for i=1:500
tmeanuni(i,1)=tmeanuni(i,1)+mean(UNI_ready{i});
tmeannrz(i,1)=tmeannrz(i,1)+mean(NRZ_ready{i});
tmeanrz(i,1)=tmeanrz(i,1)+mean(RZ_ready{i});
end

 
 figure;
 r=0:499 ;
 subplot(3,1,1);
plot(r,tmeanuni,'r');
ylim([0 5]);
title('time mean of UNIpolar Realizations ');
xlabel('realization') ;
ylabel('Mean value') ;
 r=0:499 ;
 subplot(3,1,2);
plot(r,tmeanrz,'b');
ylim([-5 5]);
title('time mean of NRZpolar Realizations ');
xlabel('realization') ;
ylabel('Mean value') ;
  r=0:499 ;
subplot(3,1,3);
plot(r,tmeanrz,'g');
ylim([-5 5]);
title('time mean of RZpolar Realizations ');
xlabel('realization') ;
ylabel('Mean value') ;


tauni=zeros(700,1);
tanrz=zeros(700,1);
tarz=zeros(700,1);
sum=0;
for tao=0:699
for x=1:700
    if((tao+x)<=700)
    y=onewaveformUNI(x,1)*onewaveformUNI(x+tao,1);
    sum=sum+y;
    end
    
    end
 sum=sum/(700-tao);
tauni(tao+1,1)=sum;
end
 
 t=0:699 ;
 figure;
 subplot(3,1,1);
plot(-t,tauni,'r',t,tauni,'r');
title('Time autocorrelation of UNipolar');
xlim([-650 650]);
xlabel('tao') ;
ylabel('Rx(tao)') ;

sum=0;
for tao=0:699
for x=1:700
    if((tao+x)<=700)
    y=onewaveformNRZ(x,1)*onewaveformNRZ(x+tao,1);
    sum=sum+y;
    end
    
    end
 sum=sum/(700-tao);
tanrz(tao+1,1)=sum;
end

 t=0:699 ;
 subplot(3,1,2);
plot(-t,tanrz,'b',t,tanrz,'b');
title('Time autocorrelation of NRZ');
xlim([-650 650]);
xlabel('tao') ;
ylabel('Rx(tao)') ;

sum=0;
for tao=0:699
for x=1:700
    if((tao+x)<=700)
    y=onewaveformRZ(x,1)*onewaveformRZ(x+tao,1);
    sum=sum+y;
    end
    
    end
 sum=sum/(700-tao);
tarz(tao+1,1)=sum;
end



 t=0:699 ;
 subplot(3,1,3);
plot(-t,tarz,'g',t,tarz,'g');
xlim([-650 650]);
title('Time autocorrelation of RZ');
xlabel('tao') ;
ylabel('Rx(tao)') ;

%% Q5: IS the random process ergodic ?? 
 % since: 
 %The time means and stastistical means are equal (1)
 %The time autocorrelation of one wave is equivalent to the 
 %the statistical autocorrelation (as it is the same waveform )(2)
 % Therefore :
 % The process is ergodic 
 

 %% PLoting the PSD of the ensemble 
    tao=700;
    %represents the size of statsical autocorrelation arrays
     figure;
 subplot(3,1,1);
 DC=mean(cell2mat(autoUNI));
 autoUNI=cell2mat(autoUNI)-DC;
 PSDUNI=fft((autoUNI));
 k=-tao/2:tao/2-1;
 hold on;
 stem(0,DC);
plot(k*fs/tao,fftshift(abs( PSDUNI)/10));
hold off;
title('power spectral density of unipolar');
  xlabel('freq(hz)') ;
  ylabel('watt/hz') ;
  axis([-50 50 0 16]);
  
 subplot(3,1,2);
 PSDNRZ=fft((cell2mat((autoNRZ))));
 k=-tao/2:tao/2-1; 
plot(k*fs/tao,fftshift(abs( PSDNRZ))/10);
title('power spectral density of Polar NRZ');
xlabel('freq(hz)') ;
ylabel('watt/hz') ;

subplot(3,1,3);
PSDRZ=fft((cell2mat((autoRZ))));
k=-tao/2:tao/2-1;
plot(k*fs/tao,fftshift(abs( PSDRZ))/10);
title('power spectral density of Polar return to zero');
xlabel('freq(hz)') ;
ylabel('watt/hz') ;

 %% Q6 : What is the bandwidth of the transmitted signal ?
 %(1)unipolar:14.2857hz 
 %(2)polarNRZ:14.2857hz
 %(3)polarRZ:25hz
 
 
 
 