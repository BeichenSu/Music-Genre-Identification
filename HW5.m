clear all; close all;clc;
% load music file of 5 second
% reduce them with one row, that is one sound track
% Guns N Roses 
% Sweet child o mine
tr = 5; % time in second 
%% Rock music
% Guns N roses: paradise city, nightrain, sweet child o mine
% each with 5 sample clips
% sweet child o mine

som1 = wavread('GNRSOM1');

som2 = wavread('GNRSOM2');
som3 = wavread('GNRSOM3');
som4 = wavread('GNRSOM4');
som5 = wavread('GNRSOM5');

% paradise city
pc1 = wavread('GNRPC1');
pc2 = wavread('GNRPC2');
pc3 = wavread('GNRPC3');
pc4 = wavread('GNRPC4');
pc5 = wavread('GNRPC5');

% nightrain
nt1 = wavread('GNRNT1');
nt2 = wavread('GNRNT2');
nt3 = wavread('GNRNT3');
nt4 = wavread('GNRNT4');
nt5 = wavread('GNRNT5');

% reduce them into mono sound track
% rescale the data into smaller size
scale = 100;
som1 = som1(:,1);
som1 = decimate(som1,scale);
som2 = som2(:,1);
som2 = decimate(som2,scale);
som3 = som3(:,1);
som3 = decimate(som3,scale);
som4 = som4(:,1);
som4 = decimate(som4,scale);
som5 = som5(:,1);
som5 = decimate(som5,scale);

pc1 = pc1(:,1);
pc1 = decimate(pc1,scale);
pc2 = pc2(:,1);
pc2 = decimate(pc2,scale);
pc3 = pc3(:,1);
pc3 = decimate(pc3,scale);
pc4 = pc4(:,1);
pc4 = decimate(pc4,scale);
pc5 = pc5(:,1);
pc5 = decimate(pc5,scale);

nt1 = nt1(:,1);
nt1 = decimate(nt1,scale);
nt2 = nt2(:,1);
nt2 = decimate(nt2,scale);
nt3 = nt3(:,1);
nt3 = decimate(nt3,scale);
nt4 = nt4(:,1);
nt4 = decimate(nt4,scale);
nt5 = nt5(:,1);
nt5 = decimate(nt5,scale);

som = [som1 som2 som3 som4 som5];
pc = [pc1 pc2 pc3 pc4 pc5];
nt = [nt1 nt2 nt3 nt4 nt5];
%% pop music 
% jay chou : qinghuaci, dongfengpo,
qhc1 = wavread('Jay_ChouQHC1');
qhc2 = wavread('Jay_ChouQHC2');
qhc3 = wavread('Jay_ChouQHC3');
qhc4 = wavread('Jay_ChouQHC4');
qhc5 = wavread('Jay_ChouQHC5');

dfp1 = wavread('Jay_Chou_DFP1');
dfp2 = wavread('Jay_Chou_DFP2');
dfp3 = wavread('Jay_Chou_DFP3');
dfp4 = wavread('Jay_Chou_DFP4');
dfp5 = wavread('Jay_Chou_DFP5');

frx1 = wavread('Jay_Chou_FRX1');
frx2 = wavread('Jay_Chou_FRX2');
frx3 = wavread('Jay_Chou_FRX3');
frx4 = wavread('Jay_Chou_FRX4');
frx5 = wavread('Jay_Chou_FRX5');

qhc1 = qhc1(:,1);
qhc1 = decimate(qhc1,scale);
qhc2 = qhc2(:,1);
qhc2 = decimate(qhc2,scale);
qhc3 = qhc3(:,1);
qhc3 = decimate(qhc3,scale);
qhc4 = qhc4(:,1);
qhc4 = decimate(qhc4,scale);
qhc5 = qhc5(:,1);
qhc5 = decimate(qhc5,scale);

dfp1 = dfp1(:,1);
dfp1 = decimate(dfp1,scale);
dfp2 = dfp2(:,1);
dfp2 = decimate(dfp2,scale);
dfp3 = dfp3(:,1);
dfp3 = decimate(dfp3,scale);
dfp4 = dfp4(:,1);
dfp4 = decimate(dfp4,scale);
dfp5 = dfp5(:,1);
dfp5 = decimate(dfp5,scale);

frx1 = frx1(:,1);
frx1 = decimate(frx1,scale);
frx2 = frx2(:,1);
frx2 = decimate(frx2,scale);
frx3 = frx3(:,1);
frx3 = decimate(frx3,scale);
frx4 = frx4(:,1);
frx4 = decimate(frx4,scale);
frx5 = frx5(:,1);
frx5 = decimate(frx5,scale);

qhc = [qhc1 qhc2 qhc3 qhc4 qhc5];
dfp = [dfp1 dfp2 dfp3 dfp4 dfp5];
frx = [frx1 frx2 frx3 frx4 frx5];
%% classical music
% Dmitri Shostakovich: Waltz No. 2,   Romance,    Piano Concerto No. 2: II. Andante
wn1 = wavread('DS-WN1');
wn2 = wavread('DS-WN2');
wn3 = wavread('DS-WN3');
wn4 = wavread('DS-WN4');
wn5 = wavread('DS-WN5');

and1 = wavread('DS-AND1');
and2 = wavread('DS-AND2');
and3 = wavread('DS-AND3');
and4 = wavread('DS-AND4');
and5 = wavread('DS-AND5');

ro1 = wavread('DS-RO1');
ro2 = wavread('DS-RO2');
ro3 = wavread('DS-RO3');
ro4 = wavread('DS-RO4');
ro5 = wavread('DS-RO5');

wn1 = wn1(:,1);
wn1 = decimate(wn1,scale);
wn2 = wn2(:,1);
wn2 = decimate(wn2,scale);
wn3 = wn3(:,1);
wn3 = decimate(wn3,scale);
wn4 = wn4(:,1);
wn4 = decimate(wn4,scale);
wn5 = wn5(:,1);
wn5 = decimate(wn5,scale);

and1 = and1(:,1);
and1 = decimate(and1,scale);
and2 = and2(:,1);
and2 = decimate(and2,scale);
and3 = and3(:,1);
and3 = decimate(and3,scale);
and4 = and4(:,1);
and4 = decimate(and4,scale);
and5 = and5(:,1);
and5 = decimate(and5,scale);

ro1 = ro1(:,1);
ro1 = decimate(ro1,scale);
ro2 = ro2(:,1);
ro2 = decimate(ro2,scale);
ro3 = ro3(:,1);
ro3 = decimate(ro3,scale);
ro4 = ro4(:,1);
ro4 = decimate(ro4,scale);
ro5 = ro5(:,1);
ro5 = decimate(ro5,scale);

wn = [wn1 wn2 wn3 wn4 wn5];
and = [and1 and2 and3 and4 and5];
ro = [ro1 ro2 ro3 ro4 ro5];

L = tr;
% all the same 2205
n = length(ro1);
t2 = linspace(0,L,n+1);t=t2(1:n);
% unit of pi is rad/s, divide by 2pi gives Hz, so use 1 instead of 2 pi
k1 = (1/L)*[0:n/2 - 1, -n/2:-1];
ks = fftshift(k1);
ks(1,end+1) = 0;




%Gabor filter
% spectrums of sweet child o mine
specSOM = [];
width = 7;                                                     
tslide =  0: 0.25 : 5;  
 for k = 1:5
    y = som(:,k).';
    Spec = [];
    for j = 1 : length(tslide)
         g = exp(-width*(t - tslide(j)).^2);
         % filter it out
         yf = g.*y;
         % Important !! fourier transform, the frequency content after the filter
         yft = ifftshift(fft(yf));
         Spec = [Spec; yft];
    end
     specSOM(:,:,k) = Spec;
 end
 
specNT = [];
width = 7;                                                     
tslide =  0: 0.25 : 5;  
 for k = 1:5
    y = nt(:,k).';
    Spec = [];
    for j = 1 : length(tslide)
         g = exp(-width*(t - tslide(j)).^2);
         % filter it out
         yf = g.*y;
         % Important !! fourier transform, the frequency content after the filter
         yft = ifftshift(fft(yf));
         Spec = [Spec; yft];
    end
     specNT(:,:,k) = Spec;
 end
 
specPC = [];
width = 7;                                                     
tslide =  0: 0.25 : 5;  
 for k = 1:5
    y = pc(:,k).';
    Spec = [];
    for j = 1 : length(tslide)
         g = exp(-width*(t - tslide(j)).^2);
         % filter it out
         yf = g.*y;
         % Important !! fourier transform, the frequency content after the filter
         yft = ifftshift(fft(yf));
         Spec = [Spec; yft];
    end
     specPC(:,:,k) = Spec;
 end
 
specQHC = [];
width = 7;                                                     
tslide =  0: 0.25 : 5;  
 for k = 1:5
    y = qhc(:,k).';
    Spec = [];
    for j = 1 : length(tslide)
         g = exp(-width*(t - tslide(j)).^2);
         % filter it out
         yf = g.*y;
         % Important !! fourier transform, the frequency content after the filter
         yft = ifftshift(fft(yf));
         Spec = [Spec; yft];
    end
     specQHC(:,:,k) = Spec;
 end
 
specFRX = [];
width = 7;                                                     
tslide =  0: 0.25 : 5;  
 for k = 1:5
    y = frx(:,k).';
    Spec = [];
    for j = 1 : length(tslide)
         g = exp(-width*(t - tslide(j)).^2);
         % filter it out
         yf = g.*y;
         % Important !! fourier transform, the frequency content after the filter
         yft = ifftshift(fft(yf));
         Spec = [Spec; yft];
    end
     specFRX(:,:,k) = Spec;
 end
 
specDFP = [];
width = 7;                                                     
tslide =  0: 0.25 : 5;  
 for k = 1:5
    y = dfp(:,k).';
    Spec = [];
    for j = 1 : length(tslide)
         g = exp(-width*(t - tslide(j)).^2);
         % filter it out
         yf = g.*y;
         % Important !! fourier transform, the frequency content after the filter
         yft = ifftshift(fft(yf));
         Spec = [Spec; yft];
    end
     specDFP(:,:,k) = Spec;
 end
 
specWN = [];
width = 7;                                                     
tslide =  0: 0.25 : 5;  
 for k = 1:5
    y = wn(:,k).';
    Spec = [];
    for j = 1 : length(tslide)
         g = exp(-width*(t - tslide(j)).^2);
         % filter it out
         yf = g.*y;
         % Important !! fourier transform, the frequency content after the filter
         yft = ifftshift(fft(yf));
         Spec = [Spec; yft];
    end
     specWN(:,:,k) = Spec;
 end
 
specRO = [];
width = 7;                                                     
tslide =  0: 0.25 : 5;  
 for k = 1:5
    y = ro(:,k).';
    Spec = [];
    for j = 1 : length(tslide)
         g = exp(-width*(t - tslide(j)).^2);
         % filter it out
         yf = g.*y;
         % Important !! fourier transform, the frequency content after the filter
         yft = ifftshift(fft(yf));
         Spec = [Spec; yft];
    end
     specRO(:,:,k) = Spec;
 end

specAND = [];
width = 7;                                                     
tslide =  0: 0.25 : 5;  
 for k = 1:5
    y = and(:,k).';
    Spec = [];
    for j = 1 : length(tslide)
         g = exp(-width*(t - tslide(j)).^2);
         % filter it out
         yf = g.*y;
         % Important !! fourier transform, the frequency content after the filter
         yft = ifftshift(fft(yf));
         Spec = [Spec; yft];
    end
     specAND(:,:,k) = Spec;
 end
 
 % Build SVD
 % order: SOM, NT, PC, QHC, DFP, FRX, WN, RO, AND
 X = [];
 for k = 1 : 5
     spec = specSOM(:,:,k);
     for j = 1:21
         temp = spec(j,:);
         X = [X; temp];
     end
 end
  for k = 1 : 5
     spec = specNT(:,:,k);
     for j = 1:21
         temp = spec(j,:);
         X = [X; temp];
     end
  end
  for k = 1 : 5
     spec = specPC(:,:,k);
     for j = 1:21
         temp = spec(j,:);
         X = [X; temp];
     end
  end
 
   for k = 1 : 5
     spec = specQHC(:,:,k);
     for j = 1:21
         temp = spec(j,:);
         X = [X; temp];
     end
 end
  for k = 1 : 5
     spec = specDFP(:,:,k);
     for j = 1:21
         temp = spec(j,:);
         X = [X; temp];
     end
  end
  for k = 1 : 5
     spec = specFRX(:,:,k);
     for j = 1:21
         temp = spec(j,:);
         X = [X; temp];
     end
  end
  for k = 1 : 5
     spec = specWN(:,:,k);
     for j = 1:21
         temp = spec(j,:);
         X = [X; temp];
     end
  end
  for k = 1 : 5
     spec = specRO(:,:,k);
     for j = 1:21
         temp = spec(j,:);
         X = [X; temp];
     end
  end
  for k = 1 : 5
     spec = specAND(:,:,k);
     for j = 1:21
         temp = spec(j,:);
         X = [X; temp];
     end
  end
 [u,s,v2] = svd(X,'econ');
 
%  figure(1)
% plot3(v2(1:105,2),v2(1:105,3),v2(1:105,4),'rs');hold on;
% plot3(v2(106:210,2),v2(106:210,3),v2(106:210,4),'ro');hold on;
% plot3(v2(211:315,2),v2(211:315,3),v2(211:315,4),'rd');hold on;
% plot3(v2(301:400,2),v2(301:400,3),v2(301:400,4),'b*');hold on;
% plot3(v2(401:600,2),v2(401:600,3),v2(401:600,4),'bo');hold on;
% plot3(v2(601:700,2),v2(601:700,3),v2(601:700,4),'m*');hold on;
% plot3(v2(701:800,2),v2(701:800,3),v2(701:800,4),'mo');hold on;
% plot3(v2(801:900,2),v2(801:900,3),v2(801:900,4),'md');
sig = diag(s);
energy = sig/sum(sig);

figure(4)
plot(energy,'o');
title('Energy of the system - Test 1')
 