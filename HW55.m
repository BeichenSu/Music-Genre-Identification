clear all; close all; clc;
genre = cell(1,9);
genre{1,1} = 'GNRWJ.mp3';
genre{1,2} = 'GNRNT.mp3';
genre{1,3} = 'GNRPC.mp3';
% genre{1,4} = 'Jay Chou QHC.mp3';
% genre{1,5} = 'Jay Chou JHT.mp3';
genre{1,6} = 'Jay Chou FRX.mp3';
genre{1,4} = 'NV_SS.mp3';
genre{1,5} = 'NV_UK.mp3';
%genre{1,6} = 'NV_LI.mp3';
% genre{1,7} = 'DS-AND.mp3';
% genre{1,8} = 'DS-RO.mp3';
genre{1,9} = 'DS-WN.mp3';
genre{1,7} = 'AIC MB.mp3';
genre{1,8} = 'AIC DF.mp3';
%genre{1,9} = 'AIC W.mp3';
% 5 second
tr = 5;
songSpec = [];
for j = 1:9
    currentSong = audioread(genre{1,j});
    currentSong = currentSong(:,1);
    %change sampling rate to pick 5-second clips and rescale
    currentSong = decimate(currentSong,100);
    currentSong = (currentSong - mean(currentSong))/var(currentSong);
    % 5 second take 2205 position in the song, start at 10000 position of
    % the song
    clips = 30;
    position = 10000;
    L = 5; % 5 second of clips
    tslide =  0: 0.25 : 5; 
    width = 6;
    allClipSpec = [];
    for k = 1:30
        currentClip = currentSong(position + 1:position + 2205,1);
        position = position + 2205;
        clipSpec = [];
        for p = 1 : length(tslide)
            % Gaussian filter
            % time span
            Fs = length(currentClip)/L;
            t = (1:length(currentClip))/Fs;
            g = exp(-width*(t - tslide(p)).^2);
            % filter it out
            currentClipf = g.*currentClip';
    
            % Important !! fourier transform, the frequency content after the filter
            currentClipft = abs(ifftshift(fft(currentClipf)));
            clipSpec = [clipSpec currentClipft'];
        end
        allClipSpec = [allClipSpec clipSpec];
    end
    songSpec = [songSpec allClipSpec];
end
[u,s,v] = svd(songSpec,'econ');
sig = diag(s);
energy = sig/sum(sig);
figure(1)
plot(energy,'o');
title('Energy of the system - Test 3')
% figure(2)
% plot3(v(1:630,2),v(1:630,3),v(1:630,4),'rs');hold on;
% plot3(v(631:1260,2),v(631:1260,3),v(631:1260,4),'rd');hold on;
% plot3(v(1261:1890,2),v(1261:1890,3),v(1261:1890,4),'ro');hold on;
% plot3(v(1891:2520,2),v(1891:2520,3),v(1891:2520,4),'bo');hold on;
% plot3(v(2521:3150,2),v(2521:3150,3),v(2521:3150,4),'bs');hold on;
% plot3(v(3151:3780,2),v(3151:3780,3),v(3151:3780,4),'bd');hold on;
% plot3(v(3781:4410,2),v(3781:4410,3),v(3781:4410,4),'md');hold on;
% plot3(v(4411:5040,2),v(4411:5040,3),v(4411:5040,4),'mo');hold on;
% plot3(v(5041:5670,2),v(5041:5670,3),v(5041:5670,4),'ms');hold on;
figure(2)
plot3(v(1:1890,2),v(1:1890,3),v(1:1890,4),'ro');hold on;
plot3(v(1891:3780,2),v(1891:3780,3),v(1891:3780,4),'gd');hold on;
plot3(v(3781:5670,2),v(3781:5670,3),v(3781:5670,4),'bs');hold on;
title('Projection of v matrix of 3 different genres');
legend('Rock','Pop','Classic')
% 
q1=randperm(1890);
q2=randperm(1890);
q3=randperm(1890);
rock = v(1:1890,2:4);
pop = v(1891:3780,2:4);
classic = v(3781:end,2:4);
xtrain=[rock(q1(1:1600),:);pop(q2(1:1600),:);classic(q3(1:1600),:)];
ctrain=[ones(1600,1);2*ones(1600,1);3*ones(1600,1)];
xtest=[rock(q1(1601:end),:);pop(q2(1601:end),:);classic(q3(1601:end),:)];



nb=fitNaiveBayes(xtrain,ctrain);
pre1=nb.predict(xtest);

gm = fitgmdist(xtrain,3);
pre2 = cluster(gm,xtest);


figure(3)
bar(pre1)
title('NaiveBayes')

figure(4)
bar(pre2)
title('Gaussian mixed')

% figure(4)
% plot3(v(1:630,2),v(1:630,3),v(1:630,4),'bs');hold on;
% plot3(v(631:1260,2),v(631:1260,3),v(631:1260,4),'rs');hold on;
% plot3(v(1261:1890,2),v(1261:1890,3),v(1261:1890,4),'gs');hold on;


