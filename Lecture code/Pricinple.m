clear all; close all;clc
load catData.mat
load dogData.mat

C = cat;
D = dog;
for j = 1:9
    subplot(3,3,j)
    imshow(reshape(D(:,j),64,64))
end

X = [double(C) double(D)];
[u,s,v] = svd(X,'econ');
figure(2),plot(diag(s)/sum(diag(s)),'ko','Linewidth',[2])
figure(3)
for j = 1:4
    subplot(2,2,j)
    mode = reshape(u(:,j),64,64);
    pcolor(flipud(mode)),shading interp, colormap(hot)
    
end
figure(4)
for j = 1:3
    subplot(3,1,j)
    bar(v(80+j,:))
end
break
figure(5)
plot3(v(2,1:80),v(3,1:80),v(4,1:80),'ro','Linewidth',[2]), grid on
hold on
plot3(v(2,81:end),v(3,81:end),v(4,81:end),'bo','Linewidth',[2]), grid on


