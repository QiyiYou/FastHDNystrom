% clc;
clear;
close all;
%% input
% [filename, user_canceled] = imgetfile; %% Uncomment this line to choose any image from your system
filename='peppers256.jpg';
I  =  double(imread(filename));
[m,n,d]=size(I);
sigs=5;
sigr=50;
Iact=I./255;
fast_flag=1;
%% direct bilateral filtering
tic,
Idirectbf=bilateral_filter(I,sigs,sigr);    %exact bilateral
Tdirect=toc;
fprintf('Direct bilateral filter : \n');
fprintf('time for direct bilateral(ms)=%3.0f \n',Tdirect*1000);
%% Kmeans filtering
% Done in two steps : Clustering and Filtering
%% Filtering
Cluster=8;
%% ICIP algorithm with change in clustering, weights for approximation for range kernel and interpolation between clusters
tic,
%Bisecting K-means clustering
% Ares=reshape(Iact(1:4:end,1:4:end),m*n/(4*4),d);
% Centre=kmeans_recursive(Ares,Cluster);
%Inbuilt Matlab clustering code     
[~,Centre] = rgb2ind(uint8(I(1:4:end,1:4:end,:)),Cluster,'nodither');

spatialtype='gaussian';     
convmethod='O1'; % Change convmethod to 'O1' for O(1) convolutions
Ikmean=fastKmeansfiltapproxnystrom(Iact,sigs,sigr/255,Centre,spatialtype,convmethod,fast_flag);      % bilateral kmeans
Ikmean=Ikmean.*255;
Ikmean(Ikmean<0)=0;
Ikmean(Ikmean>255)=255;
Tkmeans=toc;
fprintf('Fast bilateral filter by Kmeans complete with %d clusters \n',size(Centre,1));
fprintf('time for fast bilateral(ms)=%3.0f \n',Tkmeans*1000);
% 
%% Detail enhancement
% lamda=1.2;
% alpha=0.6;
% Il=Iact+(lamda.*((Iact-Ikmean).^alpha));
% Il=Il*255;
%figure;
%imshow(uint8(Il));%title('Enhanced image');
%% mse
error2 = reshape(Idirectbf-Ikmean, [d*m*n,1]);
MSE_mcbf2 = sqrt(sum(error2.^2)/(d*m*n));
PSNR2=20*log10(255/(MSE_mcbf2));
fprintf('mean sq error=%f, PSNR = %f db  \n',MSE_mcbf2,PSNR2);

% % %% output
figure;
imshow(uint8(I));%title('Original image');
figure;
imshow(uint8(Idirectbf));%title('direct bilateral filter');
figure;
imshow(uint8(Ikmean));%title('fast bilateral filter');
