%% =================================================================
% This script performs low-fibered-rank-based HSI mixed noise removal models.  
% listed as follows:
%     1. 3DTNN-based HSI mixed noise removal model
%     2. 3DLogTNN-based HSI mixed noise removal model
%
% More detail can be found in [1]
% [1] Yu-Bang Zheng, Ting-Zhu Huang*, Xi-Le Zhao*, Tai-Xiang Jiang, Tian-Hui Ma, and Teng-Yu Ji.
%     Mixed Noise Removal in Hyperspectral Image via Low-Fibered-Rank Regularization.
%
%
% Created by Yu-Bang Zheng £¨zhengyubang@163.com£©
% 8/20/2019

clear;
clc;
close all;

addpath(genpath('lib'));
addpath(genpath('data'));
EN_3DTNN    = 1;
EN_3DLogTNN = 1;
methodname  = { 'Noise', '3DTNN', '3DLogTNN'};
Mnum = length(methodname);


%%

load('cleanPavia.mat')
Ohsi=img_clean;
if max(Ohsi(:))>1
    Ohsi=my_normalized(Ohsi);
end

Nway = size(Ohsi);

%% evaluation indexes
Re_hsi  =  cell(Mnum,1);
psnr    =  zeros(Mnum,1);
ssim    =  zeros(Mnum,1);
sam     =  zeros(Mnum,1);
time    =  zeros(Mnum,1);
%%  corrupted image


sigma_n3=0.1*rand(Nway(3),1)+0.05;  % Gaussian noise
sigma= mean(sigma_n3);

p_n3=0.2*rand(Nway(3),1)+0.1;  % salt and pepper noise
p= mean(p_n3);
    
  

fprintf('=== The Gaussian noise level is  %4.3f ===\n', sigma);
fprintf('=== The impulsive noise level is %4.3f ===\n', p);
Nhsi=zeros(Nway);
for j=1:Nway(3)
    Nhsi(:,:,j) = imnoise(Ohsi(:,:,j),'salt & pepper',p_n3(j))+sigma_n3(j)*randn(Nway(1),Nway(2));
end

i  = 1;
Re_hsi{i} = Nhsi;
[psnr(i), ssim(i), sam(i)] = MSIQA(Ohsi * 255, Re_hsi{i} * 255);
enList = 1;

%% Performing 3DTNN
i = i+1;
if EN_3DTNN
    %%%%%
    opts=[];
    opts.sigma = sigma;
    opts.theta = 0.001; % controls \alpha     weight of fibered-rank
    opts.phi   = 0.004; % controls \lambda_1  \|N\|_F^2
    opts.varpi = 1;     % controls \lambda_2  \|S\|_1
    opts.omega = 100;   % controls \tau       SVT
    %opts.Xtrue = Ohsi;
    %%%%%
    fprintf('\n');
    disp(['performing ',methodname{i}, ' ... ']);
    t0= tic;
    [Re_hsi{i},~,~,~]=de_3DTNN(Nhsi,opts);
    time(i) = toc(t0);
    [psnr(i), ssim(i), sam(i)] = MSIQA(Ohsi * 255, Re_hsi{i} * 255);
    enList = [enList,i];
    fprintf('3DTNN: PSNR=%5.4f   \n',  psnr(i));
end
%% Performing 3DLogTNN
i = i+1;
if EN_3DLogTNN
    %%%%%
    opts=[];
    opts.sigma  = sigma;
    opts.theta  = 0.001;   % controls \alpha     weight of fibered-rank
    opts.phi    = 0.00005; % controls \lambda_1  \|N\|_F^2
    opts.varpi  = 0.011;   % controls \lambda_2  \|S\|_1
    opts.omega  = 10000;   % controls \tau       SVT
    opts.logtol = 80;      % \varepsilon in 3DLogTNN
    %opts.Xtrue = Ohsi;
    
    %%%%%
    fprintf('\n');
    disp(['performing ',methodname{i}, ' ... ']);
    t0= tic;
    [Re_hsi{i},~,~,~]=de_3DLogTNN(Nhsi,opts);
    time(i) = toc(t0);
    [psnr(i), ssim(i), sam(i)] = MSIQA(Ohsi * 255, Re_hsi{i} * 255);
    enList = [enList,i];
    fprintf('3DLogTNN: PSNR=%5.4f   \n',  psnr(i));
end


%% Show result

fprintf('\n');
fprintf('================== Result ==================\n');
fprintf('      %5.2s %5.3f      %5.2s %5.3f    \n', 'G:',sigma, 'S:', p);
fprintf('================== Result ==================\n');
fprintf(' %8.8s    %5.4s      %5.4s    %5.4s    \n', 'method','PSNR', 'SSIM', 'SAM');
for i = 1:length(enList)
    fprintf(' %8.8s    %5.4f    %5.4f    %5.4f    \n',...
        methodname{enList(i)},psnr(enList(i)), ssim(enList(i)), sam(enList(i)));
end
fprintf('================== Result ==================\n');

%%
figure,
showHSIResult(Re_hsi,Ohsi,0,1,methodname,enList,1,Nway(3))



