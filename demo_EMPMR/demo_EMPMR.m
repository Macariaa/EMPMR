
clc
clear all
close all

load p(Bun)

ratio=  0.01;                                              % Parameter represents the ratio of outliers
for i=1:N
   % model{i} = awgn(model{i},25);                         % Add noises (25 dB)
    ns{i}= createns(model{i});                             % Build k-d tree
    [corr,TD] = knnsearch(ns{i},model{i},'k',2);           % The search of nearest neighbor
    sigma(i)= (mean(TD(:,2)));                             % Calculate the point resolution
end
sigma= mean(sigma);
sigma= (sigma).^2;                                         % Initial guess of the isotropic covariance
for i=1:N                                                  % Transform poin sets from a set-centered frame into the model centered frame
    R0= p(i).M(1:3,1:3);    t0= p(i).M(1:3,4);
    model{i}= transform_to_global(model{i}', R0, t0);
end

%EMMR algorithm                                        
tic
[pM, Tmodel]= EMPMR(p, model, ns, sigma, ratio, N);        %  EMPMR alrogrithm 
[errR,errT]= err_cal(GrtR, GrtT, pM, N, s);               %  Compute the registration error.
toc


% Display the registration results
Model=obtain_model(shape, pM, N, s);                       % Obtain the multi-view registration shape   
crosssection(Model,1,-18.7,-18.4);                         % Obtain cross-section of registration results
axis off
axis equal









