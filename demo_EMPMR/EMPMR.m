function [pM1, model]= EMPMR(pM, model, ns, sigma, ratio, N)
%     EMPMR    Expectation-maximization per-spective for multi-view registration.
%     Input:
%        'pM'     Initial rigid transformations, which are 4*4 matrix;
%        'model'  Coarse-aligned shape
%        'ns'     K-d tree for  E-Corresponding-step
%        'sigma'  Initial isotropic covariance
%        'ratio'  The ratio of outliers
%        'N'      The number of point sets
%      Output:
%       ¡®pM1¡¯   Obtained rigid transformations
%        'model'  Well-aligned shape
%      References:
%        [1] Jihua Zhu, Jing Zhang, Jing Zhang, Shanmin Pang, Huimin Lu, and Zhongyu Li,
%            Registration of multi-view point sets under theperspective of expectation-maximization, Submitted to IEEE TIP, 2020.

iter= 0;
lambda= ((N-1)/(N))*ratio/(1-ratio);
E_cov= 0.1;
Cov= 10;

while ((iter<200)&(E_cov>(5*10^(-4))))
    iter= iter+1;
    TDss=[];
    wss=[];
    for i=1:N                                                         % Sequentially align each point set
        Datas= [];
        Mscans=[];
        ws= [];   
        TDs= [];
        % E_Step
        Data= model{i};
        for j=1:N                                                     % Utilize all availble information from other data sets
            if(j~=i)
                Mj= inv(pM(j).M);
                TData= transform_to_global(Data, Mj(1:3,1:3), Mj(1:3,4));
                [corr,TD] = knnsearch(ns{j},TData');                  % E-Corresponding-step
                Datas=[Datas,Data];
                Mscans=[Mscans,model{j}(:,corr)];
                TD2= TD.^2;
                wj= exp(-TD2/(2*sigma))/(2*pi*sigma)^1.5;             % E-Probability-Step
                ws= [ws,wj];
                TDs= [TDs,TD2];
            end
        end
        sw= sum(ws')+lambda;                                           % E-Probability-Step
        sw= repmat(sw',1,N-1);
        ws= ws./sw;                                                    % E-Probability-Step (normalization)
        TDss=[TDss;TDs];
        wss=[wss;ws];
        [rn,cn]= size(ws);
        w= reshape(ws,rn*cn,1);
        
        % M_step
        M= reg(Mscans,Datas,w);                                         % M-step (update the ith rigid transformation)
        pM= recover(pM,M,i);
        model{i}= transform_to_global(model{i}, M(1:3,1:3), M(1:3,4));
    end
   PCov= Cov;
   Cov= -sum(sum(wss.*(TDss./sigma+3*log(sigma))))'/size(TDss,1);       % The log-likelihood value ignoring constant terms
   E_cov= abs(PCov-Cov);
   sigma= sum(sum(wss.*TDss)')./(3*sum(sum(wss)'));                     % M-step (update the isotropic covariance)
end


% Attach the model centered frame to the 1st set-centered frame
for i= 1:N
    pM1(i).M= inv(pM(1).M)*pM(i).M;
end

% Update the rigid transfromaton by the optimzation of weighted least square problem
function M = reg(M,S,p)
sumP= sum(p);
p1= repmat(p',3,1);
mm= sum((M.*p1)')./sumP;
ms= sum((S.*p1)')./sumP;
Sshifted = [S(1,:)-ms(1); S(2,:)-ms(2); S(3,:)-ms(3)];
Mshifted = [M(1,:)-mm(1); M(2,:)-mm(2); M(3,:)-mm(3)];
K= Sshifted.*p1*Mshifted';
[U A V] = svd(K);
R1 = V*U';
if det(R1)<0
    B = eye(3);
    B(3,3) = det(V*U');
    R1 = V*B*U';
end
t1= sum(((M-R1*S).*p1)')./sumP;
t1= t1';
M= Rt2M(R1,t1);              




