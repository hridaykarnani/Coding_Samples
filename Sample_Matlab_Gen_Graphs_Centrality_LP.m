%% DGP - Simulation Diff Star Sup Networks
% Coded by: Hriday Karnani ---- P.I.: Jorge Miranda-Pinto
% Date: March 2023

% Introduction to this code and project: 
% This project:
% In this project we estimate unobserved networks with a new methodology
% that we name FASVAR (Factor-Augmented-Sparse-VAR). The idea is that
% unobserved networks are estimated at the same time as common factors in a
% reduced form estimator. To estimate the big amount of parameters we
% impose a bit of structure. Function FASVAR has the underlying code for
% this estimator, which is a non-linear convex optimization problem. We
% solve this optimization problem with CVX solver, which relies on
% numerical methods.

% This code:
% In this code we generate data to estimate our proposed estimator 
% following a Data Generating Process (DGP) and different degrees of 
% centrality on the network. The objective is to analyze how well the
% estimator performs under different degrees of centrality. We define a set
% of metrics to analyze the performance of the estimator and compare it to
% Lasso and factor model estimators.

clear all; clc; close all;
addpath('C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR')
addpath('C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\functions')

% Predefine variables that are going to be used in the loop
% Each iteration is a different degree of centrality. We start with a star
% supplier network and vary the amount of connections outside the first row
central=[];metrics=[];comm_fact=[];central_rmse=[];central2=[];
tic;
for j=1:50
    N=50; % nb of individuals (firms, sectors, countries)
    T=100; % nb of periods
    rng(6); % set seed

    f0=randn(T,1);%1 factor
    l0=randn(N,1);%loadings for factor
    Pi0=f0*l0';%factor component of error (note Pi0 has rank=1)
    U0=randn(T,N);%iid normal errors

    A0=rand(N); % The network
    A01 = A0(1,:);A02=A0(2:end,:); 
    p1=0.95;p2=0:0.005:0.25;
    A01(A01>p1)=0;A02(A02>p2(j))=0; % We separate connections on the first row and others
    % p1 and p2 are probab of having connections outside the first row
    A0 = vertcat(A01,A02); % re-join first row and others
    A0 = A0 - diag(diag(A0)); % set diagonals to 0 (no self-connections)
    A0(:,sum(A0,1)>0)=A0(:,sum(A0,1)>0)./(1.2.*(sum(A0(:,sum(A0,1)>0),1)));%normalize rows of A0 to sum to 1 (for all firms with at least one linkage) or zero (for firms with no linkages);

    shock = 0.15*Pi0(:,1)+0.85*U0; % has an aggregate and idiosync comp % Pi and U are the agg comp and residuals estimated by FASVAR

    Xn = shock*(eye(N) - A0)^(-1);  % simulate shock
    %Xn = randn(T,N)+ Pi0;
    X = Xn(1:end-1,:); % lag of simulated shock
    %X=Xn;
    Y = X*A0+ 0.1.*randn(T-1,N); % Simulated output based on shock and network
    shock = shock(1:end-1,:);%adapt shocks dimensions (consider lag)
    Pi0 = Pi0(1:end-1,:);%adapt dimensions (consider lag)
    %corshock = corr(shock);

    lambda2=1;
    central(:,j) = centrality(digraph(A0),'hubs');  % Centrality of the DGP network   
    %% FASVAR Estimation
    % Because of computational issues, we are going to cross-validate only
    % lambda 1
    
    % The estimation starts with a CV process, we set a grid of the
    % penalization parameters (lambda1) and select the parameter out of the
    % grid that gives the smallest MSE.
    
    % CV process
    cv=4; % Two cross-validations
    %tic;
    mse=[];
    for  i=1:cv 
        i;
        warning('off','all')
        indices = crossvalind('Kfold',X(:,1),cv) ; % Variables for CV
        test = (indices==i);train=~test;msel1l2=[];mseLl1=[];mseeLl1=[];
       parfor p=1:10
               %horzcat(i,p) % show in which iteration we are
               l1f=4.4:0.1:5.5; % Grid of lambda1 for FASVAR. We call FASVAR function
               % Note that we use a training subset and a test subset
               [A,Pi,U,a,u] = FASVAR(X(train,:),Y(train,:),l1f(p),lambda2,'no'); 
               yfit = X(test,:)*A ; % Prediction of Y based on the estimated netowrk (A) with test subset
               mse_i = mean(mean((yfit-Y(test,:)).^2)); % MSE of predicted Y and "real" Y
               msel1l2(p) = mse_i;
               l1l = 0.0:0.005:0.05;
            [AL,UL,aL,uL] = Lasso(Y(train,:),X(train,:),l1l(p)); % Same idea for Lasso estimator
            yfitL = X(test,:)*AL ;
            mseL_i = mean(mean((yfitL-Y(test,:)).^2));
            mseLl1(p) = mseL_i;
            [AeL,UeL,aeL,ueL] = Elastic_Net_Lasso(Y(train,:),X(train,:),l1l(p),0.5); % and Elastic Net Lasso
            yfiteL = X(test,:)*AeL ;
            mseeL_i = mean(mean((yfiteL-Y(test,:)).^2));
            mseeLl1(p) = mseeL_i;        
       end
       mse(:,i) = msel1l2;
       mseL(:,i) = mseLl1;
       mseeL(:,i) = mseeLl1;
    end
    %toc;

    % Select lambda1 with least MSE for FASVAR, Lasso and Elastic Net Lasso
    mmse = mean(mse,2);
    [M,I] = (min(mmse)); 
    mmseL = mean(mseL,2);
    [ML,IL] = (min(mmseL));
    mmseeL = mean(mseeL,2);
    [MeL,IeL] = (min(mmseeL));

    l1f=4.4:0.1:5.5;l1l= 0.0:0.005:0.05;
    lambda1F = l1f(I);
    lambda1L = l1l(IL);lambda1ENL=l1l(IeL);

    %% Estimation with CV param.
    % Estimate FASVAR, Lasso and EN Lasso with selected lambda1 and compare
    % the performance
    rng(6);
    test= rand(T-1,1)<.3; 
    train= ~test;

    [A,Pi,~,~,u] = FASVAR(X(train,:),Y(train,:),lambda1F,lambda2,'no');
    [AL,~,~,uL] = Lasso(Y(train,:),X(train,:),lambda1L);
    [AeL,~,~,ueL] = Elastic_Net_Lasso(Y(train,:),X(train,:),lambda1ENL,0.5);
    
    A_spa=[];MSE=[];mean_resid=[];mean_dev=[];r2aux=[];r2aux2=[];r2aux3=[];
    % check how sparse each matrix is:
    A_spa(:,1)=nnz(A)/(N*(N-1));A_spa(:,2)=nnz(AL)/(N*(N-1));A_spa(:,3)=nnz(AeL)/(N*(N-1));
    %A_spa
    % Out of sample prediction
    % Should add Pi to forecast of FASVAR, can't do it though.
    % A useful measure to consider this is in-sample residuals, which considers
    % Pi
    Yf = X(test,:)*A;
    Yf_L = X(test,:)*AL;
    Yf_eL = X(test,:)*AeL;

    % MSE
    MSE(:,1) = mean(mean((Yf-Y(test,:)).^2));
    MSE(:,2)= mean(mean((Yf_L-Y(test,:)).^2));
    MSE(:,3)= mean(mean((Yf_eL-Y(test,:)).^2));
    %MSE

    % In sample mean residuals:
    mean_resid(:,1) = mean(u);mean_resid(:,2) = mean(uL);mean_resid(:,3) = mean(ueL);
    %mean_resid

    % Mean deviation on network estimation
    mean_dev(:,1) = mean(mean((A-A0).^2));mean_dev(:,2) = mean(mean((AL-A0).^2));mean_dev(:,3) = mean(mean((AeL-A0).^2));
    %mean_dev
    % Centrality of each network (and the network by DGP)
    centralityIOS = centrality(digraph(A0),'hubs')*100;
    centralityA = centrality(digraph(A),'hubs')*100;
    centralityAL = centrality(digraph(AL),'hubs')*100;
    centralityAeL = centrality(digraph(AeL),'hubs')*100;

    central2(:,:,j) = [centralityIOS centralityA centralityAL centralityAeL];
    % RMSE of centrality
    central_rmse = horzcat(mean((centralityIOS-centralityA).^2).^0.5,mean((centralityIOS-centralityAL).^2).^0.5,mean((centralityIOS-centralityAeL).^2).^0.5);
    % Variable with the comparison of the three methodologies that has all
    % the previous metrics
    metrics(:,:,j) = vertcat(A_spa,repelem(nnz(A0)./(N*(N-1)),3),MSE,mean_resid,mean_dev,central_rmse);


    %% Variability Explained by Common Factor
    [A,Pi,U,a,u] = FASVAR(X,Y,lambda1F,lambda2,'no'); % FASVAR with full sample
    Pi_svs=svds(Pi,5);%notice that most singular values are zero because of the nuclear norm penalty (for large enough lambda2)

    [~,score,~,~,explained] = pca(zscore(Y)); % Compare FASVAR to a factor model
    varfactor = var(score(:,1));
    shar = Y./sum(Y,2); 
    w = mean(shar,1)';

    % check the agg. variability explained by common factors by factor model
    g_w = (Y*w);
    varg_w = var(g_w);
    G = score(:,1);
    bhat_agg    =   (inv(G'*G)*(G'*Y))'; % estimating loadings by OLS
    var_factorsagg2 = w'*bhat_agg*varfactor*bhat_agg'*w;
    R_F_agg2 = (var_factorsagg2/varg_w)*100;

    % check the agg. variability explained by common factors by FASVAR
    varfactor_fasvar = var(Pi(:,1));
    bhat_agg_fasvar = (inv(Pi(:,1)'*Pi(:,1))*(Pi(:,1)'*Y))';
    var_factorsagg_fasvar = w'*bhat_agg_fasvar*varfactor_fasvar*bhat_agg_fasvar'*w;
    R_F_agg_fasvar = (var_factorsagg_fasvar/varg_w)*100;

    varfactor_fasvar2 = cov(Pi);
    bhat_agg_fasvar2 = (inv(Pi'*Pi)*(Pi'*Y))';
    var_factorsagg_fasvar2 = w'*bhat_agg_fasvar2*varfactor_fasvar2*bhat_agg_fasvar2'*w;
    R_F_agg_fasvar2 = (var_factorsagg_fasvar2./varg_w)*100;

    % check the sectoral variability explained by common factors by factor
    % model and FASVAR
    for i=1:N
        r2aux(:,i) = fitlm(Pi0(:,1),Y(:,i),'intercept',false).Rsquared.Ordinary;
        r2aux2(:,i) = fitlm(Pi(:,1),Y(:,i),'intercept',false).Rsquared.Ordinary;
        r2aux3(:,i) = fitlm(score(:,1),Y(:,i),'intercept',false).Rsquared.Ordinary;
    end
    mean(r2aux);mean(r2aux2);mean(r2aux3);
    vertcat(mean((r2aux3-r2aux).^2),mean((r2aux2-r2aux).^2));
    % check the agg variab explained by the DGP
    varfactor_fasvar0 = var(Pi0(:,1));
    bhat_agg_fasvar0 = (inv(Pi0(:,1)'*Pi0(:,1))*(Pi0(:,1)'*Y))';
    var_factorsagg_fasvar0 = w'*bhat_agg_fasvar0*varfactor_fasvar0*bhat_agg_fasvar0'*w;
    R_F_agg_fasvar0 = (var_factorsagg_fasvar0/varg_w)*100;

    % compare agg variab explained by DGP, Factor model and FASVAR
    comm_fact(:,j) = vertcat(R_F_agg2, R_F_agg_fasvar, R_F_agg_fasvar0,mean((r2aux3-r2aux).^2),mean((r2aux2-r2aux).^2));
    %toc;
    j
    % show which iteration we are in
end
toc;

%% Plots with results

% save results
save('C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\sim_starsupcv_norm_crmse.mat',...
    'central','comm_fact','metrics','central_rmse','central2')
load('C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\simstarsupcv_crmse.mat')

plot(central(1,:),(comm_fact(4:5,:))')
title('MSE of Variab. explained by Common Factors','interpreter','latex')
legend('PCA Factor','FASVAR 1st Factor');xlabel('Centrality 1st Firm');ylabel('MSE')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\crmse\rmse_cf_cen','epsc')

plot(central(1,:),(abs(comm_fact(1:2,:)-comm_fact(3,:)))')
title('MSE of Variab. explained by Common Factors','interpreter','latex')
legend('PCA Factor','FASVAR 1st Factor');xlabel('Centrality 1st Firm');ylabel('MSE')

nnze = reshape(abs(metrics(1,:,:)-metrics(2,:,:)),3,50)';
plot(central(1,:),nnze)
title('Non-sparse Elements Absolute Error','interpreter','latex')
legend('FASVAR','Lasso','Elastic Net Lasso');xlabel('Centrality 1st Firm');ylabel('Abs. Error')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\crmse\nnz_abserr_cen','epsc')

mse = reshape(metrics(3,:,:),3,50)';
plot(central(1,:),mse)
title('Out of Sample MSE','interpreter','latex')
legend('FASVAR','Lasso','Elastic Net Lasso');xlabel('Centrality 1st Firm');ylabel('MSE')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\crmse\ofs_mse_cen','epsc')

mresid = reshape(abs(metrics(4,:,:)),3,50)';
plot(central(1,:),mresid)
title('In-sample Mean Residual','interpreter','latex')
legend('FASVAR','Lasso','Elastic Net Lasso');xlabel('Centrality 1st Firm');ylabel('Abs. Value Mean Residual')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\crmse\ins_meanresid_cen','epsc')

mseA = reshape(metrics(5,:,:),3,50)';
plot(central(1,:),mseA)
title('MSE from Actual $\Gamma$ matrix','interpreter','latex')
legend('FASVAR','Lasso','Elastic Net Lasso');xlabel('Centrality 1st Firm');ylabel('MSE')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\crmse\mseA_cen','epsc')

mseCent = reshape(metrics(6,:,:),3,50)';
plot(central(1,:),mseCent)
title('Centrality MSE','interpreter','latex')
legend('FASVAR','Lasso','Elastic Net Lasso');xlabel('Centrality 1st Firm');ylabel('MSE')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\crmse\mseCent_cen','epsc')

mseCent = reshape(metrics(6,:,:),3,50)';
plot(1:50,mseCent)
title('Centrality MSE','interpreter','latex')
legend('FASVAR','Lasso','Elastic Net Lasso');xlabel('Pagerank Centrality 1st Firm');ylabel('Centrality MSE')

%% FM and VAR-Lasso Plots

plot(central(1,:),(comm_fact(4,:))')
title('MSE of Variab. explained by Common Factors','interpreter','latex')
legend('PCA Factor');xlabel('Pagerank Centrality 1st Firm');ylabel('MSE')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\colnorm_cv\rmse_cf_fm_cen','epsc')

nnze = reshape(abs(metrics(1,2:3,:)-metrics(2,2:3,:)),2,50)';
plot(central(1,:),nnze)
title('Non-sparse Elements Absolute Error','interpreter','latex')
legend('Lasso','Elastic Net Lasso');xlabel('Pagerank Centrality 1st Firm');ylabel('Abs. Error')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\colnorm_cv\nnz_abserr_vlasso_cen','epsc')

mse = reshape(metrics(3,2:3,:),2,50)';
plot(central(1,:),mse)
title('Out of Sample MSE','interpreter','latex')
legend('Lasso','Elastic Net Lasso');xlabel('Pagerank Centrality 1st Firm');ylabel('MSE')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\colnorm_cv\ofs_mse_vlasso_cen','epsc')

mresid = reshape(abs(metrics(4,2:3,:)),2,50)';
plot(central(1,:),mresid)
title('In-sample Mean Residual','interpreter','latex')
legend('Lasso','Elastic Net Lasso');xlabel('Pagerank Centrality 1st Firm');ylabel('Abs. Value Mean Residual')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\colnorm_cv\ins_meanresid_vlasso_cen','epsc')

mseA = reshape(metrics(5,2:3,:),2,50)';
plot(central(1,:),mseA)
title('MSE from Actual $\Gamma$ matrix','interpreter','latex')
legend('Lasso','Elastic Net Lasso');xlabel('Pagerank Centrality 1st Firm');ylabel('MSE')
saveas(gcf,'C:\Users\Hriday\Dropbox\Hriday\Factor_Models\Replication_Prev_Methodologies\FASVAR\results\sim_starsup\plots\colnorm_cv\mseA_vlasso_cen','epsc')

