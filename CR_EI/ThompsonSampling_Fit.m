function [TS_OPT, TS_FIT] = ThompsonSampling_Fit(GP_FIT)

%     hyp_cov = GP_FIT.hyp.cov(1);
    
    hyp_cov = GP_FIT.hyp.cov(1:end-1);
    hyp_lik = GP_FIT.hyp.lik;

%     gp_x = GP_OPT.X;
%     gp_y = GP_OPT.Y;

    gp_x = GP_FIT.X_Norm;
    gp_y = GP_FIT.Y_Norm;
    
    dim = size(gp_x,2);
    GP_lengthscale = exp(hyp_cov);
    GP_noise_delta = exp(hyp_lik);

    % Parameters for Thompson Sampling
    WW_dim = 200;    % dimension of random features
    
    WW = mvnrnd(zeros(1, WW_dim),eye(WW_dim), dim);
    for i= 1:dim
        WW(i, :) = WW(i, :)/GP_lengthscale(i);
    end
    
    % WW = mvnrnd(zeros(1, WW_dim),eye(WW_dim), dim)/GP_lengthscale;
    bias = unifrnd(0, 2*pi(),1, WW_dim);

    Phi_X = sqrt(2.0/WW_dim)*[sin(gp_x*WW + bias), cos(gp_x*WW + bias)];

    A = (Phi_X')*Phi_X + eye(2*WW_dim)*GP_noise_delta;
    gx = (Phi_X')*gp_y;

    mean_theta_TS = A\gx;
    
    TS_FIT.WW_dim = WW_dim;
    TS_FIT.WW = WW;
    TS_FIT.bias = bias;
    TS_FIT.mean_theta_TS = mean_theta_TS;
    
    TS_OPT.X = gp_x;
    TS_OPT.Y = gp_y;
    
    TS_OPT.nparx  = GP_FIT.nparx;
    TS_OPT.npary  = GP_FIT.npary;
    

    
end