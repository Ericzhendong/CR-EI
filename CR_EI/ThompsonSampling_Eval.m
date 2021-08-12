function [gx] = ThompsonSampling_Eval(x, objOPTS)

    TS_FIT = objOPTS.FIT;
    TS_OPT = objOPTS.OPT;
    
    WW_dim = TS_FIT.WW_dim;
    WW = TS_FIT.WW;
    bias = TS_FIT.bias;
    mean_theta_TS = TS_FIT.mean_theta_TS;
    
    %% normalization of x
    
    mx = TS_OPT.nparx(1,:); 
    sx = TS_OPT.nparx(2,:); 

    mg=size(x,1);
    x =(x-repmat(mx,mg,1))./repmat(sx,mg,1);

    %% Calculate the thompson sample values at the point x
    phi_x = sqrt(2.0/WW_dim)*[sin(x*WW + bias), cos(x*WW + bias)];
    gx = phi_x * mean_theta_TS;
    
end