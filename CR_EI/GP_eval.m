function [ybar,vary] = GP_eval(fit,x_test)

%Normalization
mx  = fit.nparx(1,:); sx   = fit.nparx(2,:); 
my  = fit.npary(1,:); sy   = fit.npary(2,:);


mg=size(x_test,1);
x_test_Norm =(x_test-repmat(mx,mg,1))./repmat(sx,mg,1);

%Prediction
[ybar_Norm vary_Norm] = gp(fit.hyp, @infGaussLik, [], fit.covfunc, fit.likfunc, fit.X_Norm, fit.Y_Norm, x_test_Norm);

ybar = ybar_Norm*sy + my;
vary = vary_Norm*(sy*sy);
