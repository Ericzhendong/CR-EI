function [yE3I] = E3I(x, y_stars, GP_FIT)

    [mean_x,var_x] = GP_eval(GP_FIT,x);
    var2 = var_x;
%     var2 = max(var_x, repmat(1.0e-8,size(var_x,1),1));
    
    n_stars = size(y_stars,1);
    
    yE3I = zeros(size(x,1),1);
    
    for i = 1:n_stars
%         z = (mean_x - y_stars(i))./sqrt(var2);
        z = (y_stars(i)-mean_x)./sqrt(var2);
        out = ( y_stars(i)-mean_x).*normcdf(z, 0, 1) + sqrt(var2).*normpdf(z, 0, 1);
        
        yE3I = yE3I + out;  
    end
end


% u  = (yPBS - yhat)./sigma;
% EI = (yPBS - yhat).*normcdf(u, 0, 1) + sigma.*normpdf(u, 0, 1);