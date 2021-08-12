function GEI = Generalized_ExpectedImprovement(yPBS, yhat, sigma,gIndex)

u  = (yPBS - yhat)./sigma;
T(:,1) = normcdf(u,0,1);   %T0
T(:,2) = -normpdf(u,0,1);  %T1
%

if gIndex == 0
    GEI = T(:,1);   %T0
else
    for i = 3:gIndex+1
        k = i-1;
        T(:,i) = -(u.^(k-1)).*normpdf(u,0,1) + (k-1)*T(:, i-2);
    end
    
    GEI = 0;
    for i = 1:gIndex+1
        k    = i-1;
        ng   = factorial(gIndex)/(factorial(k)*factorial(gIndex-k));
        dtmp = ((-1)^(k))*ng*(u.^(gIndex-k)).*T(:,i);
        GEI  = GEI + dtmp;
    end
    GEI = GEI.*(sigma.^gIndex);
end


% GEI = GEI.^(1.0/gIndex);

% u  = (yPBS - yhat)./sigma;
% T(1) = normcdf(u,0,1);   %T0
% T(2) = -normpdf(u,0,1);  %T1
% %
% 
% if gIndex == 0
%     GEI = T(1);   %T0
% else
%     for i = 3:gIndex+1
%         k = i-1;
%         T(i) = -(u^(k-1))*normpdf(u,0,1) + (k-1)*T(i-2);
%     end
%     
%     GEI = 0;
%     for i = 1:gIndex+1
%         k    = i-1;
%         ng   = factorial(gIndex)/(factorial(k)*factorial(gIndex-k));
%         dtmp = ((-1)^(k))*ng*(u^(gIndex-k))*T(i);
%         GEI  = GEI + dtmp;
%     end
%     GEI = GEI*(sigma^gIndex);
% end
% 
% 
% GEI = GEI^(1.0/gIndex);

% % run
% u  = (yPBS - yhat)./sigma;
% EI = (yPBS - yhat).*normcdf(u, 0, 1) + sigma.*normpdf(u, 0, 1);

return


