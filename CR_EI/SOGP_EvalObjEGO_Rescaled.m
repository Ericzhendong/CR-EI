function [J C] = SOGP_EvalObjEGO_Rescaled(x, objOPTS, varargin)

if(~isempty(varargin)) 
    Xbest_Current = varargin{1};
end


switch objOPTS.OPT.SRGT
    case 'GP'        
         FIT = objOPTS.FIT;
         [Yhat PredVar] = GP_eval(FIT,x);
         
    case 'MTGP'        
        OPT = objOPTS.OPT;
        FIT = objOPTS.FIT;
        x = [x objOPTS.OPT.iTask*ones(size(x,1),1)];
        [Yhat PredVar] = MTGP_eval(OPT,FIT,x,'Default');

    case 'Swersky'       
        MTGP_OPT = objOPTS.OPT;
        MTGP_FIT = objOPTS.FIT;
        x = [x objOPTS.OPT.iTask];
        
        %function [varargout] = AveragingMTGPrediction(Fit, inf,xs)
        [Yhat PredVar] = AveragingMTGPrediction(MTGP_FIT,@MTGP_infExact,x);
               
end


% Acquisition function
switch objOPTS.InfillFlag
    case 'Fmin'
        J = Yhat;
        
    case 'maxMSE'
        J = -PredVar;
        
    case 'GEI'      
        gIndex = objOPTS.gIndex;
        %function GEI = Generalized_ExpectedImprovement(yPBS, yhat, sigma,gIndex)
        C = feval(objOPTS.infillFN, objOPTS.Ytarget, Yhat, sqrt(PredVar),objOPTS.gIndex);
        J = -C;
        
    case 'EI'   
        C = feval(objOPTS.infillFN, objOPTS.Ytarget, Yhat, sqrt(PredVar));
        J = -C; % so the "optimizer" can minimize it
        
    case 'PEI'
        
        C = feval(objOPTS.infillFN, objOPTS.Ytarget, Yhat, sqrt(PredVar));
        ex_EI_points = objOPTS.ex_EI_Points;
        %Normalization
        mx  = FIT.nparx(1,:); sx   = FIT.nparx(2,:); 
        my  = FIT.npary(1,:); sy   = FIT.npary(2,:);

        %Normalization of candidate points
        mg = size(x,1);        
        Normalized_x  =(x-repmat(mx,mg,1))./repmat(sx,mg,1);

        %Nomalization of previous EI points       
        nEI    = size(ex_EI_points,1);
        Normalized_ex_EI_points = (ex_EI_points - repmat(mx,nEI,1))./repmat(sx,nEI,1);
        
        
        nv = size(x,2);
        covfunc = 'SOGP_covSEardU';
        hyp     = FIT.hyp.cov(1:nv);
       
        correlation = zeros(size(x,1),1);

        for i =1:size(ex_EI_points,1)
            xtmp           = Normalized_ex_EI_points(i,:);
            correlation(:,i) = feval(covfunc,hyp,Normalized_x,xtmp);   
        end
        C = C.*prod(1-correlation,2);
               
        %function EI = ExpectedImprovement(yPBS, yhat, sigma)
        J = -C; % so the "optimizer" can minimize it

end

Y_Filter = prctile(objOPTS.OPT.Y, 25);
Y_Scaling = (Yhat<Y_Filter);
Distances = sqrt(sum((x - repmat(Xbest_Current, size(x,1), 1)).^2, 2));
% D_Scaling = (Distances>0.1);
Therehold = 0.1*sqrt(size(x,2));
D_Scaling = (Distances>Therehold);
J = J.*Y_Scaling.*D_Scaling; 

return


