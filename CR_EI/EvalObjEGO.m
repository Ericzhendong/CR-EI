function [J C] = EvalObjEGO(x, objOPTS)

% hyp     = objOPTS.FIT.hyp;
% covfunc = objOPTS.FIT.covfunc;
% likfunc = objOPTS.FIT.likfunc;

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
        
    case 'PEI_SOGP'
        
        C = feval(objOPTS.infillFN, objOPTS.Ytarget, Yhat, sqrt(PredVar));
        ex_EI_points = objOPTS.ex_EI_Points;
        %Normalization
        mx  = FIT.nparx(1,:); sx   = FIT.nparx(2,:); 
        my  = FIT.npary(1,:); sy   = FIT.npary(2,:);

        %Normalization of candidate points
        x_prim = x;
        mg     = size(x_prim,1);        
        x_prim_normalized =(x_prim-repmat(mx,mg,1))./repmat(sx,mg,1);
        Normalized_x      = x_prim_normalized;
        
%         Normalized_x      = [x_prim_normalized,x(:,end)];

        %Nomalization of previous EI points       
        nEI    = size(ex_EI_points,1);
        ex_tmp = ex_EI_points;
%         ex_tmp = ex_EI_points(:,1:end-1);
        Normalized_ex_tmp       = (ex_tmp - repmat(mx,nEI,1))./repmat(sx,nEI,1);
        Normalized_ex_EI_points = Normalized_ex_tmp;
%         Normalized_ex_EI_points = [Normalized_ex_tmp,ex_EI_points(:,end)];
        
%         nv = size(x,2)-1;
        nv = size(x,2);
        covfunc    = 'covSEard';
        hyp = FIT.hyp.cov;
%         num_cc_hyp = sum(1:OPT.nTasks);
%         hyp        = FIT.hyp.cov(num_cc_hyp+(1:nv));
        
        correlation = zeros(size(x,1),1);

        for i =1:size(ex_EI_points,1)
            xtmp           = Normalized_ex_EI_points(i,:);
            correlation(:,i) = feval(covfunc,hyp,Normalized_x,xtmp)/hyp(end);   
        end
        C = C.*prod(1-correlation,2);
               
        %function EI = ExpectedImprovement(yPBS, yhat, sigma)
        J = -C; % so the "optimizer" can minimize it
        
        
        
    
    
    case 'PEI'
        
        C = feval(objOPTS.infillFN, objOPTS.Ytarget, Yhat, sqrt(PredVar));
        ex_EI_points = objOPTS.ex_EI_Points;
        %Normalization
        mx  = FIT.nparx(1,:); sx   = FIT.nparx(2,:); 
        my  = FIT.npary(1,:); sy   = FIT.npary(2,:);

        %Normalization of candidate points
        x_prim = x(:,1:end-1);
        mg     = size(x_prim,1);        
        x_prim_normalized =(x_prim-repmat(mx,mg,1))./repmat(sx,mg,1);
        Normalized_x      = [x_prim_normalized,x(:,end)];

        %Nomalization of previous EI points       
        nEI    = size(ex_EI_points,1);
        ex_tmp = ex_EI_points(:,1:end-1);
        Normalized_ex_tmp       = (ex_tmp - repmat(mx,nEI,1))./repmat(sx,nEI,1);
        Normalized_ex_EI_points = [Normalized_ex_tmp,ex_EI_points(:,end)];
        
        nv = size(x,2)-1;
        covfunc    = 'MTGP_covSEardU';
        num_cc_hyp = sum(1:OPT.nTasks);
        hyp        = FIT.hyp.cov(num_cc_hyp+(1:nv));
       
        correlation = zeros(size(x,1),1);

        for i =1:size(ex_EI_points,1)
            xtmp           = Normalized_ex_EI_points(i,:);
            correlation(:,i) = feval(covfunc,hyp,Normalized_x,xtmp);   
        end
        C = C.*prod(1-correlation,2);
               
        %function EI = ExpectedImprovement(yPBS, yhat, sigma)
        J = -C; % so the "optimizer" can minimize it

end

return
