function Pnew= SOGP_IterationEGO(infillFN, targetImprov, ...
    OPT, FIT, normSpace, popsize, nInfillPointsPerCycle, NbVariables, varargin)

    objOPTS.OPT        = OPT;
    objOPTS.FIT        = FIT;
    objOPTS.infillFN   = infillFN;  
    objOPTS.Ytarget    = min(OPT.Y);
           
    X = OPT.X; 
    Pnew = [];
    ex_EI_points = [];
    
%     X_rand = rand(10,1);
%     [Yhat_EI PredVar_EI] = GP_eval(FIT,X);
%     C_EI = feval(infillFN, Ytarget, Yhat_EI, sqrt(PredVar_EI));
%     objOPTS.EI_mean = mean(C_EI);
%     objOPTS.EI_std  = std(C_EI);
    
    % First infill point
    objOPTS.InfillFlag = 'EI'; 
    Infill_criterion = @(x)SOGP_EvalObjEGO(x,objOPTS);
    %[J C] = SOGP_EvalObjEGO(x, objOPTS)
    
    poptm = SBDOInfillDEOptimization(Infill_criterion, objOPTS, normSpace, ...
        popsize,NbVariables);
    if SBDOCheckDataSet([X; poptm], size(X,1) + 1, NbVariables)   % Avoid sample overlapping      
        Pnew         = [Pnew;poptm];
        X            = [X;poptm]; 
        ex_EI_points = [ex_EI_points;poptm];
    else
        Pnew = [];
    end
    
    %2~nInfillPointsPerCycle
    if((nInfillPointsPerCycle>1) && (~isempty(ex_EI_points)))
        
        objOPTS.InfillFlag = 'PEI'; 
        objOPTS.ex_EI_Points = ex_EI_points;
        Infill_criterion = @(x)SOGP_EvalObjEGO(x,objOPTS);
        
        for i = 2:nInfillPointsPerCycle
            
            poptm = SBDOInfillDEOptimization(Infill_criterion, objOPTS, normSpace, ...
                popsize,NbVariables); 
            if SBDOCheckDataSet([X; poptm], size(X,1) + 1, NbVariables)   % Avoid sample overlapping      
                Pnew = [Pnew;poptm];
                X    = [X;poptm];
                ex_EI_points = [ex_EI_points;poptm];
                objOPTS.ex_EI_Points = ex_EI_points;
                Infill_criterion = @(x)SOGP_EvalObjEGO(x,objOPTS);
            end
            
        end
    end
    
    
%     yEI = SOGP_EvalObjEGO(Pnew,objOPTS);
end

