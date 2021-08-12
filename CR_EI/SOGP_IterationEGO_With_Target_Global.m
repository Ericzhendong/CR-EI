function [Pnew yEI]= SOGP_IterationEGO_With_Target_Global(infillFN, targetImprov, ...
    OPT, FIT, normSpace, popsize, nInfillPointsPerCycle, NbVariables, varargin)


    if(~isempty(varargin)) 
        Xbest_Current = varargin{1};
        X = varargin{2};
        Distance_Therehold = varargin{3};
    end


    objOPTS.OPT        = OPT;
    objOPTS.FIT        = FIT;
    objOPTS.infillFN   = infillFN;  
    objOPTS.Ytarget    = targetImprov;
           
%     X = OPT.X; 
    Pnew = [];
    ex_EI_points = [];
    
    
    % First infill point
    objOPTS.InfillFlag = 'EI'; 
    Infill_criterion = @(x)SOGP_EvalObjEGO_Rescaled_Global(x,objOPTS, Xbest_Current, Distance_Therehold);
%     Infill_criterion = @(x)SOGP_EvalObjEGO(x,objOPTS);
    
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
        Infill_criterion = @(x)SOGP_EvalObjEGO_Rescaled_Global(x,objOPTS, Xbest_Current, Distance_Therehold);
%         Infill_criterion = @(x)SOGP_EvalObjEGO(x,objOPTS);

        
        for i = 2:nInfillPointsPerCycle
            
            poptm = SBDOInfillDEOptimization(Infill_criterion, objOPTS, normSpace, ...
                popsize,NbVariables); 
            if SBDOCheckDataSet([X; poptm], size(X,1) + 1, NbVariables)   % Avoid sample overlapping      
                Pnew = [Pnew;poptm];
                X    = [X;poptm];
                ex_EI_points = [ex_EI_points;poptm];
                objOPTS.ex_EI_Points = ex_EI_points;
                Infill_criterion = @(x)SOGP_EvalObjEGO_Rescaled_Global(x,objOPTS, Xbest_Current, Distance_Therehold);
%                 Infill_criterion = @(x)SOGP_EvalObjEGO(x,objOPTS);
            end
            
        end
    end
    
    
%     yEI = SOGP_EvalObjEGO(Pnew,objOPTS);
end

