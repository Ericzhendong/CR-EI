function [Pnew yEI]= SOGP_IterationEGO_GEI(infillFN, targetImprov, ...
    OPT, FIT, normSpace, popsize, nInfillPointsPerCycle, NbVariables, varargin)


    if(~isempty(varargin)) 
        objOPTS.gIndex = varargin{1};
    end

    objOPTS.OPT        = OPT;
    objOPTS.FIT        = FIT;
    objOPTS.infillFN   = infillFN;  
    objOPTS.Ytarget    = min(OPT.Y);
           
    X = OPT.X; 
    Pnew = [];
     
    % First infill point
    objOPTS.InfillFlag = 'GEI'; 
    Infill_criterion = @(x)SOGP_EvalObjEGO_GEI(x,objOPTS);
    %function [J C] = SOGP_EvalObjEGO_GEI(x, objOPTS)
    
    poptm = SBDOInfillDEOptimization(Infill_criterion, objOPTS, normSpace, ...
        popsize,NbVariables);
    if SBDOCheckDataSet([X; poptm], size(X,1) + 1, NbVariables)   % Avoid sample overlapping      
        Pnew         = [Pnew;poptm]; 
    else
        Pnew = [];
    end
    
end

