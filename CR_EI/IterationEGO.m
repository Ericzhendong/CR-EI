function Pnew = IterationEGO(infillFN, targetImprov, ...
    OPT, FIT, normSpace, popsize, NbPoints, NbVariables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings

    switch OPT.SRGT
        case 'GP'
            objOPTS.Ytarget    = min(OPT.Y);
            objOPTS.OPT        = OPT;
            objOPTS.FIT        = FIT;
            objOPTS.infillFN   = infillFN;
            objOPTS.InfillFlag = 'EI';   

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            %Infill_criterion = @(x)EvalObjEGO(x,objOPTS);
            Infill_criterion = @(x)EvalObjEGO(x,objOPTS);         
            poptm = SBDOInfillDEOptimization(Infill_criterion, objOPTS, normSpace, ...
                popsize,NbVariables);
            
            if SBDOCheckDataSet([OPT.X; poptm], size(OPT.X,1) + 1, NbVariables)
                Pnew = poptm;
            else
                Pnew = [];
            end


        case 'MTGP'

            objOPTS.OPT        = OPT;
            objOPTS.FIT        = FIT;
            objOPTS.infillFN   = infillFN;
            objOPTS.InfillFlag = 'EI';   
            
            X      = OPT.X;
            i      = OPT.iTask;
            nTasks = OPT.nTasks;
            
            nIndex = find(X(:,end)==i);
            objOPTS.Ytarget = min(OPT.Y(nIndex)); 
            
            
            Infill_criterion = @(x)EvalObjEGO(x,objOPTS);
            poptm = SBDOInfillDEOptimization(Infill_criterion, objOPTS, normSpace, ...
                popsize,NbVariables);
            
            X_tmp = X(:,1:end-1);
            X_unique = unique(X_tmp,'rows');
            if SBDOCheckDataSet([X_unique; poptm], size(X_unique,1) + 1, NbVariables)
                Pnew = poptm;
            else
                Pnew = [];
            end
            
    end
        
end



