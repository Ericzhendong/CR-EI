function [X, Y, stateEGO] = SBDODriver_GEI(actualFN,designspace, normspace,... 
    OPT, FIT, targetImprov, infillFN, maxNbCycles, SBDOvariant,...
    npointspercycle, varargin)


%% Initialization
SOGP_OPT      = OPT;
SOGP_FIT      = FIT;
Ybest_Current = min(SOGP_OPT.Y);
Xbest_Current = min(SOGP_OPT.X);

[SOGP_stateEGO SOGP_npoints ndv popsize] = SBDOMakeState(SOGP_OPT);
SOGP_stateEGO.neval(1,1) = size(SOGP_OPT.Y,1);  %Number of samples consumed

%% Iteration and Serach Process
timeToStop = 0;
iter       = 1;
MSE = [];
SOGP_infill = [];


while timeToStop == 0
    
    if(iter<5) 
        gIndex = 20;
    elseif(iter<10)
        gIndex = 10;
    elseif(iter<20)
        gIndex = 5;
    elseif(iter<25)
        gIndex = 2;
    elseif(iter<35)
        gIndex = 1;
    else
        gIndex = 0;
    end
       
    
    SOGP_Pnew = [];
    SOGP_Pnew = SOGP_IterationEGO_GEI(infillFN, targetImprov, SOGP_OPT, SOGP_FIT, normspace, popsize, npointspercycle, ndv, gIndex);
   
    if(~isempty(SOGP_Pnew))
        SOGP_Xnew = ScaleVariable(SOGP_Pnew, normspace, designspace);   
        SOGP_Ynew = feval(actualFN,SOGP_Xnew);
                
        [Ybest_tmp,nbest_tmp] = min(SOGP_Ynew);              
        if(Ybest_tmp<Ybest_Current)  
            Ybest_Current = Ybest_tmp;
            Xbest_Current = SOGP_Xnew(nbest_tmp,:);
        end
                
        %update SOGP samples
        SOGP_OPT.X = [SOGP_OPT.X;SOGP_Xnew];
        SOGP_OPT.Y = [SOGP_OPT.Y;SOGP_Ynew];
        SOGP_FIT   = GPFit(SOGP_OPT);
                
        % update state
        SOGP_stateEGO.Ypbs(iter+1) = min([SOGP_stateEGO.Ypbs(iter); Ybest_Current]);
    else
        SOGP_stateEGO.Ypbs(iter+1) = SOGP_stateEGO.Ypbs(iter);
    end
    disp(['Best function value for iteration of SOGP ',num2str(iter),': ',num2str(SOGP_stateEGO.Ypbs(iter+1))]);
    
    SOGP_stateEGO.neval(iter+1,1) = size(SOGP_OPT.Y,1);  %Number of samples consumed
    iter = iter + 1;
        
    timeToStop = iter >= maxNbCycles;
          
end

%% Termination

X = SOGP_OPT.X;
Y = SOGP_OPT.Y;
[y idx] = min(Y); 
y   = y(1); 
idx = idx(1);
x   = X(idx, : );

stateEGO = SOGP_stateEGO;
        
return







