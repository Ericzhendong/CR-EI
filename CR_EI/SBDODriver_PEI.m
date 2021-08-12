function [X, Y, stateEGO] = SBDODriver_PEI(actualFN,designspace, normspace,... 
    OPT, FIT, targetImprov, infillFN, maxNbCycles, SBDOvariant,...
    npointspercycle, varargin)

if(~isempty(varargin)) 
    x_test = varargin{1};
end

%% Initialization
SOGP_OPT      = OPT;
SOGP_FIT      = FIT;
Ybest_Current = min(SOGP_OPT.Y);
Xbest_Current = min(SOGP_OPT.X);

[SOGP_stateEGO SOGP_npoints ndv popsize] = SBDOMakeState(SOGP_OPT);
SOGP_stateEGO.neval(1,1) = size(SOGP_OPT.Y,1);  %Number of samples consumed

% Ybest_Real = -10.5364;
% MaxNbPoints = maxNbCycles*3 + SOGP_stateEGO.neval(1,1);

%% Iteration and Serach Process
timeToStop = 0;
iter       = 1;
MSE = [];
SOGP_infill = [];

% y_test = feval(actualFN, x_test);
    


while timeToStop == 0
    
    %% Examine the accuracy and diversity of accuracy of the models
   
%     yPred_SOGP = GP_eval(SOGP_FIT,x_test);
%     MSE_SOGP  = (sum((yPred_SOGP-y_test).^2)/size(x_test,1));
%     MSE = [MSE;MSE_SOGP];
    
    SOGP_Pnew = [];
    SOGP_Pnew = SOGP_IterationEGO(infillFN, targetImprov, SOGP_OPT, SOGP_FIT, normspace, popsize, npointspercycle, ndv);
    
    
    
    
   
    if(~isempty(SOGP_Pnew))
        SOGP_Xnew = ScaleVariable(SOGP_Pnew, normspace, designspace);   
        SOGP_Ynew = feval(actualFN,SOGP_Xnew);
        
        SOGP_infill = [SOGP_infill; SOGP_Xnew SOGP_Ynew iter*ones(size(SOGP_Ynew,1),1)];
        
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
    
%     nevals_so_far = size(SOGP_OPT.Y,1);    
%     if((Ybest_Current>Ybest_Real)&&(nevals_so_far<MaxNbPoints))
%         timeToStop = 0;
%     else
%         timeToStop = 1;
%     end
    
    timeToStop = iter >= maxNbCycles;
    
        
end

%% Termination

% stateEGO.ReasonToStop = reason;


% X = ScaleVariable([OPT.X; Pnew], normspace, designspace);
% Y = [OPT.Y; Ynew];

X = SOGP_OPT.X;
Y = SOGP_OPT.Y;
[y idx] = min(Y); 
y   = y(1); 
idx = idx(1);
x   = X(idx, : );


% yPred_SOGP = GP_eval(SOGP_FIT,x_test);
% MSE_SOGP  = (sum((yPred_SOGP-y_test).^2)/size(x_test,1));
% MSE = [MSE;MSE_SOGP];


stateEGO = SOGP_stateEGO;
% stateEGO.MSE = MSE;
% stateEGO.SOGP_infill = SOGP_infill;
        
return







