function [X, Y, stateEGO] = SBDODriver_CR_EI(actualFN,designspace, normspace,... 
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

nBreakthrough = 0;
nFlag_to_Start_local_Rescaled_EI = 0;
Distance_Therehold = 0.1*sqrt(size(Xbest_Current,2));

while timeToStop == 0
    
    SOGP_Xnew_Collect = [];
    SOGP_Ynew_Collect = [];
    
    SOGP_Pnew = [];
    SOGP_Pnew= SOGP_IterationEGO(infillFN, targetImprov, SOGP_OPT, SOGP_FIT, normspace, popsize, npointspercycle, ndv);
    
    if(~isempty(SOGP_Pnew))
        %disp(['EI-based search:']);
        SOGP_Xnew = ScaleVariable(SOGP_Pnew, normspace, designspace);
                
        [yhat_New sigma_New] = GP_eval(SOGP_FIT,SOGP_Xnew);

        %if((0==nBreakthrough)|(yhat_New<prctile(SOGP_OPT.Y, 25))) 
        if((0==nBreakthrough)|(yhat_New<prctile(SOGP_OPT.Y, 50)))  % Setting "yhat_New<prctile(SOGP_OPT.Y, 25)" sometimes can get even better convergence rate 
            
            nFlag_to_Start_local_Rescaled_EI = 0;
            SOGP_Ynew = feval(actualFN,SOGP_Xnew);
            
            SOGP_Xnew_Collect = [SOGP_Xnew_Collect; SOGP_Xnew];
            SOGP_Ynew_Collect = [SOGP_Ynew_Collect; SOGP_Ynew];   

            if(SOGP_Ynew<Ybest_Current)  
                Ybest_Current = SOGP_Ynew;
                Xbest_Current = SOGP_Xnew;
                nBreakthrough = nBreakthrough + 1;
            end
            % update state
            SOGP_stateEGO.Ypbs(iter+1) = min([SOGP_stateEGO.Ypbs(iter); Ybest_Current]);
        else
            SOGP_stateEGO.Ypbs(iter+1) = SOGP_stateEGO.Ypbs(iter);
            nFlag_to_Start_local_Rescaled_EI = 1;
        end
        
        disp(['Best function value based on EI for iteration ',num2str(iter),': ',num2str(SOGP_stateEGO.Ypbs(iter+1))]);

    else
        SOGP_stateEGO.Ypbs(iter+1) = SOGP_stateEGO.Ypbs(iter);
    end
    

    if(nBreakthrough > 0)

        SOGP_OPT_New = SOGP_OPT;
        SOGP_FIT_New = SOGP_FIT;
        Pnew = [];
        yEI = [];
        
        X_to_Check_Overlapping = [SOGP_Xnew_Collect; SOGP_OPT.X];

        %% Set the target improvement in the Rescaled EI
        
        yPBS = min(SOGP_OPT.Y);
        if(yPBS<0) 
            yPBS_Rescaled = min(SOGP_OPT.Y) * 0.9;
        else
            yPBS_Rescaled = min(SOGP_OPT.Y) * 1.1;
        end
        
        y_Distance = sortrows([abs(SOGP_OPT.Y-yPBS_Rescaled) SOGP_OPT.Y], 1);
        
        targetImprov = y_Distance(1,2);
        if(targetImprov==yPBS)
            targetImprov = y_Distance(2,2);
        end

        %% Rescaled EI-based search
        if(~nFlag_to_Start_local_Rescaled_EI) %% Global Search only
            
            %disp(['Rescaled EI-based global search:']);
            % Rescaled EI-based global search
            Pnew_Global = SOGP_IterationEGO_With_Target_Global(infillFN, targetImprov, SOGP_OPT_New, SOGP_FIT_New, normspace, popsize, 1, ndv, Xbest_Current, X_to_Check_Overlapping, Distance_Therehold);
            
            if(~isempty(Pnew_Global)) 
                
                Xnew_Global = ScaleVariable(Pnew_Global, normspace, designspace);

                Ynew_Global = feval(actualFN,Xnew_Global);
                
                if(Ynew_Global<Ybest_Current)  
                    Ybest_Current = Ynew_Global;
                    Xbest_Current = Xnew_Global;
                end

                SOGP_Xnew_Collect = [SOGP_Xnew_Collect; Xnew_Global];
                SOGP_Ynew_Collect = [SOGP_Ynew_Collect; Ynew_Global];

                SOGP_stateEGO.Ypbs(iter+1) = min([SOGP_stateEGO.Ypbs(iter); Ybest_Current]);
                disp(['Best function value based on Rescaled Global Search for iteration ',num2str(iter),': ',num2str(SOGP_stateEGO.Ypbs(iter+1))]);
            end
            
        else %% Both Global and Local Search
            %disp(['Both global and local search:']);
                
            %% Local Search
            Pnew_Local = SOGP_IterationEGO_With_Target_Local(infillFN, targetImprov, SOGP_OPT_New, SOGP_FIT_New, normspace, popsize, 1, ndv, Xbest_Current, X_to_Check_Overlapping);
            
            if(~isempty(Pnew_Local))  % Evaluate if and only if the corresponding sample location close to the presetn best solution
                Xnew_Local = ScaleVariable(Pnew_Local, normspace, designspace);
                Distance_to_PBS = sqrt(sum((Xnew_Local - Xbest_Current).^2, 2));
   
                if(Distance_to_PBS<Distance_Therehold)
                    Ynew_Local = feval(actualFN,Xnew_Local);
                    
                    if(Ynew_Local<Ybest_Current)  
                        Ybest_Current = Ynew_Local;
                        Xbest_Current = Xnew_Local;
                         
                    end

                    SOGP_Xnew_Collect = [SOGP_Xnew_Collect; Xnew_Local];
                    SOGP_Ynew_Collect = [SOGP_Ynew_Collect; Ynew_Local];

                    SOGP_stateEGO.Ypbs(iter+1) = min([SOGP_stateEGO.Ypbs(iter); Ybest_Current]);
                    disp(['Best function value based on Rescaled local Search for iteration ',num2str(iter),': ',num2str(SOGP_stateEGO.Ypbs(iter+1))]);
                    
                end
            end
            
            %% Global Search
            
            X_to_Check_Overlapping = [SOGP_Xnew_Collect; SOGP_OPT.X];
            Pnew_Global = SOGP_IterationEGO_With_Target_Global(infillFN, targetImprov, SOGP_OPT_New, SOGP_FIT_New, normspace, popsize, 1, ndv, Xbest_Current, X_to_Check_Overlapping, Distance_Therehold);
            
            if(~isempty(Pnew_Global)) 
 
                Xnew_Global = ScaleVariable(Pnew_Global, normspace, designspace);

                Ynew_Global = feval(actualFN,Xnew_Global);
                if(Ynew_Global<Ybest_Current)  
                    Ybest_Current = Ynew_Global;
                    Xbest_Current = Xnew_Global;
                end

                SOGP_Xnew_Collect = [SOGP_Xnew_Collect; Xnew_Global];
                SOGP_Ynew_Collect = [SOGP_Ynew_Collect; Ynew_Global];

                SOGP_stateEGO.Ypbs(iter+1) = min([SOGP_stateEGO.Ypbs(iter); Ybest_Current]);
                disp(['Best function value based on Rescaled Global Search for iteration ',num2str(iter),': ',num2str(SOGP_stateEGO.Ypbs(iter+1))]);
            end

        end
  
    end
    
    %update SOGP
    SOGP_OPT.X = [SOGP_OPT.X;SOGP_Xnew_Collect];
    SOGP_OPT.Y = [SOGP_OPT.Y;SOGP_Ynew_Collect];
    SOGP_FIT   = GPFit(SOGP_OPT);

    SOGP_stateEGO.neval(iter+1,1) = size(SOGP_OPT.Y,1);  %Number of samples consumed
    iter = iter + 1;
       
    timeToStop = iter >= maxNbCycles;
    
    disp(['Number of samples consumed: ',num2str(size(SOGP_OPT.Y,1))]);
        
end

%% Termination

X = SOGP_OPT.X;
Y = SOGP_OPT.Y;

stateEGO = SOGP_stateEGO;
     
return







