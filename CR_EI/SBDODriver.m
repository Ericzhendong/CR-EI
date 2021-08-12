function [X, Y, stateEGO] = SBDODriver_MPSE(actualFN,designspace, normspace,... 
    OPT, FIT, targetImprov, infillFN, maxNbCycles, SBDOvariant,...
    npointspercycle, varargin)

if(~isempty(varargin)) 
    DataFolds = varargin{1};
end

%% EGO parameters

% if(OPT.SRGT=='MTGP')
if(strcmp(OPT.SRGT,'MTGP'))
    
    nTasks = OPT.nTasks;
    X = OPT.X;
    Y = OPT.Y;
          
    [Ybest_Current,nIndicate] = min(OPT.Y);
    nTask_Best    = OPT.X(nIndicate,end);
    Xbest_Current = OPT.X(nIndicate,1:end-1);

    %add in best solution to the related task before the iteration
    X_in = [];
    Y_in = [];
    for h = 1:nTasks
        if(h~=nTask_Best)
            X_in = [X_in;Xbest_Current,h];
            Y_in = [Y_in;Ybest_Current,h];
            OPT.nIndex(h) = OPT.nIndex(h) + 1;
        end
    end

    X_temp = [X;X_in];
    Y_temp = [Y,X(:,end);Y_in];

    X_temp = sortrows(X_temp,size(X_temp,2));
    Y_temp = sortrows(Y_temp,size(Y_temp,2));

    OPT.X  = X_temp;
    OPT.Y  = Y_temp(:,1:end-1);

    FIT    =  MTGP_Fit(OPT,'Default');  %% Directly use the hyper-parameters of the latest training
end

[stateEGO npoints ndv popsize] = SBDOMakeState(OPT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iteration and Serach Process
timeToStop = 0;
iter       = 1;


while timeToStop == 0
    
    switch OPT.SRGT
        
        case 'GP'
            
            Pnew = [];
            Pnew = IterationEGO(infillFN, targetImprov, OPT, FIT, normspace, popsize, npoints, ndv);
            if(~isempty(Pnew))
                Xnew = ScaleVariable(Pnew, normspace, designspace);   
                Ynew = feval(actualFN, Xnew);

                % update state
                stateEGO.Iteration(iter+1,1) = iter;
                stateEGO.Ypbs(iter+1,1)      = min([stateEGO.Ypbs(iter); Ynew]);
                
             else
                Xnew = [];
                Ynew = [];
                stateEGO.Ypbs(iter+1,1)=stateEGO.Ypbs(iter);
            end
            disp(['Best function value for iteration ',num2str(iter),': ',num2str(stateEGO.Ypbs(iter+1,1))]); 
            
        case 'MTGP'
            
            disp(['Iteration',num2str(iter),':']);
            
            X_add = [];
            Y_add = [];
            for i = 1:nTasks            
                nBreakThrough = 0;
                Pnew = [];
                OPT.iTask = i;                              
%                 Pnew = IterationEGO(infillFN, targetImprov, OPT, FIT, normspace, popsize, npoints, ndv);
                Pnew = IterationMultipleSamplesEGO(infillFN, targetImprov, OPT, FIT, normspace, popsize, npoints, ndv);
                
                if(~isempty(Pnew))
                    Xnew = ScaleVariable(Pnew, normspace, designspace);
                    %eval(sprintf('Ynew = Toy_Func%s(Xnew);',num2str(i)));
                    Ynew = feval(actualFN, Xnew);
                    
                    [YnewOpt nNewOpt] = min(Ynew);
                    XnewOpt = Xnew(nNewOpt,:);
                    
                    if(YnewOpt<Ybest_Current)
                        nBreakThrough = 1;
                        Ybest_Current = YnewOpt;
                        nTask_Best    = i;
                    end
                    
                    stateEGO.Ypbs(iter+1,i) = min([stateEGO.Ypbs(iter,i); YnewOpt]);   
                    disp(['Best function value for Task ',num2str(i),': ',num2str(stateEGO.Ypbs(iter+1,i))]); 

                    % Update MTGP
                    X_add   = [Xnew,ones(size(Xnew,1),1)*i];
                    Y_add   = [Ynew,ones(size(Xnew,1),1)*i];
                    OPT.nIndex(i) = OPT.nIndex(i)+size(Xnew,1);  %%Prepare for infill more than one samples
                                        
                    %leak_in strategy
                    if(nBreakThrough)
                        for k = 1:nTasks
                            if(k~=nTask_Best)
                                X_add   = [X_add;XnewOpt,k];
                                Y_add   = [Y_add;YnewOpt,k];
                                OPT.nIndex(k) = OPT.nIndex(k) + 1;  %only the best optimal leak in 
                            end
                        end
                    end


                    X      = [OPT.X;X_add];
                    Y_Prim = [OPT.Y,OPT.X(:,end);Y_add];

                    X      = sortrows(X,size(X,2));
                    Y_Prim = sortrows(Y_Prim,size(Y_Prim,2));
                    OPT.X  = X;
                    OPT.Y  = Y_Prim(:,1:end-1);
                                                           
                    FIT    =  MTGP_Fit(OPT,'Default');  %% Directly use the hyper-parameters of the latest training
                else                      
                    stateEGO.Ypbs(iter+1,i) = stateEGO.Ypbs(iter,i);
                end
                  
               
            end 
          
            stateEGO.Iteration = iter;  
            
            xtmp    = 0:0.01:1;
            [X1,Y1] = meshgrid(xtmp);
            x_tmp   = [X1(:),Y1(:)];
            nTest   = size(x_tmp,1);
            xtest   = repmat(designspace(1,:),nTest,1) + repmat(designspace(2,:)-designspace(1,:),nTest,1).*x_tmp;

            ytest   = feval(actualFN,xtest);    
            nTest   = size(xtest,1);
            yPred   = GP_eval(FIT,xtest);

            RMSE    = sqrt(sum((ytest-yPred).^2)/nTest);
            R2      = 1 - sum((ytest-yPred).^2)/sum((ytest-mean(ytest)).^2);
            fprintf(' current RMSE: %f, Current R2: %f\n', RMSE, R2);
            
            
    end
   
   
    if ~strcmp(OPT.SRGT, 'KRG')
        FIT.KRG_DACEModel = [];
    end
    
    
    [timeToStop reason] = SBDOIsItTimeToStop(FIT.KRG_DACEModel, iter, maxNbCycles);
    
    if ~timeToStop
        % update surrogates
        switch OPT.SRGT
            
            case 'GP'
                
                OPT.X = [OPT.X; Xnew];
                OPT.Y = [OPT.Y; Ynew];               

                FIT = GPFit(OPT);
                
%                 xTest = (0:0.01:1)';
%                 nTest   = size(xTest,1);
%                 yTest = feval(actualFN,xTest);
%                 [Yhat_Test PredVar] = GP_eval(FIT,xTest);
%                 RMSE = sqrt(sum((yTest-Yhat_Test).^2)/nTest);
%                 R2   = 1 - sum((yTest-Yhat_Test).^2)/sum((yTest-mean(yTest)).^2);
%                 
%                 figure();
%                 plot(xTest,yTest,'k-');
%                 hold on;
%                 scatter(OPT.X,OPT.Y);
%                 plot(xTest,Yhat_Test,'b--');
% 
%                 figure();
%                 yPBS = min(OPT.Y);
%                 pEI  = ExpectedImprovement(yPBS, Yhat_Test, sqrt(PredVar));
%                 plot(xTest,pEI);
                
            case 'MTGP'
                
                %X = OPT.X
                %Y = OPT.Y
                FIT = MTGP_Fit(OPT,'Default');
               
        end
        
            
%         npoints = npoints + npointspercycle;
        npoints = npoints + size(Pnew,1);


    end
    
    iter = iter + 1;
end

%% Termination

stateEGO.ReasonToStop = reason;


% stateEGRA.ReasonToStop = reason;
% 
% X = srgtsScaleVariable([srgtOPT.P; Pnew], normspace, designspace);
% Y = [srgtOPT.T; Ynew];
% 
% stateEGRA.X = X;
% stateEGRA.Y = Y;
% 
% if isequal(display, 'ON')
%     strFinal = 'Sequential sampling successfully completed.';
%     fprintf('%s\n\n',strFinal);
% end



switch OPT.SRGT
   case 'GP' 
       X = ScaleVariable([OPT.X; Pnew], normspace, designspace);
       Y = [OPT.Y; Ynew];

       stateEGO.X = X;
       stateEGO.Y = Y;

       [y idx] = min(Y); y = y(1); idx = idx(1);
       x = X(idx, : );
       
    case 'MTGP'   
        X = OPT.X;
        Y = OPT.Y;
        stateEGO.X = X;
        stateEGO.Y = Y;
end
        
return


