function [Pnew]= SOGP_Iteration_E3I(GP_OPT, GP_FIT, normSpace, popsize, NbVariables, varargin)

    %% Generate the M Thompson sample optima
    M = 100;
    y_stars_normalized = zeros(M,1);
    x_stars = zeros(M, size(normSpace, 2));
    
    for i = 1:M
        [TS_OPT, TS_FIT] = ThompsonSampling_Fit(GP_FIT);
        
        %% plot
        
        
        objOPTS.OPT = TS_OPT;
        objOPTS.FIT = TS_FIT;
        
%         x_test = (0:0.01:1)';
%         y_pred_TS_Norm = ThompsonSampling_Eval(x_test, objOPTS);
        
%         my = GP_FIT.npary(1,:); sy = GP_FIT.npary(2,:);   
%         y_pred_TS = y_pred_TS_Norm*sy + my;
%         
%         figure();
%         plot(x_test, y_pred_TS);
%         hold on;
%         scatter(GP_OPT.X', GP_OPT.Y');

        [y_stars_normalized(i), x_stars(i,:)] = Iteration_Thompson_Sampling(TS_OPT, TS_FIT, normSpace);
        
%         y_tmp_norm = y_stars_normalized(i);
%         y_tmp = y_tmp_norm*sy + my;
           
    end
    
    
    %% denormalization
    my = GP_FIT.npary(1,:); sy = GP_FIT.npary(2,:);   
    y_stars = y_stars_normalized*sy + my;


    %% 
    %Infill_criterion = @(x)E3I(x, y_stars, GP_FIT)
    
    
    X = GP_OPT.X; 
    
    Infill_criterion = @(x)EvalObj_E3I(x, y_stars, GP_FIT)
   
    poptm = SBDOInfillDEOptimization(Infill_criterion, [], normSpace, popsize, NbVariables);
    
    if SBDOCheckDataSet([X; poptm], size(X,1) + 1, NbVariables)   % Avoid sample overlapping      
        Pnew  = poptm;
 
    else
        Pnew = [];
    end

end

