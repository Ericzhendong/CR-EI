function [y_star, x_star]= Iteration_Thompson_Sampling(TS_OPT, TS_FIT, normSpace)


objOPTS.OPT = TS_OPT;
objOPTS.FIT = TS_FIT;

ndv = size(TS_OPT.X,2);

% Infill_criterion = @(x)ThompsonSampling_Eval(x, objOPTS);
Infill_criterion = @(x)EvalObj_ThompsonSampling(x, objOPTS);

if ndv < 10
    popsize = max(20*ndv, 50);
else
    popsize = 100;
end

x_star = SBDOInfillDEOptimization(Infill_criterion, [], normSpace, popsize,ndv);

y_star = ThompsonSampling_Eval(x_star, objOPTS);

end