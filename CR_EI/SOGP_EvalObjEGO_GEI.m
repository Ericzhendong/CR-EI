function [J C] = SOGP_EvalObjEGO_GEI(x, objOPTS)

FIT = objOPTS.FIT;
[Yhat PredVar] = GP_eval(FIT,x);

% gIndex = objOPTS.gIndex;
%function GEI = Generalized_ExpectedImprovement(yPBS, yhat, sigma,gIndex)
C = feval(objOPTS.infillFN, objOPTS.Ytarget, Yhat, sqrt(PredVar),objOPTS.gIndex);
J = -C;

return


