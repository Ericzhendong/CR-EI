function [J C] = EvalObj_ThompsonSampling(x, objOPTS)

C = ThompsonSampling_Eval(x, objOPTS);
J = C;  %so the "optimizer" can minimize it

end
