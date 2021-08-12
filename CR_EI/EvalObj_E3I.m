function [J C] = EvalObj_E3I(x, y_stars, GP_FIT)

C = E3I(x, y_stars, GP_FIT);
J = -C;  %so the "optimizer" can minimize it

end