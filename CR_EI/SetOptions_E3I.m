function OPT = SetOptions_E3I(Flag, X, Y)

    OPT.SRGT = Flag;
              
    OPT.X = X;
    OPT.Y = Y;   
    OPT.covfunc = 'covSEiso';       
    OPT.hyp.cov = [0.5*log(0.5);0.1]; 

    OPT.hyp.lik = log(sqrt(1.0e-3));

end

