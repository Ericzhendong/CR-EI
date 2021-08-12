function nFlag = Check_Overlapping(Pnew,Points)
    
    D = Points - repmat(Pnew,size(Points,1),1);
    % nFlag = min(sum(abs(D),2)) ~= 0;
    nFlag = (min(sum(abs(D),2)) >= 0.01);

end

