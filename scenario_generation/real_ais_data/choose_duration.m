function [inds_chosen,duration] = choose_duration(data,mode)
    durationFun = @(A) A(end,5);
    durationMatrix = cellfun(durationFun,data);
    if mode == "max"
        duration = max(durationMatrix);
    elseif mode=="mean"
        duration = mean(durationMatrix);
    elseif mode=="quartile"
        [r,q] = iqr(durationMatrix);
        duration = q(1);
    else
        disp("Erreur: choose mode among 'mean' and 'max'")
    end 
    inds_chosen = find(durationMatrix > duration);  
end 
