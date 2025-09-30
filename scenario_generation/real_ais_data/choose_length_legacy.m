function [inds_chosen,duration] = choose_length_legacy(data)
    size1 = @(A) size(A,1);
    length_matrix = cellfun(size1,data);
    [inds_chosen,duration] = most_occuring(length_matrix);
end 
