function [errorMat,scoreTable,scoreMean] = evaluate_associations(estAssocs,realAssoc,weightsTable)
    errorCell = cellfun(@(assocTable) evaluate_table(assocTable,realAssoc),estAssocs,'UniformOutput',false);
    %errorMat = cell2mat(errorCell);
    errorMat = [errorCell{:}];
    scoreTable = sum(errorMat,1);
    scoreMean = sum(scoreTable.*weightsTable');
end

function errorTable = evaluate_table(assocTable,realAssoc)
    diffTable = realAssoc - assocTable;
    errorTable = (diffTable ~= 0);
end