estAssoc = cell(1,3);
estAssoc{1} = [1 0 1 1];
estAssoc{2} = [0 1 0 1];
estAssoc{3} = [1 1 1 0];

realAssoc = [0 1 0 1];

weightsTable = [0.5; 0.3; 0.2]; % Define weights for the associations

[errorMat,scoreTable,scoreMean] = evaluate_associations(estAssoc,realAssoc,weightsTable);