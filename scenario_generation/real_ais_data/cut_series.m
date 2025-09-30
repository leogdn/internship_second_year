function dataCut = cut_series(data,duration)
    %cut_mat = @(A) A(A(:,5)<=duration,:);
    dataCutTemp = cellfun(@(mat) cut_mat(mat,duration),data,'UniformOutput',false);
    keepArray = inds_to_keep(dataCutTemp);
    dataCut = dataCutTemp(keepArray);
    % for i=1:nCut
    %     if size(dataCutTemp(i),1)<=2
    %         dataCut(i) = [];
    %     end
    % end
    
end

function keepArray = inds_to_keep(data)
    length_short_mat = (@(A) size(A,1) > 2);
    keepArray = cellfun(length_short_mat, data);
end

function cuttedMat = cut_mat(mat,duration)
    indMax = find(mat(:,5)>duration, 1 );
    cuttedMat = mat(1:indMax,:);
end