function rowMatrix = computeRowsOfToeplitzHankelMatrix(rowNumber,...
        nColumns, crossCorrelationVectors, hankelMatrixIsAdded, dcIsIncluded)
    if rowNumber == 1
        toeplitzRows = crossCorrelationVectors(1:nColumns,:);
    else
        toeplitzRows = ...
            [flip(crossCorrelationVectors(2:rowNumber,:),1);...
            crossCorrelationVectors(1:nColumns-rowNumber+1,:)];
    end
    
    
    if dcIsIncluded && hankelMatrixIsAdded
        hankelOffset = 1;
    else
        hankelOffset = 3;
    end
    
    hankelRows = crossCorrelationVectors((0:nColumns-1)+...
                                         hankelOffset+rowNumber-1,:);

    if hankelMatrixIsAdded
        rowMatrix = toeplitzRows + hankelRows;
    else
        rowMatrix = toeplitzRows - hankelRows;
    end

end