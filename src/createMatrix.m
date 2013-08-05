function globalMatrix = createMatrix(localMatrix, geomMatrix)
    n = max(max(geomMatrix));
    
    globalMatrix = zeros(n);
    
    for i = 1:n
        for j = 1:n
            [rows1, cols1] = find(geomMatrix == i);
            [rows2, cols2] = find(geomMatrix == j);
            
            rows = intersect(rows1, rows2);
            
            if length(rows) > 0
                value = 0;

                for row = rows'
                    iindex = find(geomMatrix(row,:) == i);
                    jindex = find(geomMatrix(row,:) == j);

                    value = value + localMatrix(iindex, jindex);
                end

                globalMatrix(i, j) = value;
            end
        end
    end
end

