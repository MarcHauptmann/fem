function globalMatrix = createMatrix(localMatrix, geomMatrix, h)
% CREATEMATRIX erstellt globale Matrix
%
%   globalMatrix = CREATEMATRIX(localMatrix, geomMatrix)
%   globalMatrix = CREATEMATRIX(localMatrix, geomMatrix, h)
%
%
    if not(exist('h', 'var'))
        h = 1;
    end
    
    n = max(max(geomMatrix));               % Anzahl der globalen Basisfunktionen
    elementCount = size(geomMatrix, 1);     % Anzahl der Elemente
    functionCount = size(geomMatrix, 2);    % Anzahl der lokalen Basisfunktionen
    
    % Faktoren (Länge) der Elemente
    if length(h) == 1
        hVector = h * ones(elementCount,1);
    else
        hVector = h;
    end
    
    globalMatrix = zeros(n);
    
    
    
    for elementIndex = 1:elementCount
        for localPhi1 = 1:functionCount
            for localPhi2 = 1:functionCount
                phi1 = geomMatrix(elementIndex, localPhi1);
                phi2 = geomMatrix(elementIndex, localPhi2);

                value = hVector(elementIndex) * localMatrix(localPhi1, localPhi2);

                globalMatrix(phi1, phi2) = globalMatrix(phi1, phi2) + value;
            end
        end
    end
end
