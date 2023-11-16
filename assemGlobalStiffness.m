function [globalKmatrix] = assemGlobalStiffness(elementData, inputData)

numberofElement = inputData.numberofElement;
totalDOF = inputData.totalDOF;
globalKdata = elementData.globalKdata;

% Set size of golbal stiffness matrix
globalKmatrix = zeros(totalDOF);

% Assemble golbal stiffness of elements
for en = 1:numberofElement                 % en: a certain element number
    global_idx = cell2mat(globalKdata(en,1));
    globalKmatrix_element = cell2mat(globalKdata(en,2));
    
    for local_row = 1:length(global_idx)
        global_row = global_idx(local_row);
        for local_col = 1:length(global_idx)
            global_col = global_idx(local_col);
            globalKmatrix(global_row,global_col) = globalKmatrix(global_row,global_col)...
                + globalKmatrix_element(local_row,local_col);
        end
    end
    
end

