function [elementData] = getElementStiffness(inputData)

numberofElement = inputData.numberofElement;
nodeCoordinate = inputData.nodeCoordinate;
nodeConnectivity = inputData.nodeConnectivity;
DOFperNode = inputData.DOFperNode;
materialData = inputData.materialData;

for en = 1:numberofElement
    % Generate global stiffness matrix of element
    elementNodes = nodeConnectivity(en,:);
    for ii = 1:3
       p(ii,:) = nodeCoordinate(elementNodes(ii),:);
    end
    
    % Construct B matrix
    xi = p(1,1);
    xj = p(2,1);
    xm = p(3,1);
    yi = p(1,2);
    yj = p(2,2);
    ym = p(3,2);
    
    A = (xi*(yj-ym) + xj*(ym-yi)+xm*(yi-yj))/2;
    ai = xj*ym - yj*xm;
    aj = yi*xm - xi*ym;
    am = xi*yj - yi*xj;
    bi = yj-ym;
    bj = ym-yi;
    bm = yi-yj;
    ci = xm-xj;
    cj = xi-xm;
    cm = xj-xi;
    
    B = [bi 0 bj 0 bm 0;
        0 ci 0 cj 0 cm;
        ci bi cj bj cm bm]/(2*A);
    
    % Construct D matrix
    E = materialData(en,1);
    nu = materialData(en,2);
    D = (E/(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];    
    
    localKmatrix = B'*D*B*1*A;
    stress_DB = D*B;
    
    globalKmatrix_element = localKmatrix;
    
    % Index matrix with global index
    global_idx = [];
    for i = 1:length(elementNodes)
        node = elementNodes(i);
        global_idx = [ global_idx DOFperNode*(node-1)+1 : DOFperNode*node ];
    end
    globalKdata(en,:) = { global_idx, globalKmatrix_element };
    localKdata(en,:) = { elementNodes, localKmatrix global_idx stress_DB};
end

elementData.globalKdata = globalKdata;
elementData.localKdata = localKdata;


