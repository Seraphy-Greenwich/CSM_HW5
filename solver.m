function [outputData] = solver(globalKmatrix, elementData, inputData)

numberofElement = inputData.numberofElement;
constrainedDispl_globalidx = inputData.constrainedDispl_globalidx;
constrainedDispl = inputData.constrainedDispl;
constrainedForce = inputData.constrainedForce;
totalDOF = inputData.totalDOF;
localKdata = elementData.localKdata;
globalKdata= elementData.globalKdata;
central_Elements = inputData.central_Elements;
initalData = inputData.initalData;

% get reduced matrix      
reduced_Force = constrainedForce;
reduced_Force(constrainedDispl_globalidx,:) = [];

reduced_K = globalKmatrix;
reduced_K(constrainedDispl_globalidx,:) = [];
nozero = find(constrainedDispl);
reduced_Force = reduced_Force - reduced_K(:,nozero)*constrainedDispl(nozero);
reduced_K(:,constrainedDispl_globalidx) = [];
freeDispl = reduced_K \ reduced_Force;

% get global displacement and force
globalDispl = constrainedDispl;
freeDispl_globalidx = setdiff((1:totalDOF), constrainedDispl_globalidx);
globalDispl(freeDispl_globalidx) = freeDispl;
globalForce = globalKmatrix*globalDispl;

% get stress data
for en = 1:numberofElement
    elementNodes = cell2mat(localKdata(en,1));
    stress_DB = cell2mat(localKdata(en,4));
    
    global_idx = cell2mat(globalKdata(en,1));
    globalDispl_element = globalDispl(global_idx);
    
    stress(en,:) = stress_DB*globalDispl_element;
end

% get a stress concentration factor
[stress_centralmax,I] = max(stress(central_Elements,1));
stress_elementmax = central_Elements(I); %[N/mm^2]

F = initalData(4);
d = initalData(2);
rho = initalData(3);
stress_nominal = F/d; %[N/mm^2]

S_cc = stress_centralmax/stress_nominal;
ratio = rho/d;

% save output data
% outputData.globalDispl = globalDispl;
% outputData.globalForce = globalForce;
% outputData.globalKmatrix = globalKmatrix;
outputData.stress = stress;
outputData.stress_centralmax = stress_centralmax;
outputData.stress_elementmax = stress_elementmax;
outputData.stress_nominal = stress_nominal;
outputData.S_cc = S_cc;
outputData.ratio = ratio;



