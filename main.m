% Spring 2023 Computer Simulation and Machine Learning (Prof. Do-Nyun Kim)
% Simulation-driven Structure Design Lab.

% TA session #1

clear; close all; clc;
fprintf('%s\n','%%%%%%%%%% Start %%%%%%%%%%');
fprintf('\n');

C = [0.1:0.1:1];
for ii = 1:length(C)
%% Set input data
const = C(ii);
[inputData]=input_a_hole(const);

%% Get element stiffness data
[elementData] = getElementStiffness(inputData);

%% Assemble global stiffness matrix
[globalKmatrix] = assemGlobalStiffness(elementData, inputData);

%% Reduce matrix and solve the equation
[outputData] = solver(globalKmatrix, elementData, inputData);

totalData(ii,:) = outputData;
saveas(gcf,sprintf('%.1f %s',const,'figure.fig'));
close all;
end
save('totalData_a_hole.mat','totalData');

fprintf('\n');
fprintf('%s\n','%%%%%%%%%% End %%%%%%%%%%');
