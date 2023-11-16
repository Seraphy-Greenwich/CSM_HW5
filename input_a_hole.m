function [inputData]=input_a_hole(const)

% Meshing
D = 20;             %[mm]
d = D/(2*const+1);  %[mm]
L = 56;             %[mm]
F_total = 500;      %[N]
rho = const*d;
initalData = [ D d rho F_total];

t = pi:-pi/12:0;
pgon = polyshape({[-L/2 -L/2 rho*cos(t) L/2 L/2]}, ...
                 {[D/2 0 rho*sin(t) 0 D/2]});
tr = triangulation(pgon);
model = createpde;
tnodes = tr.Points';
telements = tr.ConnectivityList';
geometryFromMesh(model,tnodes,telements);
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.4,'Hmin',0.1);
pdemesh(model)
hold on

% Set model data
nodeCoordinate = mesh.Nodes';   %[mm]
nodeConnectivity = mesh.Elements';  
numberofNode = length(nodeCoordinate);                   % number of node
numberofElement = length(nodeConnectivity);              % number of element
DOFperNode = 2;                                          % DOF per node
totalDOF = DOFperNode*numberofNode;                      % total DOF

E(1:numberofElement,:) = 7e4;   %[N/mm^2]
nu(1:numberofElement,:) = 0.25;
materialData = [E nu];

% Find boundary and loading condition
boundary_nodes = findNodes(mesh,'region','Edge',[2 5])';
boundaryRigid_nodes = findNodes(mesh,'box',[-rho 0],[-0.1 0.1])';
loadingR_nodes = findNodes(mesh,'region','Edge',[1])';
loadingHalfR_nodes = [findNodes(mesh,'nearest',[L/2;0]); findNodes(mesh,'nearest',[L/2;D/2])];
loadingL_nodes = findNodes(mesh,'region','Edge',[4])';
loadingHalfL_nodes = [findNodes(mesh,'nearest',[-L/2;0]); findNodes(mesh,'nearest',[-L/2;D/2])];
central_nodes = findNodes(mesh,'box',[-5 5], [-D/2 D/2])';
central_Elements = findElements(mesh,'box',[-5 5], [-D/2 D/2])';


[C,ia,ib] = intersect(boundary_nodes,loadingHalfR_nodes);
boundary_nodes(ia) = [];
[C,ia,ib] = intersect(loadingR_nodes,loadingHalfR_nodes);
loadingR_nodes(ia) = [];
[C,ia,ib] = intersect(boundary_nodes,loadingHalfL_nodes);
boundary_nodes(ia) = [];
[C,ia,ib] = intersect(loadingL_nodes,loadingHalfL_nodes);
loadingL_nodes(ia) = [];


boundary_nodes_idx = [];
for i = 1:length(boundary_nodes)
    node = boundary_nodes(i);
    boundary_nodes_idx = [boundary_nodes_idx DOFperNode*(node-1)+1 : DOFperNode*node];
end
boundaryRigid_nodes_idx = [];
for i = 1:length(boundaryRigid_nodes)
    node = boundaryRigid_nodes(i);
    boundaryRigid_nodes_idx = [boundaryRigid_nodes_idx DOFperNode*(node-1)+1 : DOFperNode*node];
end


loadingR_nodes_idx = [];
for i = 1:length(loadingR_nodes)
    node = loadingR_nodes(i);
    loadingR_nodes_idx = [loadingR_nodes_idx DOFperNode*(node-1)+1 : DOFperNode*node];
end
loadingHalfR_nodes_idx = [];
for i = 1:length(loadingHalfR_nodes)
    node = loadingHalfR_nodes(i);
    loadingHalfR_nodes_idx = [loadingHalfR_nodes_idx DOFperNode*(node-1)+1 : DOFperNode*node];
end


loadingL_nodes_idx = [];
for i = 1:length(loadingL_nodes)
    node = loadingL_nodes(i);
    loadingL_nodes_idx = [loadingL_nodes_idx DOFperNode*(node-1)+1 : DOFperNode*node];
end
loadingHalfL_nodes_idx = [];
for i = 1:length(loadingHalfL_nodes)
    node = loadingHalfL_nodes(i);
    loadingHalfL_nodes_idx = [loadingHalfL_nodes_idx DOFperNode*(node-1)+1 : DOFperNode*node];
end


% Draw the plot
% figure
% pdemesh(model,'NodeLabels','on')
% pdemesh(model,'ElementLabels','on')
hold on
plot(mesh.Nodes(1,loadingR_nodes),mesh.Nodes(2,loadingR_nodes),'>','MarkerFaceColor','r')
plot(mesh.Nodes(1,loadingL_nodes),mesh.Nodes(2,loadingL_nodes),'<','MarkerFaceColor','r')
plot(mesh.Nodes(1,loadingHalfR_nodes),mesh.Nodes(2,loadingHalfR_nodes),'>','MarkerFaceColor','b')
plot(mesh.Nodes(1,loadingHalfL_nodes),mesh.Nodes(2,loadingHalfL_nodes),'<','MarkerFaceColor','b')
plot(mesh.Nodes(1,central_nodes),mesh.Nodes(2,central_nodes),'or','MarkerFaceColor','y')
plot(mesh.Nodes(1,boundary_nodes),mesh.Nodes(2,boundary_nodes),'or','MarkerFaceColor','g')
plot(mesh.Nodes(1,boundaryRigid_nodes),mesh.Nodes(2,boundaryRigid_nodes),'or','MarkerFaceColor','k')
axis equal;

title(['Configuration (rho/d: ', num2str(const),')']);

% Set boundary and loading condition
constrainedDispl = zeros(totalDOF,1);
idx1 = [1:2:length(boundary_nodes_idx)];
boundary_nodes_idx(idx1) = [];
constrainedDispl_globalidx = [boundary_nodes_idx boundaryRigid_nodes_idx(1)];       % global index having constrained displacement

constrainedForce = zeros(totalDOF,1);
% Right part
fR = F_total/(2*(length(loadingR_nodes)+length(loadingHalfR_nodes)-1));
idx2 = [2:2:length(loadingR_nodes_idx)];
loadingR_nodes_idx(idx2) = [];
idx3 = [2:2:length(loadingHalfR_nodes_idx)];
loadingHalfR_nodes_idx(idx3) = [];
constrainedForce(loadingR_nodes_idx) = fR;           %[N]
constrainedForce(loadingHalfR_nodes_idx) = fR/2;     %[N]
% Left part
fL = F_total/(2*(length(loadingL_nodes)+length(loadingHalfL_nodes)-1));
idx4 = [2:2:length(loadingL_nodes_idx)];
loadingL_nodes_idx(idx2) = [];
idx5 = [2:2:length(loadingHalfL_nodes_idx)];
loadingHalfL_nodes_idx(idx3) = [];
constrainedForce(loadingL_nodes_idx) = -fL;           %[N]
constrainedForce(loadingHalfL_nodes_idx) = -fL/2;     %[N]

% Save input data
inputData.numberofNode = numberofNode;
inputData.numberofElement = numberofElement;
inputData.DOFperNode = DOFperNode;
inputData.totalDOF = totalDOF;
inputData.nodeCoordinate = nodeCoordinate;
inputData.nodeConnectivity = nodeConnectivity;
inputData.constrainedDispl = constrainedDispl;
inputData.constrainedDispl_globalidx = constrainedDispl_globalidx;
inputData.constrainedForce = constrainedForce;
inputData.materialData = materialData;
inputData.initalData = initalData;
inputData.const = const;
inputData.central_Elements = central_Elements;