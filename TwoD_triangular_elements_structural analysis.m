function TwoD_triangular_elements()
%{
This is the main program, where a 2D triangular element model can be set,
starting at Line 114.

Test case statement: 
Rectangular plate (L = 6 in.; H = 1 in; t = 1 in.) fixed on the left and
loaded in right-ward horizontal (positive x) direction at the free tip
by two 10,000-lbf forces (one at each node).
Young's modulus = 30*10^6 psi.
Poisson'ratio = 0.3.
Find displacements at all nodes & stresses in all elements for plane strains.

                  y       |\ |\ |\ |\ |
Schematic:        |_ x    | \| \| \| \|

% Test case solution:
% U = 
%    (1,1)       0.0000
%    (2,1)       0.0000
%    (3,1)       0.0007
%    (4,1)       0.0000
%    (5,1)       0.0017
%    (6,1)      -0.0003
%    (7,1)       0.0026
%    (8,1)      -0.0005
%    (9,1)       0.0035
%   (10,1)      -0.0008
%   (11,1)       0.0000
%   (12,1)      -0.0000
%   (13,1)       0.0009
%   (14,1)      -0.0003
%   (15,1)       0.0018
%   (16,1)      -0.0005
%   (17,1)       0.0027
%   (18,1)      -0.0008
%   (19,1)       0.0036
%   (20,1)      -0.0010
          
%}

% To run the test case, simply add a "%" at the beginning of the line below
%%{
clc
clear;
ndof_per_node = 2;      % 2 local degrees of freedom: u, v
nnode_per_elem = 3;     % 3 nodes in triangular element

% User input required in this section

% Building 2D finite element model
L = 7;                  % Length (in)
H = 1;                  % Height (in)
t = 0.02;                  % Thickness (in)
E = 2.07*10^11;            % Young's modulus of Element (psi)
nu = 0.3;               % Poisson's ratio (-)
F = 125000/H;              % Point force (lbf)

nnodeL = 6;             % number of seed nodes in x-direction
nnodeH = 2;             % number of seed nodes in y-direction
nnode = nnodeL*nnodeH;  % number of nodes in model
nelem = 2*(nnodeL-1)*(nnodeH-1);% number of 2D triangular elements in model

% Nodal coordinates
Nodal_coordinates = zeros(nnode,2);
node_index = 0;
for node_indexH = 1:nnodeH
    for node_indexL = 1:nnodeL
        node_index = node_index + 1;
        Nodal_coordinates(node_index,1) = (node_indexL-1)*L/(nnodeL-1);
        Nodal_coordinates(node_index,2) = (node_indexH-1)*H/(nnodeH-1);
    end
end

% Connectivity table + section and material properties per element
Connect_table = zeros(nelem,6);
elem_index = 0;
for node_indexH = 1:nnodeH - 1
    for node_indexL = 1:nnodeL - 1
        elem_index = elem_index + 1;
    Connect_table(elem_index,1) = node_indexL + (node_indexH-1)*nnodeL;
    Connect_table(elem_index,2) = node_indexL + 1 + (node_indexH-1)*nnodeL;
    Connect_table(elem_index,3) = node_indexL + (node_indexH)*nnodeL;      
        elem_index = elem_index + 1;
    Connect_table(elem_index,1) = node_indexL + (node_indexH)*nnodeL;   
    Connect_table(elem_index,2) = node_indexL + 1 + (node_indexH-1)*nnodeL;
    Connect_table(elem_index,3) = node_indexL + 1 + (node_indexH)*nnodeL; 
    end
end

for elem_index = 1:nelem
    Connect_table(elem_index,4) = E;        % Young's modulus of element
    Connect_table(elem_index,5) = nu;       % Poisson's ratio of element
    Connect_table(elem_index,6) = t;        % thickness of element
end

% Known external loads
Fa = zeros(nnode*ndof_per_node,1);
Fa((nnodeL-1)*2+1,1) = F;
Fa((nnodeL*2-1)*2+1,1) = -F;

% Boundary conditions (1: dof is fixed)
BC = zeros(nnode,ndof_per_node);

for node_index = 1:nnodeL:nnode
    node_number_with_BC = node_index;
    BC(node_number_with_BC,1) = 1;
    BC(node_number_with_BC,2) = 1;
end

%}

% -- User input required in section below 
% (no need to re-invent the wheel: copy and edit code from above as needed)



% -- End of user input required; leave the rest alone

% Global element stiffness matrices
Ke = Global_element_stiffness_matrices(nelem,Nodal_coordinates,Connect_table);

% Assemblage matrix
Ka = Assemblage_stiffness_matrix(nelem,nnode,nnode_per_elem,ndof_per_node,...
                                 Connect_table,Ke);

% Modifying Ka into KaSol and Fa into FaSol to get solution 
coef = max(max(Ka))*10E19;
KaSol = Ka;
FaSol = Fa;
for node_index = 1:nnode
    for dof_index = 1:ndof_per_node
        if BC(node_index,dof_index) == 1
            KaSol(dof_index + (node_index-1)*ndof_per_node,...
                  dof_index + (node_index-1)*ndof_per_node) = ...
                  Ka(dof_index,dof_index) + coef;
            FaSol(dof_index + (node_index-1)*ndof_per_node,1) = 0;
        end
    end
end

% Solution without post-processing
U = sparse(KaSol)\sparse(FaSol);  
disp(U);

Sig = Stresses(nelem,Nodal_coordinates,Connect_table,U);
disp(Sig);

% Representation of undeformed and deformed shapes
Undeformed = triangulation(Connect_table(:,1:3),Nodal_coordinates);
Updated_nodal_coordinates = zeros(nnode,2);

Scale_factor = 50;
for elem_index = 1:nelem     
 index = Connect_table(elem_index,1);
 Updated_nodal_coordinates(index,1) = ...
     Nodal_coordinates(index,1) + Scale_factor*U((index-1)*2 + 1);
 Updated_nodal_coordinates(index,2) = ...
     Nodal_coordinates(index,2) + Scale_factor*U(index*2);
 
 index = Connect_table(elem_index,2);
 Updated_nodal_coordinates(index,1) = ...
     Nodal_coordinates(index,1) + Scale_factor*U((index-1)*2 + 1);
 Updated_nodal_coordinates(index,2) = ...
     Nodal_coordinates(index,2) + Scale_factor*U(index*2);   
 
 index = Connect_table(elem_index,3);
 Updated_nodal_coordinates(index,1) = ...
     Nodal_coordinates(index,1) + Scale_factor*U((index-1)*2 + 1);
 Updated_nodal_coordinates(index,2) = ...
     Nodal_coordinates(index,2) + Scale_factor*U(index*2);    
end 

Deformed = triangulation(Connect_table(:,1:3),Updated_nodal_coordinates);
figure
triplot(Undeformed,"-+b"); hold on;
triplot(Deformed,"-+r"); hold on;
axis equal

end

function Ke = Global_element_stiffness_matrices(nelem,Nodal_coordinates,Connect_table)
   Ke = zeros(6,6,nelem);
   for elem_index = 1:nelem     
        index = Connect_table(elem_index,1);
        xi = Nodal_coordinates(index,1);
        yi = Nodal_coordinates(index,2);
        index = Connect_table(elem_index,2);
        xj = Nodal_coordinates(index,1);
        yj = Nodal_coordinates(index,2);
        index = Connect_table(elem_index,3);
        xk = Nodal_coordinates(index,1);
        yk = Nodal_coordinates(index,2);   
        E = Connect_table(elem_index,4);
        nu = Connect_table(elem_index,5);
        t = Connect_table(elem_index,6);
        
        % Course notes, 4.3.1
        A = 0.5*det([1 xi yi;...
                     1 xj yj;...
                     1 xk yk]);
        
        m21 = (yj-yk)/(2*A);
        m31 = (xk-xj)/(2*A);
        
        m22 = (yk-yi)/(2*A);
        m32 = (xi-xk)/(2*A);
        
        m23 = (yi-yj)/(2*A);
        m33 = (xj-xi)/(2*A);

        % Course notes, 5.1.2
        B = [m21 0 m22 0 m23 0;...
             0 m31 0 m32 0 m33;...
             m31 m21 m32 m22 m33 m23];
        
         % Course notes, 5.1.3 Plane stress state
         D = E*[1 nu 0;...
              nu 1 0;...
              0 0 (1-nu)/2]/(1-nu^2);
         
%        % Course notes, 5.2 Plane strain state
%        D = E*[1-nu nu 0;...
%              nu 1-nu 0;...
%              0 0 (1-2*nu)/2]/((1+nu)*(1-2*nu));         
        
        % Course notes, 5.1.5
         Ke(:,:,elem_index) = B'*D*B*t*A;
            
    end   

end

function Ka = Assemblage_stiffness_matrix(nelem,nnode,nnode_per_elem,ndof_per_node,...
                                          Connect_table,Ke)
    Ka = zeros(nnode*ndof_per_node,nnode*ndof_per_node);
    for k = 1:nelem 
        for ii = 1:nnode
            for jj = 1:nnode           
                if ii == Connect_table(k,1) && jj == Connect_table(k,1)
                    for i = 1:ndof_per_node
                        for j = 1:ndof_per_node
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) = ...
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) + Ke(i,j,k);                       
                        end
                    end
                elseif ii == Connect_table(k,2) && jj == Connect_table(k,2)
                    for i = 1:ndof_per_node
                        for j = 1:ndof_per_node
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) = ...
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) + ...
    Ke(i + ndof_per_node,j + ndof_per_node,k);
                        end
                    end
                elseif ii == Connect_table(k,3) && jj == Connect_table(k,3)
                    for i = 1:ndof_per_node
                        for j = 1:ndof_per_node
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) = ...
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) + ...
    Ke(i + (nnode_per_elem-1)*ndof_per_node,j + (nnode_per_elem-1)*ndof_per_node,k);
                        end
                    end                    
                elseif ii == Connect_table(k,1) && jj == Connect_table(k,2)
                    for i = 1:ndof_per_node
                        for j = 1:ndof_per_node
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) = ...
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) + ...
    Ke(i,j + ndof_per_node,k);
                        end
                    end
                elseif ii == Connect_table(k,2) && jj == Connect_table(k,1)
                    for i = 1:ndof_per_node
                        for j = 1:ndof_per_node 
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) = ...
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) + ...
    Ke(i + ndof_per_node,j,k);
                        end
                    end                    
                elseif ii == Connect_table(k,1) && jj == Connect_table(k,3)
                    for i = 1:ndof_per_node
                        for j = 1:ndof_per_node
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) = ...
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) + ...
    Ke(i,j + (nnode_per_elem-1)*ndof_per_node,k);
                        end
                    end                    
                elseif ii == Connect_table(k,3) && jj == Connect_table(k,1)
                    for i = 1:ndof_per_node
                        for j = 1:ndof_per_node
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) = ...
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) + ...
    Ke(i + (nnode_per_elem-1)*ndof_per_node,j,k);
                        end
                    end            
                elseif ii == Connect_table(k,2) && jj == Connect_table(k,3)
                    for i = 1:ndof_per_node
                        for j = 1:ndof_per_node 
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) = ...
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) + ...
    Ke(i + ndof_per_node,j + (nnode_per_elem-1)*ndof_per_node,k);
                        end
                    end         
                elseif ii == Connect_table(k,3) && jj == Connect_table(k,2)
                    for i = 1:ndof_per_node
                        for j = 1:ndof_per_node
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) = ...
    Ka((ii - 1)*ndof_per_node + i,(jj - 1)*ndof_per_node + j) + ...
    Ke(i + (nnode_per_elem-1)*ndof_per_node,j + ndof_per_node,k);
                        end
                    end
                end                
            end            
        end
    end
end
function Sig = Stresses(nelem,Nodal_coordinates,Connect_table,U)
   Sig = zeros(3,nelem);
   B = zeros(3,6,nelem);
   ui = zeros(1,nelem);
   vi = zeros(1,nelem);
   uj = zeros(1,nelem);
   vj = zeros(1,nelem);
   uk = zeros(1,nelem);
   vk = zeros(1,nelem);   
   for elem_index = 1:nelem     
        index = Connect_table(elem_index,1);
        xi = Nodal_coordinates(index,1);
        yi = Nodal_coordinates(index,2);
        ui(1,elem_index) = U((index-1)*2 + 1);
        vi(1,elem_index) = U(index*2); 
        
        index = Connect_table(elem_index,2);
        xj = Nodal_coordinates(index,1);
        yj = Nodal_coordinates(index,2);
        uj(1,elem_index) = U((index-1)*2 + 1);
        vj(1,elem_index) = U(index*2); 
        
        index = Connect_table(elem_index,3);
        xk = Nodal_coordinates(index,1);
        yk = Nodal_coordinates(index,2);   
        E = Connect_table(elem_index,4);
        uk(1,elem_index) = U((index-1)*2 + 1);
        vk(1,elem_index) = U(index*2); 
        
        nu = Connect_table(elem_index,5);
                
        % Course notes, page 37
        A = 0.5*det([1 xi yi;...
                     1 xj yj;...
                     1 xk yk]);
        
        m21 = (yj-yk)/(2*A);
        m31 = (xk-xj)/(2*A);
        
        m22 = (yk-yi)/(2*A);
        m32 = (xi-xk)/(2*A);
        
        m23 = (yi-yj)/(2*A);
        m33 = (xj-xi)/(2*A);

        % Course notes, page 47
        B(:,:,elem_index) = [m21 0 m22 0 m23 0;...
                             0 m31 0 m32 0 m33;...
                             m31 m21 m32 m22 m33 m23];     
            
   end   
  % Course notes, 5.1.3 Plane stress state
   D = E*[1 nu 0;...
          nu 1 0;...
          0 0 (1-nu)/2]/(1-nu^2);
         
%   % Course notes, 5.2 Plane strain state
%    D = E*[1-nu nu 0;...
%           nu 1-nu 0;...
%           0 0 (1-2*nu)/2]/((1+nu)*(1-2*nu));    
       
   % Course notes, page 53
   for elem_index = 1:nelem
       Ue = [ui(1,elem_index);vi(1,elem_index);...
             uj(1,elem_index);vj(1,elem_index);...
             uk(1,elem_index);vk(1,elem_index)];
       Sig(:,elem_index) = D*B(:,:,elem_index)*Ue;  
   end
    
end
