# Finite-Element-Analysis
##In this project you will see the structural analysis of the beamm element using ANSYS APDL and MATLAB for plane183 and beam elements.

                         y         |\ |\ |\ |\ |
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
L = 7;                  % Length (m)
H = 1;                  % Height (m)
t = 0.02;                  % Thickness (m)
E = 2.07*10^11;            % Young's modulus of Element (Pa)
nu = 0.3;               % Poisson's ratio (-)
F = 125000/H;              % Point force (N)
% we applied our parameters in this code

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
% In our problem we have two forces in the opposite direction, so we refined the direction of the forces

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
%   we need to use plane stress, so we activated the above-highlighted lines       
         
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
%   we need to use plane stress, so we activated the above-highlighted lines
         
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

 
MATLAB Results for Deliverable 2:
BLUE HIGHLIGHTED ARE THE FINAL RESULTS 
Y-Displacement
   (1,1)       0.0000
   (2,1)       0.0000
   (3,1)       0.0000
   (4,1)       0.0001
   (5,1)       0.0001
   (6,1)       0.0003
   (7,1)       0.0001
   (8,1)       0.0006
   (9,1)       0.0002
  (10,1)       0.0010
  (11,1)       0.0002
  (12,1)       0.0016
  (13,1)      -0.0000
  (14,1)      -0.0000
  (15,1)      -0.0000
  (16,1)       0.0001
  (17,1)      -0.0001
  (18,1)       0.0003
  (19,1)      -0.0001
  (20,1)       0.0006
  (21,1)      -0.0002
  (22,1)       0.0010
  (23,1)      -0.0002
  (24,1)       0.0016
Stresess:
   1.0e+06 *
    7.4140   -7.4140    7.4140   -7.4140    7.4140   -7.4140    7.4151   -7.4151    7.6372   -7.6372
    2.2242   -2.2242    2.2242   -2.2242    2.2242   -2.2241    2.2246   -2.2012    2.3145    2.4810
    3.6329   -3.6329    3.6329   -3.6329    3.6329   -3.6329    3.6321   -3.6321    3.4735   -3.4735


Deliverable 3: ANSYS script for the same number of the element as MATLAB
/clear
/PREP7     

!Defining material property               
et,1,plane182,,0                                           
keyopt,1,3,3
r,1,0.02                    
mp,ex,1,2.07e11               
mp,prxy,1,0.3  

!Defining keypoints to build a geometry        
K,1,0,0,0
K,2,0,1,0
K,3,7,1,0
k,4,7,0,0

!create lines using keypoints to create a geometry
lstr,1,2
lstr,2,3
lstr,3,4                  
lstr,4,1

!Creating an area using lines created
al,1,2,3,4 

!Dividing lines and prepare the geometry for meshing                   
lesize,1,,,1
lesize,2,,,5
lesize,3,,,1  
lesize,4,,,5 

!Meshing the geometry  
mshape,1,2D                   
mshkey,1        
amesh,all                     

!Apply boundary conditions by fixing
dk,1,all,0                    
dk,2,all,0

!Apply forces
fk,4,fx,125000        
fk,3,fx,-125000    

finish
/solu
solve                         
/post1
set
/output
prnsol,s 
prnsol,u,y     
ANSYS Result for deliverable 3: Displacement and Stress
BLUE HIGHLIGHTED ARE THE FINAL RESULTS 
Displacement in ANSYS:
  ***** POST1 NODAL DEGREE OF FREEDOM LISTING *****                            
 
  LOAD STEP=     1  SUBSTEP=     1                                             
   TIME=    1.0000      LOAD CASE=   0                                         
 
  THE FOLLOWING DEGREE OF FREEDOM RESULTS ARE IN THE GLOBAL COORDINATE SYSTEM  
 
    NODE       UY    
       1   0.0000     
       2  0.15943E-002
       3  0.63882E-004
       4  0.25553E-003
       5  0.57494E-003
       6  0.10221E-002
       7  0.16174E-002
       8   0.0000     
       9  0.10222E-002
      10  0.57494E-003
      11  0.25553E-003
      12  0.63882E-004

 MAXIMUM ABSOLUTE VALUES
 NODE          7
 VALUE   0.16174E-002
 
Stresses in ANSYS:
PRINT S    NODAL SOLUTION PER NODE
 
  ***** POST1 NODAL STRESS LISTING *****                                       
  PowerGraphics Is Currently Enabled                                           
 
  LOAD STEP=     1  SUBSTEP=     1                                             
   TIME=    1.0000      LOAD CASE=   0                                         
  NODAL RESULTS ARE FOR MATERIAL   1                                           
 
  THE FOLLOWING X,Y,Z VALUES ARE IN GLOBAL COORDINATES                         
 
    NODE     SX           SY           SZ           SXY          SYZ          SXZ     
       1  0.74140E+007 0.22242E+007  0.0000      0.36329E+007  0.0000       0.0000     
       2 -0.46566E-009 0.23978E+007  0.0000     -0.23283E-009  0.0000       0.0000     
       3  0.24713E+007 0.74140E+006  0.0000      0.12110E+007  0.0000       0.0000     
       4  0.24713E+007 0.74140E+006  0.0000      0.12110E+007  0.0000       0.0000     
       5  0.24717E+007 0.74159E+006  0.0000      0.12107E+007  0.0000       0.0000     
       6  0.25457E+007 0.77931E+006  0.0000      0.11578E+007  0.0000       0.0000     
       7 -0.76372E+007 0.24810E+007  0.0000     -0.34735E+007  0.0000       0.0000     
       8   0.0000      0.23283E-009  0.0000      0.23283E-009  0.0000       0.0000     
       9 -0.24717E+007 0.86478E+006  0.0000     -0.12107E+007  0.0000       0.0000     
      10 -0.24713E+007-0.73354E+006  0.0000     -0.12110E+007  0.0000       0.0000     
      11 -0.24713E+007-0.74136E+006  0.0000     -0.12110E+007  0.0000       0.0000     
      12 -0.24713E+007-0.74140E+006  0.0000     -0.12110E+007  0.0000       0.0000     
 MINIMUM VALUES
 NODE          7           12            1            7            1            1
 VALUE  -0.76372E+007-0.74140E+006  0.0000     -0.34735E+007  0.0000       0.0000     

 MAXIMUM VALUES
 NODE          1            7            1            1            1            1
 VALUE   0.74140E+007 0.24810E+007  0.0000      0.36329E+007  0.0000       0.0000

 

 





 


















Conclusion Step 1:

•	In the first part of deliverable 1, we calculated deflection based on the strength of the material and we calculated y_max= 8.87 E-3, and maximum stress= 37.5 MPa.  

Answering questions of deliverable 3:

1-	YES, THE RESULTS ARE THE SAME, we changed the parameters in the provided MATLAB code and the results are as follows: 
Maximum deflection: 0.16174E-002
Maximum stress: 0.74140E+007
In ANSYS we get the exact same results as MATLAB which is different from our hand calculation. We expected the same results because the same mesh types are used in these deliverables. 

2-	In this type of question because of the nature of the problem, hand calculation is accurate. Also, we did not use enough meshes in deliverables 2 and 3, and we expect different results from deliverable 1. 
 
Deliverable 4: Beam188 mesh refinement
Code:
/clear
i = 1
*do,i,1,10
/prep7

!Defining material property
et,1,beam188                                       				!Type of element
keyopt,1,5,1	                                   				!
l=7                                                				!Lenght
mp,ex,1,2.07e11                                    				!Elastic modulus
mp,prxy,1,0.3                                      				!Define nu

!Defining keypoints to use to build a geometry
k,1,0,0                                            
k,2,l,0
k,3,0,0,100

!Creating line using the defined keypoints to build a geometry
lstr,1,2

!Defining cross sectional area of the geometry
sectype,1,beam,rect,,0
secdata,0.02,1,2,2
secnum,1

!Selecting lines, dividing them and meshing
lsel,s,line,,1
latt,1,,1,,3,,1

!Dividing lines i times in respect to our do loop
lesize,all,,,i

!Mesh the geometry
lmesh,all

!Defining degree of freedom and forces applied
dk,1,all
fk,2,Mz,125000 

eplot
allsel

finish
/solu
solve
finish

/post1
set,last
max_y=Uy(2)

!Create a table based on what needed in our problem
etable,sbzt,smisc,34,39

!Get information that we need
*get,max_stress,elem,1,etab,sbzt

!Store information that we need in specified folder and send it to folder as a massage
/output,C:\Users\Daniyel\Desktop\FEA_project\Deliverable4_result,txt,,append  
*msg,info,i,max_stress,max_y
nelem = %I, max_stress = %G, max_y = %G 

/output
finish
/clear

*enddo

Results in Table for Y-Displacement at the tip and Maximum bending Stress at fixed point:
Nelem = 1, max_stress = -37500000, max_y = 8.876811594E-03.             
 Nelem = 2, max_stress = -37500000, max_y = 8.876811594E-03.             
 Nelem = 3, max_stress = -37500000, max_y = 8.876811594E-03.             
 Nelem = 4, max_stress = -37500000, max_y = 8.876811594E-03.             
 Nelem = 5, max_stress = -37500000, max_y = 8.876811594E-03.             
 Nelem = 6, max_stress = -37500000, max_y = 8.876811594E-03.             
 Nelem = 7, max_stress = -37500000, max_y = 8.876811594E-03.             
 Nelem = 8, max_stress = -37500000, max_y = 8.876811594E-03.             
 Nelem = 9, max_stress = -37500000, max_y = 8.876811594E-03.             
 Nelem = 10, max_stress = -37500000, max_y = 8.876811594E-03.







Plot for mesh refinement:

 


How many elements are enough for this simulation?
One element is enough for this simulation because we cannot see any improvement in the mesh refinement process.






Deliverable 5: The same analysis as deliverables 2 and 3 but with quadratic mesh

/PREP7     

!Defining material property               
et,1,plane182,,0                                           
keyopt,1,3,3
r,1,0.02                    
mp,ex,1,2.07e11               
mp,prxy,1,0.3  

!Defining keypoints            
K,1,0,0,0
K,2,0,1,0
K,3,7,1,0
k,4,7,0,0

!create lines using keypoints
lstr,1,2
lstr,2,3
lstr,3,4                  
lstr,4,1

!Creating an area
al,1,2,3,4 

!Dividing lines and prepare them for meshing                   
lesize,1,,,1
lesize,2,,,5
lesize,3,,,1  
lesize,4,,,5

!Meshing the geometry in the quadratic shape
mshape,0,2D                   
mshkey,1        
amesh,all                     

!Apply boundary conditions by fixing
dk,1,all,0                    
dk,2,all,0

!Apply forces
fk,4,fx,125000        
fk,3,fx,-125000    

finish
/solu
solve
finish

/post1
set
/output
prnsol,s 
prnsol,u,y








ANSYS RESULTS, DISPLACEMENT

 PRINT S    NODAL SOLUTION PER NODE
 
  ***** POST1 NODAL STRESS LISTING *****                                       
  PowerGraphics Is Currently Enabled                                           
 
  LOAD STEP=     1  SUBSTEP=     1                                             
   TIME=    1.0000      LOAD CASE=   0                                         
  NODAL RESULTS ARE FOR MATERIAL   1                                           
 
  THE FOLLOWING X,Y,Z VALUES ARE IN GLOBAL COORDINATES                         
 
    NODE     SX           SY           SZ           SXY          SYZ          SXZ     
       1  0.22242E+008 0.66726E+007  0.0000      0.10899E+008  0.0000       0.0000     
       2 -0.22242E+008-0.66726E+007  0.0000      0.10899E+008  0.0000       0.0000     
       3 -0.22242E+008-0.66726E+007  0.0000     -0.10899E+008  0.0000       0.0000     
       4 -0.22242E+008-0.66726E+007  0.0000      0.93132E-009  0.0000       0.0000     
       5 -0.22242E+008-0.66726E+007  0.0000      0.93132E-009  0.0000       0.0000     
       6 -0.22242E+008-0.66726E+007  0.0000      0.93132E-009  0.0000       0.0000     
       7 -0.22242E+008-0.66726E+007  0.0000      0.93132E-009  0.0000       0.0000     
       8  0.22242E+008 0.66726E+007  0.0000     -0.10899E+008  0.0000       0.0000     
       9  0.22242E+008 0.66726E+007  0.0000      0.93132E-009  0.0000       0.0000     
      10  0.22242E+008 0.66726E+007  0.0000      0.93132E-009  0.0000       0.0000     
      11  0.22242E+008 0.66726E+007  0.0000      0.93132E-009  0.0000       0.0000     
      12  0.22242E+008 0.66726E+007  0.0000      0.93132E-009  0.0000       0.0000     

 MINIMUM VALUES
 NODE          2            2            1            3            1            1
 VALUE  -0.22242E+008-0.66726E+007  0.0000     -0.10899E+008  0.0000       0.0000     

 MAXIMUM VALUES
 NODE          1            1            1            1            1            1
 VALUE   0.22242E+008 0.66726E+007  0.0000      0.10899E+008  0.0000       0.0000










ANSYS RESULTS, STRESS


 PRINT U    NODAL SOLUTION PER NODE
 
  ***** POST1 NODAL DEGREE OF FREEDOM LISTING *****                            
 
  LOAD STEP=     1  SUBSTEP=     1                                             
   TIME=    1.0000      LOAD CASE=   0                                         
 
  THE FOLLOWING DEGREE OF FREEDOM RESULTS ARE IN THE GLOBAL COORDINATE SYSTEM  
 
    NODE       UY    
       1   0.0000     
       2   0.0000     
       3  0.47912E-002
       4  0.19165E-003
       5  0.76659E-003
       6  0.17248E-002
       7  0.30663E-002
       8  0.47912E-002
       9  0.30663E-002
      10  0.17248E-002
      11  0.76659E-003
      12  0.19165E-003

 MAXIMUM ABSOLUTE VALUES
 NODE          3
 VALUE   0.47912E-002












Deliverable 5: Plane182 with *do loop and mesh refinement
Code:
/clear
i = 1
*do,i,1,16,1
/PREP7     

!Defining material property               
et,1,plane182,,0                                           
keyopt,1,3,3
r,1,0.02                    
mp,ex,1,2.07e11               
mp,prxy,1,0.3  

!Defining keypoints            
K,1,0,0,0
K,2,0,1,0
K,3,7,1,0
k,4,7,0,0

!create lines using keypoints
lstr,1,2
lstr,2,3
lstr,3,4                  
lstr,4,1

!Creating an area
al,1,2,3,4

!Dividing lines and prepare them for meshing                   
lesize,1,,,i
lesize,2,,,7*i
lesize,3,,,i  
lesize,4,,,7*i 

!Meshing the geometry  
mshape,0,2D                   
mshkey,1        
amesh,all                     

!Apply boundary conditions by fixing
!dk,1,all,0                    
!dk,2,all,0
dl,1,1,all

!Apply forces
fk,4,fx,125000        
fk,3,fx,-125000    

finish
/solu
solve
finish

/post1
set,last
max_y=Uy(i+2)

etable,s:x,s,x
*get,max_stress,elem,1,etab,s:x
/output,C:\Users\Daniyel\Desktop\FEA_project\Deliverable5_result,txt,,append  
*msg,info,i,max_stress,max_y
nelem = %I, max_stress = %G, max_y = %G 

/output
finish
/clear

*enddo


Results and Table for Y-displacement and maximum bending stress:
Nelem = 1, max_stress = 0, max_y = 5.983628556E-03.                     
 Nelem = 2, max_stress = 16938047.6, max_y = 7.901850752E-03.            
 Nelem = 3, max_stress = 23831041.5, max_y = 8.421346558E-03.            
 Nelem = 4, max_stress = 27435187.2, max_y = 8.626614892E-03.            
 Nelem = 5, max_stress = 29659795, max_y = 8.729092034E-03.              
 Nelem = 6, max_stress = 31204860.5, max_y = 8.788448479E-03.            
 Nelem = 7, max_stress = 32367744.5, max_y = 8.826529173E-03.            
 Nelem = 8, max_stress = 33295607, max_y = 8.852846592E-03.              
 Nelem = 9, max_stress = 34068459, max_y = 8.872089984E-03.              
 Nelem = 10, max_stress = 34733363.5, max_y = 8.886794977E-03.           
 Nelem = 11, max_stress = 35319715, max_y = 8.898434851E-03.             
 Nelem = 12, max_stress = 35846789, max_y = 8.90791627E-03.              
 Nelem = 13, max_stress = 36327755.5, max_y = 8.91582409E-03.            
 Nelem = 14, max_stress = 36771923.5, max_y = 8.922550691E-03.           
 Nelem = 15, max_stress = 37186065.5, max_y = 8.928368228E-03.           
 Nelem = 16, max_stress = 37575241.5, max_y = 8.933470938E-03.
Nelem = 17, max_stress = 37943310, max_y = 8.938000931E-03.             
 Nelem = 18, max_stress = 38293277.5, max_y = 8.942064461E-03.           
 Nelem = 19, max_stress = 38627521, max_y = 8.945742514E-03.             
 Nelem = 20, max_stress = 38947961.5, max_y = 8.949097879E-03.







Plot for mesh refinement in deliverable 5:
 


How many elements are enough for this simulation? (mesh sensitivity analysis)
We think that iterations that will be 500-600 elements would be enough which is in i = 12 because both displacement and stress are converging constant points.

Step 2_Conclusion:

•	In Deliverable 4, our result in the plots shows, stresses at a fixed point and deformations at the tip are exactly the same as hand calculation in deliverable 1. Because the right answer is calculated with one mesh, in mesh refinement we do not have any improvement.
•	In Deliverable 5, as shown in the plot, stress at a fixed point, and deflection at the tip converged to our hand calculation result.
•	We reached the conclusion that our hand calculation and deliverables 4 and 5 are accurate but deliverables 2 and 3 are not reliable and the answer in deliverables 2 and 3 cannot be correct. 
Deliverable 6: Bracket Optimization
Table:
D1	              Von Mises, Mpa	       Mass, kg
55 mm	657.539	 3.63806
60 mm	659.666	 3.61499
65 mm	665.015	 3.59185
70 mm	664.268	       3.5686
75 mm	668.537	 3.54539
		
Chart:

 





Deliverable 7: Bracket Geometry Optimization
First, we checked each of the D properties to have initial guesses for each geometry parameter. Then at the end, we combined them all together for the final optimal geometry.
 
Optimal D1 is between 100 to 110 because mass is descending and Von Misses is not exceeding 800MPa
 
Optimal D2 is the minimum amount because mass is minimum at this point and Von Misses is not exceeding 800MPa

 
Optimal D4 is the maximum amount because mass is descending and Von Misses is not exceeding 800MPa

 
Optimal D6 is between 27 to 30 because mass is low at this point and Von Misses is not exceeding 800MPa


 
Optimal D8 is between 53 to 55 because mass is low at this point and Von Misses is not exceeding 800MPa
 
Optimal D10 is the minimum amount because mass is minimum at this point and Von Misses is not exceeding 800MPa

So, based on the optimization of each D we can combine it and get the final optimal geometry:
The optimal geometry which is less than 800MPa as follows:

D1	D2, D3	D4, D5	D6, D7	D8, D9	D10, D11	Von Misses	Mass
107	10	29	26	53	2	795.249MPa	2.1865kg





 






Changing the element size and reporting the results:

Is the optimization correct if we change the element size by 2?
•	No it will not be the same, by changing the element size of the first optimization in deliverable 7, Von Misses stress exceeds 800MPa and it will be 943.9MPa
To have a geometry that works with both element sizes, we must change the D values as follows:

D1	D2, D3	D4, D5	D6, D7	D8, D9	D10, D11	Von Misses	Von Misses, element divided by 2	Mass
104	10	29	30	55	2	694.383	796.823MPa	2.3917kg


when the element size is divided by two in the second optimization (to be less than 800MPa in element sizes (elsize) 2 and 1) Von Misses stress is 796.823MPa and less than 800MPa

 


In the new geometry when the size of the element is not divided by 2, Von Misses stress is 694.383MPa and less than 800MPa
 















