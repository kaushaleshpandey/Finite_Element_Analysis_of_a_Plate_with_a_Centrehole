%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 INDIAN INSTITUTE OF TECHNOLOGY GUWAHATI                 %
%                  DEPARTMENT OF MECHANICAL ENGINEERING                   %
%                                                                         %
%                          2022-23 2ND SEMESTER                           %
%                                                                         %
%               ME 682 - NONLINEAR FINITE ELEMENT METHODS                 %
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%The code is written for solving a finding the displacement, strains and  % 
% stresses for a plate with a hole subjected to biaxial tractions .       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalization Parameters
Lo = 1 ;  % Length
Eo = 1 ; % Youngs Modulous
To = 1 ;  % Time

filename = 'Input/normalization_parameters.txt' ;
fidNP = fopen(filename,'w') ;
fprintf(fidNP,'%25.15f \t %25.15f  \t %25.15f ',Lo,Eo,To);
fclose(fidNP);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fix the origin of the problem

Xo = 0   ; 
Yo = 0   ;

filename = 'Input/origin.txt' ;
fid = fopen(filename,'w') ;
fprintf(fid,'%20.15f \t %20.15f',Xo,Yo);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Dimension of the geometry; Here it is a Plate
% 
% Lx = Length of the Plate in the x direction
% Ly = Length of the Plate in y direction
% thickness_of_beam = this is thickness of the beam in the z direction.
% Note that this has to be changed for plane stress cases 

% Lx = 10 ;
% Ly = 1 ;
% changing this as per plate with hole problem (0.8m * 0.4m)
 Lx = 0.4 ;
 Ly = 0.4 ;

thickness_of_beam =  0.08 ; % for plane strain thickness is set to 1.

filename = 'Input/geometric_data.txt' ;
fid = fopen(filename,'w') ;
fprintf(fid,'%20.15f \t %20.15f  \t %20.15f',Lx,Ly,thickness_of_beam);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Material properties of the Plate
%
% E = Young's modulus of the Plate
% nu = Poissions ratio of the Plate

 Eo = 210e09 ;
 nu = 0.3 ;


%
% option_type_2D --> plane strain
%
% option_type_2D --> plane strain

option_type_2D = 1 ;

get_material_data

filename = 'Input/material_data.txt' ;
fid = fopen(filename,'w') ;
fprintf(fid,'%20.15f \t %20.15f \t %g',Eo,nu,option_type_2D);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Finite element data

get_connectivity_coordinate_data

% Displacement rearragement arrays - this is needed at the end when the
% displacement is added to the current coordinates to get the deformed
% coordinates
i1 = zeros(nno,1) ; i2 = zeros(nno,1) ;
for i = 1:nno
    i1(i) = 2*i-1 ; % x dof of all the nodes
    i2(i) = 2*i ; % y dofs of all the nodes
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get the essential boundary condition array
%

get_boundary_conditon_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% External Force Vector - this has to be computed for point load and
% distributed load
F_ext = zeros(ndof,1) ;  

% % Point load
% F_tip = 2000 ;


% Body Forces - these take care of gravity load.
ldv = 1 ; % Flag to enable gravity load
g = [0 ; 0] ; % specify acc due to gravity value with direction here.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial displacements. U is total displacement whereas DU is incremental
% displacement
U = zeros(ndof,1) ; DU = zeros(ndof,1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quadature parameters for volume intergation
jgaus = 2 ; fequad ; xigv = xig2 ; ngpv = jgaus^2 ; 

% Quadature parameters for surface intergation. This is needed when the
% surface integrals in cases like UDL are required.
jgaus = 5 ; fequad ; xigs = xig1 ; ngps = jgaus ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Call the function to find the size of arrays needed for parallel
% computing so that code can used for meshes with extremely large number of
% elements.

parallel_computing_array


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Factor by which the displacement must be multiplied to get a visible
% deformed geometry at the end of the deformation process
factor = 1e6 ;

filename = 'Input/factor.txt' ;
fid = fopen(filename,'w') ;
fprintf(fid,'%g',factor);
fclose(fid);

% Initialize the total time of the analysis (only for dynamic problems)
ttime = 0 ;

 % FLag for error catching. 0 means no error. It will be set to nonzero 
 % value in case an error is detected
jterm = 0 ;


% Computation Options
jinfo = 2 ;
jcomp = 1 ;