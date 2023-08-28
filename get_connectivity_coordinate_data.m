

% % Number of elements in y direction
ny = 24 ;
% 
% % Number of elements in x direction
nx = ny; % <--- change this for your own problem
% 
% % Number of Elements
nel = nx*ny ;
% 
% % Number of Nodes
nno = (nx+1)*(ny+1) ;
% 
% % Number of dofs per element 
ndoel = 8 ;
% 
% % Number of dofs in the FE mesh, Number of dofs per element
 ndof = 2*nno ; ndoelo = 8 ;
% 
% filename = 'Input/fe_data.txt' ;
% fid = fopen(filename,'w') ;
% fprintf(fid,'%g \t %g \t %g \t %g \t %g',ny,nx,nel,ndoel,ndof);
% fclose(fid);

% Element Connectivity
%
CON = zeros(nel,4) ; % since 4 nodes per element
iel = 0 ;
for j = 1:ny
    for i = 1:nx 
        iel = iel + 1 ; 
        CON(iel,:) = [ (j-1)*(nx+1)+i  (j-1)*(nx+1)+i+1  j*(nx+1)+i+1  j*(nx+1)+i ] ;    
    end
end

% filename = 'Input/connectivity.txt' ;
% fid = fopen(filename,'w') ;
% for iel = 1:nel
%     if iel < nel
%         fprintf(fid,'%g \t %g \t %g \t %g \n',CON(iel,:));   
%     else
%         fprintf(fid,'%g \t %g \t %g \t %g',CON(iel,:));
%     end
% end
% fclose(fid);


% Initial Node Coordinates
 Xn = zeros(nno,2) ; % 2 dof per node i.e. x coordinate and y coordinate.

% Providing Initial Node Coordinates
np = ny/2 ;  % no. of angular divisions on each half of the plate
n = nx ;     % no. of divisions along x
Lx = 0.4 ; Ly = 0.4 ; %  dimensions of one quarter of the plate
r = 0.03 ; % Radius of circular hole
theta = atan(Ly/Lx);  
delta1 = theta/np;   
delta2 = ((pi/2) - theta)/np; 

% Code to get the coordinates

j = 0;
while(j*delta1 <= theta)
    h1(j+1) = ((Lx/cos(j*delta1)) - r)/n;
    for i = 0:n
        Xn(j*(n+1)+i+1,1) = (r+i*h1(j+1))*cos(j*delta1);
        Xn(j*(n+1)+i+1,2) = (r+i*h1(j+1))*sin(j*delta1);
    end
    j = j+1;
end

j = 1;
while(j*delta2 > 0 && j*delta2 <= (pi/2 - theta))
    h2(j) = ((Ly/cos(pi/2 - j*delta2 - theta)) - r)/n;
    for i = 0:n
        Xn((j+np)*(n+1)+i+1,1) = (r+i*h2(j))*cos(j*delta2 + theta);
        Xn((j+np)*(n+1)+i+1,2) = (r+i*h2(j))*sin(j*delta2 + theta);
    end
    j = j+1;
end
Xn;


% filename = 'Input/coordinate.txt' ;
% fid = fopen(filename,'w') ;
% for ino = 1:nno
%     if ino < nno
%         fprintf(fid,'%20.15f \t %20.15f \n',Xn(ino,:));   
%     else
%         fprintf(fid,'%20.15f \t %20.15f',Xn(ino,:));   
%     end
% end
% fclose(fid);


% Current node coordinates initialized as initial coordinates to start
% with.
xn = Xn ; % Xn contains the coordinates of the reference configuration all the time
          % xn contains the coordinates of the current configuration all
          % the time        