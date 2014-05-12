% Filename:	caller.m
%
% Purpose:  Caller file for testing Langlet's algorithm (Langlet, 1979) in the
%           context of finding channels in proteins. This file was produced during a
%           bachelor's thesis project on Chalmers University of Technology in 2009.
%
% Authors:  Fredrik Boulund
%           Viktor Jonsson
%           Erik Sterner
%
% Date:     2009-06-03
%
% Usage:    The program reads a Protein Data Bank (.pdb) file on line 21
%           and after running Langlet 3D2D and Langlet 3D on all 
%           combinations of all atoms contained in Model(1) of the file, 
%           it outputs two BILD files for viewing the results in UCSF
%           Chimera. 
%
% Known limitations: 
%           It is very slow when handling protein files with large number 
%           of atoms. Normally do not run on PDB files containing more than
%           ~200 atoms. Largest protein ever run to date is Model(1) of 
%           1GRM.pdb (260 atoms).

clear all
clc
protein = pdbread('60atoms.pdb');           %Read the protein structure
start = 1;                                  %Start iteration from atom #
stop  = length(protein.Model(1).Atom) ;     %End iteration at atom #
MAXRADIUS = 3.5 ;   %maximum radius of all solutions, in Å
MINRADIUS = 1.15 ;  %minimum radius of all solutions, (HOLE uses 1.15, water molecule)

% Reads Model(1) from the PDB file and stores all atom coordinates in a
% matrix
atomxyz = [[protein.Model(1).Atom.X]',...
           [protein.Model(1).Atom.Y]',...
           [protein.Model(1).Atom.Z]'];
       
%% List of radii associated with atoms of different types
% This cell associates each atom type with a atom radius, read from
% 'atom_types.dat'. 
[A B C] = textread('atom_types.dat', '%s %s %s');   %Read file with atom types
AB = [char(A) char(B)];
C = char(C);
D = ['C3H0';'C3H1';'C4H1';'C4H2';'C4H3';'N3H0';'N3H1';'N3H2';'N4H3';'O1H0';'O2H1';'S2H0';'S2H1'];
D = char(D);
E = [1.61 1.76 1.88 1.88 1.88 1.64 1.64 1.64 1.64 1.42 1.46 1.77 1.77];
radius = zeros(length(protein.Model(1).Atom),1);
%Creates a vector with radius for each atom in protein
for i = 1:length(protein.Model(1).Atom)
    temp1 = strmatch(char([protein.Model(1).Atom(i).resName protein.Model(1).Atom(i).AtomName]),AB) ;
    temp2 = strmatch(C(temp1(1),:),D) ;
    radius(i) = E(temp2) ;
end

%% Langlet "3D2D"
% Calculates the void spheres between combinations of three atoms

tic
counter = 0; %initialization
sphereid2d = zeros((stop*(stop-1)*(stop-2))/(2*3),3);   %prealloc
solutions2d = zeros((stop*(stop-1)*(stop-2))/(2*3),4);  %prealloc
for i=start:stop-2
    for j=i+1:stop-1
        v1 =[(atomxyz(j,1)-atomxyz(i,1))...
             (atomxyz(j,2)-atomxyz(i,2))...
             (atomxyz(j,3)-atomxyz(i,3))];
        for k=j+1:stop
            
            %These operations are to determine the normal to the plane
            %spanned by the three atoms i j and k
            v2 =[(atomxyz(k,1)-atomxyz(i,1))...
                 (atomxyz(k,2)-atomxyz(i,2))...
                 (atomxyz(k,3)-atomxyz(i,3))];
            kryss = [v1(2).*v2(3)-v1(3).*v2(2)
                     v1(3).*v2(1)-v1(1).*v2(3)
                     v1(1).*v2(2)-v1(2).*v2(1)]; %cross(v1,v2);
            enhet = norm(kryss);                 %sqrt(dot(kryss,kryss));
            flytta = (kryss/enhet).*1e-10;

            %create a new imaginary sphere (atom) a very small distance
            %along the direction of the normal
            x = atomxyz(i,1) + flytta(1);
            y = atomxyz(i,2) + flytta(2);
            z = atomxyz(i,3) + flytta(3);

            %Calls the Langlet-function to find the sphere inscribed 
            %by the three (four) spheres/atoms
            temp = langlet([atomxyz(i,1) ...
                            atomxyz(j,1) ...
                            atomxyz(k,1) ...
                            x],...
                            [atomxyz(i,2) ...
                            atomxyz(j,2) ...
                            atomxyz(k,2) ...
                            y],...
                            [atomxyz(i,3) ...
                            atomxyz(j,3) ...
                            atomxyz(k,3) ...
                            z],...
                            [radius(i),radius(j),radius(k),radius(i)] );
            
            %Determine if the atoms are on opposing sides of 
            %the channel - The Same Side Check
            u1 = temp(1:3) - [atomxyz(i,1) atomxyz(i,2) atomxyz(i,3)]';
            u2 = temp(1:3) - [atomxyz(j,1) atomxyz(j,2) atomxyz(j,3)]';
            u3 = temp(1:3) - [atomxyz(k,1) atomxyz(k,2) atomxyz(k,3)]';
            u = u1 + u2 + u3;
            vlength = norm(u); %sqrt(dot(v,v));
           
            if vlength < 1*norm(u1) %sqrt(dot(v1,v1))
                counter = counter + 1;
                
                %store the solution sphere 
                solutions2d(counter,:) = temp;
                
                %Store the three atom serial numbers of the atoms that
                %created this solution sphere
                sphereid2d(counter,:) = [i j k];
            end
        end
    end
end
toc

%% Checks the validity of the "2D" solutions
%removes solutions that overlap with any atoms and put the valid ones into
%the matrix "realsolutions"

tic
fel= 0; %initialization
counter2 = 0; %initialization
validsphereid2d=[]; %initialization
for i = 1:counter
    
    %rough conditions for removing totally incorrect solution spheres
    if isnan(solutions2d(i,1)) == 0 && solutions2d(i,4) < MAXRADIUS && solutions2d(i,4) > MINRADIUS
        
        %Check for steric overlap with all of the atoms in the protein
        for j=start:stop
            v = solutions2d(i,1:3)-[atomxyz(j,1) ...
                                    atomxyz(j,2) ...
                                    atomxyz(j,3)];
            
             if  v*v' + 0.001 < (solutions2d(i,4)+radius(j))^2
                 fel = 1;
                 break
             end
             
        end

        %If the sphere is valid, store it together with its sphereid
        if fel == 0
            counter2 = counter2 + 1;
            realsolutions2d(counter2,:) = solutions2d(i,:);
            validsphereid2d(counter2,:) = sphereid2d(i,:);
        end
        fel = 0;
    end
end
toc

%% Writeout to Chimera BILD-file
%Writes the real solutions to a .bild file viewable in Chimera
tic
fid = fopen('chimeraoutput2d.bild', 'wt') ;
fprintf(fid,'.comment --- This file produced by MATLAB for use as a BILD-file in Chimera\n.color red\n');
for i = 1:counter2
    fprintf(fid,'.sphere %8.3f %8.3f %8.3f %8.3f\n',...
        realsolutions2d(i,1),...
        realsolutions2d(i,2),...
        realsolutions2d(i,3),...
        realsolutions2d(i,4)) ;
end
fclose(fid) ;
toc


%% Langlet 3D
% -------------------------------------------------------------------------
% THIS IS THE SECOND PART OF THE PROGRAM WHICH CALCULATES REAL LANGLET
% SOLUTION SPHERES WITH COMPLETE COMBINATIONS OF FOUR ATOMS - VERY SLOW

%calculating solutions using langlet 3D
tic
counter=0; %initialization
solutions3d = zeros((stop*(stop-1)*(stop-2)*(stop-3))/(2*3*4),4); %prealloc
sphereid3d = zeros((stop*(stop-1)*(stop-2)*(stop-3))/(2*3*4),4);  %prealloc
for i=start:stop-3
    for j=i+1:stop-2
        for k=j+1:stop-1
            for l=k+1:stop
                counter = counter + 1;
                solutions3d(counter,:) = langlet([atomxyz(i,1) ...
                                                  atomxyz(j,1) ...
                                                  atomxyz(k,1) ...
                                                  atomxyz(l,1)],...
                                                  [atomxyz(i,2) ...
                                                  atomxyz(j,2) ...
                                                  atomxyz(k,2) ...
                                                  atomxyz(l,2)],...
                                                  [atomxyz(i,3) ...
                                                  atomxyz(j,3) ...
                                                  atomxyz(k,3) ...
                                                  atomxyz(l,3)],...
                                                  [radius(i),radius(j),radius(k),radius(l)]) ;
                sphereid3d(counter,:) = [i j k l];
            end
        end
    end
end
disp('Calculation of all combinations of four spheres is finished')
toc
%% Checks the validity of the 3D solutions
%removes solutions that overlap with any atoms and put the valid ones into
%the matrix "realsolutions"
tic
fel= 0; %initialization
counter = 0; %initalization
for i = 1:size(solutions3d,1)
    
    %rough conditions for removing totally incorrect solution spheres
    if isnan(solutions3d(i,1)) == 0 && solutions3d(i,4) < MAXRADIUS && solutions3d(i,4) > MINRADIUS
                
        %Check for steric overlap with all of the atoms in the protein
        for j=start:stop
            v = solutions3d(i,1:3)-[atomxyz(j,1) ...
                                    atomxyz(j,2) ...
                                    atomxyz(j,3)];
            
             if  v*v' + 0.001 < (solutions3d(i,4)+radius(j))^2
                 fel = 1;
                 break
             end
             
        end

        %If the sphere is valid, store it together with its
        %sphereid
        if fel == 0
            counter = counter + 1;
            realsolutions3d(counter,:) = solutions3d(i,:);
            validsphereid3d(counter,:) = sphereid3d(i,:);
        end
        fel = 0;
    end
end
toc

%% Writeout to Chimera BILD-file
%Writes the real 3D solutions to a .bild file viewable in Chimera
tic
fid = fopen('chimeraoutput3d.bild', 'wt') ;
fprintf(fid,'.comment --- This file produced by MATLAB for use as a BILD-file in Chimera\n.color blue\n');
for i = 1:counter
    fprintf(fid,'.sphere %8.3f %8.3f %8.3f %8.3f\n',...
        realsolutions3d(i,1),...
        realsolutions3d(i,2),...
        realsolutions3d(i,3),...
        realsolutions3d(i,4)) ;
end
fclose(fid) ;
toc
