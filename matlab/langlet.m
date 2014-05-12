function solutionsphere=langlet(x,y,z,r)
% Filename: langlet.m
%
% Version:  Langlet algorithm in 3D v1.0
%
% Purpose:  Function to run Langlet's proposed algorithm for finding
%           the largest sphere that can be inscribed between four others
%           Usage comments are in the code below.
%
% Authors:  Fredrik Boulund
%           Viktor Jonsson
%           Erik Sterner
%
% Date:     2009-06-03
%
% Input:    Vectors x (1x4), y (1x4), z (1x4), with coordinates to the
%           centers of the four spheres along with a vector, r (1x4), 
%           containing the radius of each sphere.
%
% Example:  
%           >> s = langlet([0 1 2 3],[4 3 5 1],[6 4 8 4],[1 1 1 1])
%           
%              s =
% 
%                  2.5000
%                  2.5000
%                  6.5000
%                  1.9580
%

% These are the matrices M, P, S as defined by Langlet (1979) 
% WITH CORRECTION FROM ERROR IN THE ARTICLE! 
M = [ x(2)-x(1) y(2)-y(1) z(2)-z(1) ; ...
      x(3)-x(2) y(3)-y(2) z(3)-z(2) ; ...
      x(4)-x(3) y(4)-y(3) z(4)-z(3) ] ;
Minv = inv(M) ;
P = [ r(1)-r(2) ; ...
      r(2)-r(3) ; ...
      r(3)-r(4) ] ;
S = .5 .* [ (r(1)^2-r(2)^2) + (x(2)^2-x(1)^2) + (y(2)^2-y(1)^2) + (z(2)^2-z(1)^2) ; ...
            (r(2)^2-r(3)^2) + (x(3)^2-x(2)^2) + (y(3)^2-y(2)^2) + (z(3)^2-z(2)^2) ; ...
            (r(3)^2-r(4)^2) + (x(4)^2-x(3)^2) + (y(4)^2-y(3)^2) + (z(4)^2-z(3)^2) ] ;
%% Solution
Vi = [x(1) ; y(1); z(1)] ;          %'Vi' as defined by Langlet
P_ = Minv * P ;                     %'P_' from Langlet
S_ = Minv * S ;                     %'S_' from Langlet
D = S_ - Vi ;                       %'D' from Langlet
% These are the coefficients for the second grade equation,
% as proposed by Langlet
a = P_'*P_ - 1 ;                    
b = P_'* D - r(1) ;
c = D'*D - r(1)^2 ;
% These are the two solutions to the previous mentioned equation. 
% Only r2 is of interest to us.
r2 = -b/a - sqrt(b^2-a*c)/a ;
V = S_ + P_ .* r2 ;
% The solution sphere is stored in variable 'solutionsphere' 
% in the following format:
%        V = [ x-pos ; y-pos ; z-pos ; radius ]
solutionsphere = [V ; r2] ;