function AdjacencyMatrix = lineSegmentIntersect(xy1,xy2)
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Author: Chiranjib Saha, Harpreet S. Dhillon
% Email: csaha@vt.edu, hdhillon@vt.edu
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% This code was used to make coverage plots in the following paper
% Bibentry goes here ----
%- ------


% Run this script to generate the association probability for a given
% configuration (lambda_bl,L_bl) of blockage distribution. Set these values
% in "parameters.m". 
%   AdjacencyMatrix = lineSegmentIntersect(xy1,xy2) finds the number of
%   intersection points between the set of line segments given in XY1 and XY2.
%   
%   Input:
%   xy1 and xy2 are N1x4 and N2x4 matrices. Rows correspond to line segments. 
%   Each row is of the form [x1 y1 x2 y2] where (x1,y1) is the start point and 
%   (x2,y2) is the end point of a line segment:
%
%                  Line Segment
%       o--------------------------------o
%       ^                                ^
%    (x1,y1)                          (x2,y2)
%
%   Output:
%
%   'AdjacencyMatrix' : N1xN2 indicator matrix where the entry (i,j) is 1 if
%       line segments XY1(i,:) and XY2(j,:) intersect.
%
% The code is a modified version of the code available in ---
% https://www.mathworks.com/matlabcentral/fileexchange/27205-fast-line-segment-intersection
% Author:  U. Murat Erdem
%-------------------------------------------------------------------------------
[n_rows_1,~] = size(xy1);
[n_rows_2,~] = size(xy2);
%%% Prepare matrices for vectorized computation of line intersection points.
%-------------------------------------------------------------------------------
X1 = repmat(xy1(:,1),1,n_rows_2);
X2 = repmat(xy1(:,3),1,n_rows_2);
Y1 = repmat(xy1(:,2),1,n_rows_2);
Y2 = repmat(xy1(:,4),1,n_rows_2);
xy2 = xy2';
X3 = repmat(xy2(1,:),n_rows_1,1);
X4 = repmat(xy2(3,:),n_rows_1,1);
Y3 = repmat(xy2(2,:),n_rows_1,1);
Y4 = repmat(xy2(4,:),n_rows_1,1);
X4_X3 = (X4-X3);
Y1_Y3 = (Y1-Y3);
Y4_Y3 = (Y4-Y3);
X1_X3 = (X1-X3);
X2_X1 = (X2-X1);
Y2_Y1 = (Y2-Y1);
numerator_a = X4_X3 .* Y1_Y3 - Y4_Y3 .* X1_X3;
numerator_b = X2_X1 .* Y1_Y3 - Y2_Y1 .* X1_X3;
denominator = Y4_Y3 .* X2_X1 - X4_X3 .* Y2_Y1;
u_a = numerator_a ./ denominator;
u_b = numerator_b ./ denominator;
% Find the adjacency matrix A of intersecting lines.
INT_X = X1+X2_X1.*u_a;
INT_Y = Y1+Y2_Y1.*u_a;
INT_B = (u_a >= 0) & (u_a <= 1) & (u_b >= 0) & (u_b <= 1);
AdjacencyMatrix = INT_B;
end
