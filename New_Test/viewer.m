clc
clear all
close all

%VIEWER Visualizes the structure - all Chapters
%   
%   Usage:  viewer('plate.mat') or
%           viewer('plate') or
%           viewer plate
%
%   Copyright 2002 AEMM. Revision 2002/03/05 Chapter 2
path = regexprep(mfilename('fullpath'), mfilename(), ''); % path to current script

% filename1 = 'E168_ppQua.txt';
% filename2 = 'E168_ttQua.txt';

filename1 = 'pp_tri_S1.txt';
filename2 = 'tt_tri_S1.txt';

p = dlmread(filename1,'', 0, 0)'; %read data. %amplitude
t = dlmread(filename2,'', 0, 0)'; %phase 

 
mm = length(t);
%Load the data
%load('sphere', '-mat');

%  p = [path, 'pp', '.txt'];
%  t = [path, 'tt', '.txt'];

Color=[0.15 0.55 0.75]; 
%Color='g';
%colormap jet
% Cut=find (t(4,:)>1); %t(4, N) --array of node numbers for each triangle; the number of triangles is N;
% t(:,Cut)=[];  

for m= 1: mm%length(t)
    n=t(1:8, m)
    X(1:8,m)=p(1, n)';  %array of cartesian node coordinates x, y and z; the number of nodes is P;
    Y(1:8,m)=p(2, n)';
    Z(1:8,m)=p(3, n)';
end
%save('E:\Boundary integral equation method and MoM\Spheres with triangle mesh\Output_X.txt', 'X', '-ascii');     
%save ('SphereFu', p, t, '-mat');

g=fill3(X,Y,Z,Color);


axis('equal')
rotate3d on
