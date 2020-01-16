function f = readmoltxt(filename);

A = importdata(filename,'\t');

%   1       2   3   4   5   6       7       8       9   10  11  12
% Cas44178	X	Y	Xc	Yc	Height	Area	Width	Phi	Ax	BG	I	
%   13  14      15      16      17  18
% Frame	Length	Link	Valid	Z	Zc
mol.cat = A.data(:,1);
mol.x = A.data(:,2);
mol.y = A.data(:,3);
mol.xc = A.data(:,4);
mol.yc = A.data(:,5);
mol.h = A.data(:,6);
mol.area = A.data(:,7);
mol.width = A.data(:,8);
mol.phi = A.data(:,9);
mol.Ax = A.data(:,10);
mol.bg = A.data(:,11);
mol.I = A.data(:,12);
mol.frame = A.data(:,13);
mol.length = A.data(:,14);
mol.link = A.data(:,15);
mol.valid = A.data(:,16);
mol.z = A.data(:,17);
mol.zc = A.data(:,18);

clear A;

f = mol;