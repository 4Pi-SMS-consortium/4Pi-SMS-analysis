function [stepnum]=getstepnum(filename)
inded=find(filename=='.',1,'last');
indst=find(filename=='_',1,'last');
stepnum=str2double(filename(indst+1:inded-1));