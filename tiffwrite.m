function tiffwrite(A,FileStr,idx)
%% Writting a row*column*N matrix to a stack of tiff image
% ****************************************************************
% tiffwrite(A);
% tiffwrite(A,FileStr);
% tiffwrite(A,FileStr,idx);
% ****************************************************************
% Author:ZYD,IBP,CAS 11/30/2009
% Revise:WY,12/2/2014. Converted the 'imwrite' into  'imwritestack', because of the
% imwrite always has no permission to write.

if nargin==1
[filename, pathname] = uiputfile( ...
    {'*.tif;*.tiff', 'All TIF-Files (*.tif,*.tiff)'; ...
        '*.*','All Files (*.*)'}, ...
    'Save Image File');
if isequal([filename,pathname],[0,0])
    return
else
    if ~contains(filename,'.tif')
        filename=strcat(filename,'.tif');
    end
    FileStr = fullfile(pathname,filename);
end
end

if exist(FileStr,'file')
    delete(FileStr);
end

if nargin==3
    index=idx;
    if length(index)==1
        imwritestack(A(:,:,index),FileStr);
    else
        v=index(1):index(2);
        imwritestack(A(:,:,v),FileStr);
    end
else
    imwritestack(A,FileStr);
end
