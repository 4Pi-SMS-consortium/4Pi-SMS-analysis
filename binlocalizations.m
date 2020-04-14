% Binomialnoise Converts position data into a 2d histogram image
%     BINIMAGE = binlocalizations(POSITIONS,XSIZE,YSIZE) returns an image that is made by binning
%     the positions in the Mx2 matrix POSITIONS in pixel bins. The
%     resulting image has dimensions XSIZE x YSIZE. The top left corner of
%     the image corresponds to the coordinate (0,0).

function binimage = binlocalizations(varargin)

if nargin == 1
    positions = varargin{1};
    xsize = 256;
    ysize = 256;
end

if nargin == 2
    positions = round(varargin{1});
    xsize = varargin{2};
    ysize = xsize;
end

if nargin == 3
    positions = round(varargin{1});
    xsize = varargin{2};
    ysize = varargin{3};
end

% Test if emitter positions have been defined as input
if isempty(positions)
    binimage = zeros(xsize,ysize);
    return;
end

% Test if 2D positions have been defined
if (size(size(positions),2) ~= 2)
    binimage = zeros(xsize,ysize);
    return;
end

% Bin emitter locations in pixel bins using a 2D histogram
if (isempty(xsize) && isempty(ysize))
    maxx = ceil(max(positions(:,1)));
    maxy = ceil(max(positions(:,2)));
    binimage = hist3(positions,'CTRS',{1:max(maxx,maxy),1:max(maxx,maxy)})';
else
    %binimage = hist3(positions,'CTRS',{0:(xsize-1),0:(ysize-1)})';
    binimage = cHistRecon(xsize,ysize,single(positions(:,1)),single(positions(:,2)),0);
end