function imwritestack(stack, filename)
t = Tiff(filename, 'w');
n=size(stack, 3);
tagstruct.ImageLength = size(stack, 1);
tagstruct.ImageWidth = size(stack, 2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
h=waitbar(0,'Writting images,please wait...');
for k = 1:n
    waitbar(k/n);
	t.setTag(tagstruct)
	t.write( uint16(stack(:, :, k)));
	t.writeDirectory();
end
t.close();
close(h);
pause(eps)
