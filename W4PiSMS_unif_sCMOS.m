function smoothim=W4PiSMS_unif_sCMOS(im,varim,sz)

if nargin<3
    sz=3;
end
smoothim=varunif(im,varim,sz)-varunif(im,varim,2*sz+3);