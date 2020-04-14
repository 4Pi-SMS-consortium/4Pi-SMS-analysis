function [dmap]=build_dmap(mtc_seg,zang_seg,msz,sigma)

normphi=(zang_seg+pi).*msz./2./pi;
normmtc=(mtc_seg+pi)*msz./2./pi;
[dmap]=binlocalizations([normphi normmtc],msz);
dmap=double(gaussf(dmap,sigma));