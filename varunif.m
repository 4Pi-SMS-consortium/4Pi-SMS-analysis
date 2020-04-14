function unifim=varunif(im,varmap,sz)
im=single(im);
varmap=single(varmap);
if size(varmap,3)==1
    varmap=repmat(varmap,[1 1 size(im,3)]);
end
tmp1=im./varmap;
if sz==5
    kernel=[1/16,1/4,3/8,1/4,1/16];
    km=kernel'*kernel;
elseif sz==9
    kernel=[1/16,0,1/4,0,3/8,0,1/4,0,1/16];
    km=kernel'*kernel;
else
    km=ones(sz,sz);
end
cim=convn(tmp1,km,'same');
wim=convn(1./varmap(:,:,1),km,'same');
wim=repmat(wim,[1 1 size(cim,3)]);
unifim=cim./wim;