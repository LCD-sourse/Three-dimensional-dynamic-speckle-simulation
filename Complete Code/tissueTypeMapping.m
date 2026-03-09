function [indx,indy,indz]=tissueTypeMapping(x,y,z,X,Y,L)   

indx = round(x+X/2);
if (indx<=0) 
    indx=1;
end
if (indx>=X) 
    indx=X;
end
indy = round(y+Y/2);
if (indy<=0) 
    indy=1;
end
if (indy>=Y) 
    indy=Y;
end
indz = round(z);
if (indz<=0) 
    indz=1;
end
if (indz>=L) 
    indz=L;
end