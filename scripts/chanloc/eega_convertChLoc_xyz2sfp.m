% This functions converts the layout xyz to sfp
function eega_convertChLoc_xyz2sfp(fileIn, fileOut)

chanlocs = readlocs(fileIn);
fID = fopen(fileOut,'wt');
for i=1:size(chanlocs,2)
   fprintf(fID,'%s\t%f\t%f\t%f\n',chanlocs(i).labels,chanlocs(i).X,chanlocs(i).Y,chanlocs(i).Z);
end
fclose(fID);