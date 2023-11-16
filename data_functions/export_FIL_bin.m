function export_FIL_bin(data,filename)
%__________________________________________________________________________
% Function to export data in .bin format same as FIL OPM data. Useful for
% making a file after pre-processing or for appending data between runs or
% conditions.
% 
% Authors:  Robert Seymour      (rob.seymour@ucl.ac.uk) 
%__________________________________________________________________________

fileID = fopen(filename,'w');

A = data;
A = single(A);
fwrite(fileID,A,'single',0,'b');
fclose(fileID);
end
