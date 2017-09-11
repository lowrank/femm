subfolders = {'./class/TriangleMesh/private',...
    './class/QuadMode/private', ...
    './class/FunctionSpace/private', ...
    './class/FormBuilder/private',...
    './class/BC/private',...
    './utility/MeshPartition'};
compiled = 1;

for i = 1:size(subfolders, 2)
    f = subfolders{i};
    subf = dir(f);
    if sum([subf.isdir]) == 0
        disp(f)
        mkdir(f);
        compiled = 0;
    end
end

if compiled == 0
!make all
end

clear;
