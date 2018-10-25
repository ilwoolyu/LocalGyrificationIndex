% Ilwoo Lyu, ilwoolyu@gmail.com
% Release: Apr 26, 2017
% Update: Jul 3, 2017

function write_vtk(fn, v, f)
    fp = fopen(fn,'w');
    fprintf(fp, '# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\n');
    fprintf(fp, 'POINTS %d float\n', size(v, 1));
    fprintf(fp, '%f %f %f\n', v');
    fprintf(fp, 'POLYGONS %d %d\n', size(f, 1), size(f, 1) * 4);
    fprintf(fp, '3 %d %d %d\n', f');
    fclose(fp);
end