% Ilwoo Lyu, ilwoolyu@gmail.com
% Release: Apr 26, 2017
% Update: Jul 3, 2017

function [v, f] = read_vtk(filename)
    lines = textread(filename, '%s', 'delimiter', '\n');

    b = find(strncmpi(lines,'POINTS',6)==1) + 1;
    nVert = sscanf(lines{b-1}, 'POINTS %d');
    point_str = lines(b:b+nVert-1, :);
    point_str = sprintf('%s ', point_str{:});
    v=sscanf(point_str, '%f %f %f', [3, nVert])';
    
    if (nargout == 2)
        b = find(strncmpi(lines,'POLYGONS',8)==1) + 1;
        nFaces = sscanf(lines{b-1}, 'POLYGONS %d');
        face_str = lines(b:b+nFaces-1, :);
        face_str = sprintf('%s ', face_str{:});

        f=sscanf(face_str, '%d %d %d %d', [4, nFaces])';
        f=f(:, 2:end);
%         f=uint32(f(:, 2:end));
    end
end