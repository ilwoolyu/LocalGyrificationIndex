% Ilwoo Lyu, ilwoolyu@gmail.com
% Release: OCt 25, 2018
% Update: Nov 7, 2018

function OuterHull(vtk_input, vtk_output)
    %% check args
    if nargin ~= 2
        fprintf('Usage: OuterHull(vtk_input_surface, vtk_output_hull)\n');
        return;
    end
    
    %% ensure types
    assert(isa(vtk_input, 'char'))
    assert(isa(vtk_output, 'char'))
    
    %% read vtk and set volume size
    fprintf('read vtk: %s\n', vtk_input);
    [v,f] = read_vtk(vtk_input);
    vox = max(v)-min(v);
    vox = ceil(vox / max(vox) * 256);
    
    %% voxelize mesh using Aitkenhead's tool 
    fprintf('voxelization.. ');
    fv = struct('vertices', v, 'faces', f+1);
    [volume,x,y,z] = VOXELISE(vox(1), vox(2), vox(3), fv, 'xyz');
    volume = padarray(volume, [20, 20, 20]);
    volume = imfill(volume, 26, 'holes');
    fprintf('done\n');

    %% morphological closing
    % convert to double type
    fprintf('morphological closing.. ');
    volume = double(volume);
    volume(volume == 1) = 255;

    % Guassian smoothing to sufficiently envelop the mesh
    if verLessThan('matlab', '8.4')
        volume = gauss3filter(volume,2);
    else
        volume = imgaussfilt3(volume,2);
    end

    % Closing operation
    se = strel('ball', 15, 15);
    volume = imclose(volume, se);

    % Dilation operation
    se = strel('ball', 1, 1);
    volume = imdilate(volume, se);

    % conversion to binary volume
    volume = double(volume > 25) * 255;
    fprintf('done\n');

    %% outer hull
    % Iso-surface generation
    fprintf('outer hull creation.. ');
    [f1, v1] = isosurface(volume, 1);
    
    % Surface adjustment
    v2 = [v1(:,2) v1(:,1) v1(:,3)]; % yxz
    v2 = v2 + repmat(-min(v2), size(v1,1), 1);  % translation to the origin
    v2 = v2 .* repmat((max([max(x) max(y) max(z);max(v)])-min([min(x) min(y) min(z);min(v)]))./max(v2),size(v1,1),1);   % scaling
    v2 = v2 + repmat(min([min(x) min(y) min(z);min(v)]) + [0.5 0 -0.5], size(v1,1),1);  % origin adjustment
    v2 = v2 * 1.02; % ensure full converage
    
    % largest connected component
    A = adjacency(f1);
    [p,~,r] = dmperm(A'+speye(size(A)));
    bins = cumsum(full(sparse(1,r(1:end-1),1,1,size(A,1))));
    bins(p) = bins;

    tab1 = find(bins == 1);
    tab2 = zeros(max(f1(:)),1);
    tab2(tab1) = 1: length(tab1);
    v3 = v2(tab1, :);
    f3 = f1;
    valid = sum(ismember(f1, tab1),2);
    f3(valid < 3, :) = [];
    f3 = tab2(f3);
    fprintf('done\n');

    fprintf('repositioning.. ');
    if verLessThan('matlab', '8.5')
        Options.Verbose = false;
        v4 = ICP_finite(v,v3,Options);
    else
        [~, v4] = pcregrigid(pointCloud(v3), pointCloud(v));
        v4 = v4.Location;
    end
    fprintf('done\n');

    % write output
    fprintf('write vtk: %s\n', vtk_output);
    write_vtk(vtk_output, v4, f3 - 1);
end

function A = adjacency(f)
    n = max(f(:));

    % remove duplicated edges
    rows = [f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)];
    cols = [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)];
    rc = unique([rows,cols], 'rows','first');

    % fill adjacency matrix
    A = sparse(rc(:,1),rc(:,2),1,n,n);    
end
