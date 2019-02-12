%
%  Retrieve quadratic triangulation
%

geom_type = 3;
filename_geo = 'sphere320.q.tri';

fid = fopen(filename_geo,'r');

nverts=0;
nfaces=0;
[nverts] = fscanf(fid,'%d',1);
[nfaces] = fscanf(fid,'%d',1);

verts=zeros(3,nverts);
iqfaces=zeros(6,nfaces);

[verts] = fscanf(fid,'%f',[3,nverts]);
[iqfaces] = fscanf(fid,'%d',[6,nfaces]);

%%%nverts,nfaces

fclose(fid);


%
%  Ellipsoid
%


scale_geo=[25/2, 25/2, 75/2]';
shift_geo=[-40/2, 0, 0]';

verts = [verts .* repmat(scale_geo,1,nverts) - repmat(shift_geo,1,nverts)];


filename_out = 'ellipsoid-25x25x75.q.tri';

fid = fopen(filename_out,'w');

fprintf(fid,'%d %d\n',nverts,nfaces);
fprintf(fid,'%22.15e %22.15e %22.15e\n',verts);
fprintf(fid,'%d %d %d %d %d %d\n',iqfaces);

fclose(fid);


filename_out = 'ellipsoid-25x25x75.a.tri';

fid = fopen(filename_out,'w');

fprintf(fid,'%d %d\n',nverts,nfaces);
fprintf(fid,'%22.15e %22.15e %22.15e\n',verts);
fprintf(fid,'%d %d %d\n',iqfaces(1:3,:));

fclose(fid);


