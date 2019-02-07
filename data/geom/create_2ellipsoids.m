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
%  Construct two separated ellipsoids
%


sep=40;

scale_geo=[25/2, 25/2, 75/2]';
shift_geo=[-sep/2, 0, 0]';

verts1 = [verts .* repmat(scale_geo,1,nverts) - repmat(shift_geo,1,nverts)];
verts2 = [verts .* repmat(scale_geo,1,nverts) + repmat(shift_geo,1,nverts)];

iqfaces1 = iqfaces;
iqfaces2 = iqfaces + nverts;

nverts = 2*nverts;
nfaces = 2*nfaces;
verts = [verts1 verts2];
iqfaces = [iqfaces1 iqfaces2];


filename_out = ['2ellipsoids-25x25x75-sep' num2str(sep) '.q.tri'];

fid = fopen(filename_out,'w');

fprintf(fid,'%d %d\n',nverts,nfaces);
fprintf(fid,'%22.15e %22.15e %22.15e\n',verts);
fprintf(fid,'%d %d %d %d %d %d\n',iqfaces);

fclose(fid);


filename_out = ['2ellipsoids-25x25x75-sep' num2str(sep) '.a.tri'];

fid = fopen(filename_out,'w');

fprintf(fid,'%d %d\n',nverts,nfaces);
fprintf(fid,'%22.15e %22.15e %22.15e\n',verts);
fprintf(fid,'%d %d %d\n',iqfaces(1:3,:));

fclose(fid);


