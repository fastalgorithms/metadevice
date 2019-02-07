%
%  Retrieve flat triangulation
%

geom_type = 2;
filename_geo = 'sphere320.a.tri';
filename_geo = 'sphere720.a.tri';
filename_geo = '2ellipsoids-25x25x75-sep40.a.tri';
%%%filename_geo = 'sphere11520.a.tri';

fid = fopen(filename_geo,'r');

nverts=0;
nfaces=0;
[nverts] = fscanf(fid,'%d',1);
[nfaces] = fscanf(fid,'%d',1);

verts=zeros(3,nverts);
ifaces=zeros(3,nfaces);

[verts] = fscanf(fid,'%f',[3,nverts]);
[ifaces] = fscanf(fid,'%d',[3,nfaces]);

fclose(fid);

nverts,nfaces

%
%  create triangle vertex and normal arrays
%

ntri = nfaces;
triangles = zeros(3,3,ntri);

for i=1:ntri
        
%triangles(1:3,1,i) = verts(1:3,ifaces(1,i));
%triangles(1:3,2,i) = verts(1:3,ifaces(2,i));
%triangles(1:3,3,i) = verts(1:3,ifaces(3,i));

triangles(1:3,1:3,i) = verts(1:3,ifaces(1:3,i));

end


trisurf(ifaces', verts(1,:)', verts(2,:)', verts(3,:)')
axis equal
