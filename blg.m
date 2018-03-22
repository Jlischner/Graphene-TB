a0 = 2.46; # lattice constant in angstrom
a1 = [sqrt(3)/2 -1/2]*a0;
a2 = [sqrt(3)/2 +1/2]*a0;
aNN = a0/sqrt(3);
Nmax = 10;
tol = 10^-10;

%# rotation integers;
nr = 2;
mr = 3;

theta = acos( (nr^2 + 4*nr*mr + mr^2)/2/(nr^2 + nr*mr + mr^2));
printf("for n=%d and m=%d, theta = %f degrees \n",nr,mr,theta*360/(2*pi));

matr = [cos(theta) sin(theta); -sin(theta) cos(theta)]; # note that this was transposed (compared to wikipedia) because it acts on ROW vectors

%# unit cell vectors
t1 = nr*a1 + mr*a2;
t2 = -mr*a1 + (nr+mr)*a2;
T = [t1; t2];

N = 4*(nr^2 + mr^2 + nr*mr);
printf("number of atoms in supercell Nat = %d \n",N);

%# find unit cells with origin in supercell
Cs = [];
Cr = [];
for nc = -Nmax : Nmax;
  for mc = -Nmax : Nmax;

    R1 = nc*a1 + mc*a2;
    tau = R1 * inv(T);

    if( tau(1) < (1-tol) && tau(1) >= -tol && tau(2) < (1-tol) && tau(2) >= -tol);
      Cs = [Cs; R1];
    endif

    R2 = R1 + [aNN 0];
    tau = R2 * inv(T);
    
    if( tau(1) < (1-tol) && tau(1) >= -tol && tau(2) < (1-tol) && tau(2) >= -tol);
      Cs = [Cs; R2];
    endif

    R1r = R1 * matr;
    tau = R1r * inv(T);

    if( tau(1) < (1-tol) && tau(1) >= -tol && tau(2) < (1-tol) && tau(2) >= -tol);
      Cr = [Cr; R1r];
    endif

    R2r = R2 * matr;
    tau = R2r * inv(T);
    
    if( tau(1) < (1-tol) && tau(1) >= -tol && tau(2) < (1-tol) && tau(2) >= -tol);
      Cr = [Cr; R2r];
    endif
    
  endfor;
endfor;

C1s = [];
C2s = [];  
for nx = -Nmax : Nmax ;
  for ny = -Nmax: Nmax ;
    pos1 = nx*a1 + ny*a2;
    pos2 = pos1 + [aNN 0];
    C1s = [C1s; pos1; pos2];
    C2s = [C2s; pos1*matr; pos2*matr];
  endfor;
endfor;

plot(C1s(:,1),C1s(:,2),'bo',Cs(:,1),Cs(:,2),'bo',C2s(:,1),C2s(:,2),'go',Cr(:,1),Cr(:,2),'go',[0 t1(1)],[0 t1(2)],'m-',[0 t2(1)],[0 t2(2)],'m-',[t1(1) t1(1)+t2(1)],[t1(2) t1(2)+t2(2)],'m-',[t2(1) t1(1)+t2(1)],[t2(2) t1(2)+t2(2)],'m-');axis("equal");axis([-10 15 -10 10])
#plot(C1s(:,1),C1s(:,2),'bo',C2s(:,1),C2s(:,2),'ro');axis([30 70 -20 20]);axis("equal");

fid = fopen("positions.dat","w")
fprintf(fid,"%d \n",size(Cs)(1)+size(Cr)(1));

fprintf(fid,"%20.12f %20.12f   0.000000 \n",t1(1),t1(2));
fprintf(fid,"%20.12f %20.12f   0.000000 \n",t2(1),t2(2));
fprintf(fid,"%20.12f %20.12f  20.000000 \n",0.0, 0.0);

for ii = 1:size(Cs)(1);
  fprintf(fid,"C %20.12f %20.12f  1.0000000 \n",Cs(ii,1),Cs(ii,2));
endfor;

for ii = 1:size(Cr)(1);
  fprintf(fid,"C %20.12f %20.12f  4.3500000 \n",Cr(ii,1),Cr(ii,2));
endfor;

fclose(fid);

save output t1 t2 Cs Cr
