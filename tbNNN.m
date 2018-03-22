%# tight-binding for twister bilayer graphene

load output; # load output from blg.m run (t1, t2, Cs, Cr)
C1s = Cs;
C2s = Cr;
N = size(C1s)(1) + size(C2s)(1);

t = -2.7; # nearest neighbor hopping in eV
aNN = 2.46/sqrt(3); # nearest neighbor distance in angstrom
tol = 10^-6;
a1 = 3.35; # distance between planes in angstrom
gamma1 = 0.48; # interlayer hopping energy
qsig = 2.218 * a1;
qpi  = 2.218 * aNN;
rc = 4.14; # cutoff for hopping (6.14 in paper)
lc = 0.265;

A = [t1; t2];
B = 2*pi*inv(A);
b1 = B(:,1)';
b2 = B(:,2)';

%# set up TB Hamiltonian 
%# find nearest neigbohr first
C1b = [];
C2b = [];
indx = [];
partner = [];
for ii = -1:1;
  for jj = -1:1;
    for nn = 1:N/2;
      C1b = [C1b; C1s(nn,:)+ii*t1+jj*t2];
      C2b = [C2b; C2s(nn,:)+ii*t1+jj*t2];
      indx = [indx; ii jj];
      partner = [partner; nn];
    endfor;
  endfor;
endfor;

H = zeros(N,N);
klat = [0.1 0.3];
kx = klat(1);
ky = klat(2);
entries = [];

for ii = 1:N/2;
  
  C1ii = C1s(ii,:);

  # find hoppings from layer 1 to layer 1
  for jj = 1:length(C1b);
    rij = norm( [C1ii 0] - [C1b(jj,:) 0]);
    if( rij < rc + tol && rij > tol);
      Fc = 1/(1+exp((rij-rc)/lc) );
      Vpi = t*exp(qpi * (1-rij/aNN) ) * Fc;
      t12 = Vpi;

      if( norm(indx(jj,:)) > tol);
	H(ii,partner(jj)) += t12*exp(2*pi*I*(kx*indx(jj,1) + ky*indx(jj,2) ));
	entries = [entries; ii partner(jj) t12 indx(jj,1) indx(jj,2)];
      else;
	H(ii,partner(jj)) += t12;
	entries = [entries; ii partner(jj) t12 indx(jj,1) indx(jj,2)];
      endif;
      
    endif;
 
  
  endfor;

  # find hoppings from layer 1 to layer 2  
  for jj = 1:length(C2b);
    rij = norm( [C1ii 0] - [C2b(jj,:) a1]);
    if( rij < rc + tol);
      ndir = a1/rij;
      Fc = 1/(1+exp((rij-rc)/lc) );
      Vpi = t*exp(qpi * (1-rij/aNN) ) * Fc;
      Vsig= gamma1 * exp(qsig * (1-rij/a1) ) * Fc; 
      t12 = ndir^2 * Vsig + (1-ndir^2) * Vpi;

      if( norm(indx(jj,:)) > tol);
	H(ii,partner(jj)+N/2) += t12*exp(2*pi*I*(kx*indx(jj,1) + ky*indx(jj,2) ));
	entries = [entries; ii partner(jj)+N/2 t12 indx(jj,1) indx(jj,2)];
      else;
	H(ii,partner(jj)+N/2) += t12;
	entries = [entries; ii partner(jj)+N/2 t12 indx(jj,1) indx(jj,2)];
      endif;
      
    endif;
    
  endfor;


  # find hoppings from layer 2 to layer 2
  C2ii = C2s(ii,:);
  for jj = 1:length(C2b);
    rij = norm( [C2ii 0] - [C2b(jj,:) 0]);
    if( rij < rc + tol && rij > tol);
      Fc = 1/(1+exp((rij-rc)/lc) );
      Vpi = t*exp(qpi * (1-rij/aNN) ) * Fc;
      t12 = Vpi;

      if( norm(indx(jj,:)) > tol);
	H(ii+N/2,partner(jj)+N/2) += t12*exp(2*pi*I*(kx*indx(jj,1) + ky*indx(jj,2) ));
	entries = [entries; ii+N/2 partner(jj)+N/2 t12 indx(jj,1) indx(jj,2)];
      else;
	H(ii+N/2,partner(jj)+N/2) += t12;
	entries = [entries; ii+N/2 partner(jj)+N/2 t12 indx(jj,1) indx(jj,2)];
      endif;
      
    endif;

  endfor;

  # find hoppings from layer 2 to layer 1
  for jj = 1:length(C1b);
    rij = norm( [C2ii 0] - [C1b(jj,:) a1]);
    if( rij < rc + tol);
      ndir = a1/rij;
      Fc = 1/(1+exp((rij-rc)/lc) );
      Vpi = t*exp(qpi * (1-rij/aNN) ) * Fc;
      Vsig= gamma1 * exp(qsig * (1-rij/a1) ) * Fc; 
      t12 = ndir^2 * Vsig + (1-ndir^2) * Vpi;

      if( norm(indx(jj,:)) > tol);
	H(ii+N/2,partner(jj)) += t12*exp(2*pi*I*(kx*indx(jj,1) + ky*indx(jj,2) ));
	entries = [entries; ii+N/2 partner(jj) t12 indx(jj,1) indx(jj,2)];
      else;
	H(ii+N/2,partner(jj)) += t12;
	entries = [entries; ii+N/2 partner(jj) t12 indx(jj,1) indx(jj,2)];
      endif;
      
    endif;
    
  endfor;
  
endfor;

save Hdat entries;


