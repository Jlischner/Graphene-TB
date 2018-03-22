kxs = linspace(0,1,20);
kys = kxs;

load output;
Nat = size(Cs)(1) + size(Cr)(1);
load Hdat;


ws = [-3:0.1:3]';
dos = zeros(size(ws));
sigma2 = 0.15^2;

ks = [];
eigs = [];

for ix = 1:length(kxs)-1;
  for iy = 1:length(kys)-1;

    kx = kxs(ix);
    ky = kys(iy);

    H = zeros(Nat,Nat);
    
    for nn = 1:size(entries)(1);
      H(entries(nn,1),entries(nn,2)) += entries(nn,3)*exp(2*pi*I*(kx*entries(nn,4) + ky*entries(nn,5) ));
    endfor;

    ens = eig(H);
    for ii = 1:Nat;
      dos += exp(-(ws-ens(ii)).^2/sigma2);
    endfor;

    ks = [ks; kx ky];
    eigs = [eigs real(ens)];
    
  endfor;
endfor;
