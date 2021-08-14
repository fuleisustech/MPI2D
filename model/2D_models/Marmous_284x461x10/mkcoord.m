% Author : Lei Fu, March. 2015
% is  is   ig   xs   zs   xg   zg   t 
clear;close all ; clc

addpath('~/utility_function/');
%%
fn_coord = 'coord_csg_full.dat';  % coord 
fn_coord_Q = 'coord_csg_full_Q.dat';  % coord 

vp_str = 'vp.bin';
rho_str = 'den.bin';
vs_str = 'vs.bin';
vp0_fine_str = 'vp0_fine.bin';
vp0_bad_str = 'vp0_bad.bin';
q_str = 'Q.bin';
mask_str = 'bottom_mask.bin';

dx = 10; 
dy=dx;


dt=1e-3;
nt = 4000;
fc = 15;
sou = ricker(fc,dt,nt);
sou1 = d_integerN(sou,1);

write_bin(['source.bin'],sou1);

%%
vp = read_bin('vel_nz284_nx461_new.bin',284,461);
Q = read_bin('Q_nz284_nx461_new.bin',284,461);


[nz,nx] = size(vp);

rho = 0.31*vp.^(0.25);
vs = vp./sqrt(3);
vs(1:37,:)=0;


vp00 = linspace(1500,4000,nz)';
vp0_bad = repmat(vp00,1,nx);
vp0_bad(1:37,:) = 1500;

[vp0_fine,~] = vel_smooth(vp,13,9,15);
vp0_fine(1:37,:)=1500;

tapering = zeros(nz,1);
for j = 1:nz
    if j<=50    
       tapering(j) = exp(-(j-50)^2/2/(4)^2) ;
    else 
        tapering(j) =1;
    end
end

mask=repmat(tapering,1,nx);
mask = ones(size(vp));
mask(1:37,:)=0;

% er_v = (vp - vp0);
% er = ceil(100*sqrt(sum(er_v(:).^2))/sqrt(sum(vp(:).^2)));
% 
vp_homo = 1500*ones(size(vp));
%%
write_bin(vp_str,vp);
write_bin(vs_str,vs); 
write_bin(rho_str,rho); 
write_bin(q_str,Q); 
write_bin(vp0_fine_str,vp0_fine);   
write_bin(vp0_bad_str,vp0_bad);    

write_bin(mask_str,mask); 
write_bin('vp_homo.bin',vp_homo); 


%%

cw = 0.01;
ch = 0.4;

figure
subplot('position',[0.2 0.3 0.6 0.4]);
imagesc((1:nx)*dx/1000,(1:nz)*dx/1000,vp/1000);
title('Marmousi Model','fontsize',14)
xlabel('X Distance (km)','fontsize',12)
ylabel('Z Distance (km)','fontsize',12)
set(gca,'fontsize',12)
h = colorbar('position',[0.82 0.3 cw ch]);
set(get(h,'title'),'string','km/s');
set(gca,'tickdir','out');

saveas(gca,['model.eps'],'epsc')
print('-dpng','-r300',['model.png']);



%%

cw = 0.01;
ch = 0.3;

figure
subplot('position',[0.2 0.55 0.5 0.3]);
imagesc((0:nx-1)*dx/1000,(0:nz-1)*dx/1000,vp/1000);
title('Marmousi Model','fontsize',14)
% xlabel('X Distance (km)','fontsize',12)
ylabel('Z Distance (km)','fontsize',12)
set(gca,'ytick',[0:1:3])
set(gca,'fontsize',12)
h = colorbar('position',[0.71 0.55 cw ch]);
set(get(h,'title'),'string','km/s');
set(gca,'tickdir','out');
set(h,'ytick',[1:1:4])

subplot('position',[0.2 0.15 0.5 0.3]);
imagesc((0:nx-1)*dx/1000,(0:nz-1)*dx/1000,log10(Q));
title('Q Model','fontsize',14)
xlabel('X Distance (km)','fontsize',12)
ylabel('Z Distance (km)','fontsize',12)
set(gca,'fontsize',12)
set(gca,'ytick',[0:1:3])
h = colorbar('position',[0.71 0.15 cw ch]);
set(get(h,'title'),'string','log_{10}Q');
set(gca,'tickdir','out');
set(h,'ytick',[1:1:4])

saveas(gca,['model.eps'],'epsc')
print('-dpng','-r300',['model.png']);









%% make coord.dat
t = 0.0 ;
offset = 0;     % shot to the first receiver
xs0 = 0;        % 1st shot x
zs0 = 10;        % 1st shot z
nshot = 230;     % shot number 
dsx = 2*dx;       % shot interval
ntrace = 230;   % trace number per shot 
drx = 2*dx;       % receiver interval
xg0 = xs0 + offset; 
zg0 = zs0;  

% for Seis_F90 acquistion
fid = fopen(fn_coord,'wt');
for kk = 1:nshot
    for ii = 1:ntrace         
        
         is = kk;
         ig = ii;
         xs1 = xs0 + (kk-1)*dsx;
         zs = zs0;
         xg = xg0 + (ii-1)*drx;
         zg = zs;
         fprintf(fid,'%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n',is,is,ig,xs1,zs,xg,zg,t);        
        
    end
end
fclose(fid);

% for Dutta acquistion
t = 1.0;
fid = fopen(fn_coord_Q,'wt');
for kk = 1:nshot
    for ii = 1:ntrace
         is = kk;
         ig = ii;
         xs1 = xs0 + (kk-1)*dsx;
         zs = zs0;
         xg = xg0 + (ii-1)*drx;
         zg = zs;
         fprintf(fid,'%d\t%d\t%f\t%f\t%f\t%f\t%f\n',is,ig,xs1,zs,xg,zg,t);        
        
    end
end
fclose(fid);

