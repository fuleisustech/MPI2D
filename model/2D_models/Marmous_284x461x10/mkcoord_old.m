% Author : Lei Fu, March. 2015
% is  is   ig   xs   zs   xg   zg   t 
clear;close all ; clc

addpath('~/utility_function/');
%%
fn_coord_Q = 'coord_csg_full.dat';  % coord 
fn_coord = 'coord_csg.dat';  % coord 
vtrue_str = 'vel.bin';
rho_str = 'den.bin';
vs_str = 'vs.bin';
q_str = 'Q.bin';

vsmooth_fine = 'vp0_fine.bin';
vsmooth_bad = 'vp0_bad.bin';
vsmooth_homo = 'vsmooth_homo.bin';


nz = 284;
nx = 461;
dx = 10; 
dy=dx;

%%

vp = read_bin('vel_nz284_nx461_new.bin',nz,nx);


% figure;imagesc(vp);
% hold on
% plot(Qbd(:,1),Qbd(:,2),'ow')
% % Qbd = ginput();
% save Qbd Qbd
% load 
% [ X, Z ] = meshgrid(1:nx,1:nz);
% [in ] = inpolygon(X,Z,Qbd(:,1),Qbd(:,2));

% Q
% figure
% imagesc(in)
% Q = 10000*ones(size(vp));
% Q(in)=10;
% 
% figure
% imagesc(Q)
% write_bin('Q2.bin',Q); 





Q = read_bin('Q2.bin',nz,nx);
% Q(Q==20)=5;

rho = 0.31*vp.^(0.25);

vs = vp./sqrt(3);
vs(1:37,:) =0;


vp_homo = 1500*ones(size(vp));


write_bin(vtrue_str,vp); 
write_bin(vs_str,vs); 
write_bin(rho_str,rho); 
write_bin(q_str,Q); 
write_bin('vp_homo.bin',vp_homo); 


[vp0_1,~] = vel_smooth(vp,13,9,11);
vp0_1(1:37,:) = 1500;
write_bin(vsmooth_fine,vp0_1); 
er_v = (vp - vp0_1);
er1 = ceil(100*sqrt(sum(er_v(:).^2))/sqrt(sum(vp(:).^2)));


vp00 = linspace(1300,4000,nz)';
vp0_2 = repmat(vp00,1,nx);
vp0_2(1:37,:) = 1500;
write_bin(vsmooth_bad,vp0_2); 
er_v = (vp - vp0_2);
er2 = ceil(100*sqrt(sum(er_v(:).^2))/sqrt(sum(vp(:).^2)));

% 
% vp0_3 = 2500*ones(nz,nx);
% vp0_3(1:37,:) = 1500;
% write_bin(vsmooth_homo,vp0_3); 
% er_v = (vp - vp0_3);
% er3 = ceil(100*sqrt(sum(er_v(:).^2))/sqrt(sum(vp(:).^2)));


% mask=ones(size(vp));
% mask(1:37,:)=0;
% write_bin('mask_bottom.bin',mask); 

% figure
% imagesc(mask);

tapering = zeros(nz,1);
for j = 1:nz
    if j<=50    
       tapering(j) = exp(-(j-50)^2/2/(4)^2) ;
    else 
        tapering(j) =1;
    end
end
figure
plot(tapering);


mask_e=repmat(tapering,1,nx);
mask_e(1:37,:)=0;
figure
imagesc(mask_e);

write_bin('mask_bottom.bin',mask_e); 


%%
figure
subplot(211)
imagesc(vp)
axis image
subplot(212)
imagesc(vp0_2)
axis image

figure
plot(vp(:,200),'k');
hold on
plot(vp0_1(:,200),'b');
plot(vp0_2(:,200),'m');


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
zs0 = 100;        % 1st shot z
nshot = 116;     % shot number 
dsx = 4*dx;       % shot interval
ntrace = 461;   % trace number per shot 
drx = 1*dx;       % receiver interval
xg0 = xs0 + offset; 
zg0 = zs0;  


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
