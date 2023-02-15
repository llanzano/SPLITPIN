
function []=simpleSPLIT_PIN_v1(~);

% User-fiendly code based on the following articles:
% D'Amico, Di Franco, Cerutti et al, A phasor-based approach to improve optical sectioning in any confocal microscope with a tunable pinhole. Microsc Res Tech 2022 
% Di Franco et al, SPLIT-PIN software enabling confocal and super-resolution imaging with a virtually closed pinhole.  Scientific Reports 2023

% Short description:
% The code opens a stack of N frames acquired at tunable (decreasing) pinhole size

% open file (stack of N images, with N>1)
[filename,pathname, filterindex] = uigetfile({'*.tif'},'Please select an image stack');
filenamefull = [pathname, filename(1:end-4)];   
A=simpleSPLIT_readfiletif(filenamefull);
T=size(A,3);
X=size(A,1);
Y=size(A,2);

%initialize parameters
Tg=T; % first frame to sum (default: last frame only)
sm=0.3; % smoothing factor
smp=1;  % smoothing repetitions
smr=0;

% visualize image and estimate threshold
Ngat1=sum(A(:,:,Tg:T),3);
intlist=sort(Ngat1(Ngat1>0));
intmin=intlist(round(length(intlist)*0.5));
Thrg=intmin;

% answer2 = questdlg('Load previously saved parameters?');
% if strcmp(answer2,'Yes')==1

[filenamepatt,pathnamepatt, filterindex] = uigetfile({'*_par.mat'},'Initialize parameters from file');
if filenamepatt ~= 0
filenamefullpatt = [pathnamepatt, filenamepatt];   
load( [filenamefullpatt], 'Tg','sm','smp','smr','Thrg', 'Min','Mout' ) ;
end    

% if N>2 calculate phasor by FFT
if T>2
[gfraw ]=simpleSPLIT_phasor_raw(A(:,:,1:end));   
end

% start processing with interactive menu
imagetoplot='elongation';
figcheck=1;
while figcheck==1
    
%define images 
Ntot=sum(A,3);
% maxc=round(max(max(Ntot)));
% max1=max(max(A(:,:,1)));
Ngat1=sum(A(:,:,Tg:T),3);
Nlast=A(:,:,T);
    
close all
if max(max(Ngat1))<Thrg
    Thrg=0;
end
figure
B=simpleSPLIT_Threshold(Ngat1,Thrg);

%calculate modulation image (mod_ns) from 2-image series or phasor;
%estimate M-in and M-out from modulation histogram
if T==2
[mod_ns, mod, Mmin, maxmod]=simpleSPLIT_calc_mod(A,Tg,Thrg,sm,smp, imagetoplot);
else
[mod_ns , Mmin, maxmod ]=simpleSPLIT_phasor_plot_only(A,Tg, gfraw,B,Thrg,1,sm,smp, imagetoplot)     ;
end
% mod_ns with no threshold, used for SPLIT-PIN
% in the menu, modulation is visualized with out-of-threshold-pixels set to maxmod

prompt2 = {'First frame to sum:','Smoothing factor:','Smoothing cycles phasor/modulation:','Threshold value:'}; 
dlg_title2 = 'Processing parameters (1)'; 
num_lines = 1;
def2 = {num2str(Tg),num2str(sm),num2str(smp),num2str(Thrg)  };
answer = inputdlg(prompt2,dlg_title2,num_lines,def2);
figcheck=~isempty(answer); 
if figcheck==1
Tg=str2num(answer{1});
sm=str2num(answer{2});
smp=str2num(answer{3}); % repetition smooth phasor plot  
Thrg=str2num(answer{4});
end
end


if filenamepatt ~= 0
% modepar='file';
% filenamefullpatt = [pathnamepatt, filenamepatt];   
% load( [filenamefullpatt],  'Min','Mout'  ) ;
else
Min=Mmin;
Mout=maxmod;
end
modepar='manual';

figure
figcheck3=1;
while figcheck3==1 
Ngat1=sum(A(:,:,Tg:T),3);
Ntot=simpleSPLIT_smooth_simple(Ntot,sm,smr);
Ngat1=simpleSPLIT_smooth_simple(Ngat1,sm,smr);
Nlast=simpleSPLIT_smooth_simple(Nlast,sm,smr);
    
%choice of Min and Mout
switch modepar
    case 'file'
    if filenamepatt ~= 0
    filenamefullpatt = [pathnamepatt, filenamepatt];   
    load( [filenamefullpatt],  'Min','Mout'  ) ;
    end
%     case 'calib'
%     PH1size=M(1);
%     PH2size=M(2);
%     [Min, Mout] = simpleSPLIT_2PH_Get_const(PH1size,PH2size) ;
    case 'histo'
    Min=Mmin;
    Mout=maxmod;
    case 'manual'
    Min=Min;
    Mout=Mout;
    
end  

%calculate fractions based on Min and Mout
f=zeros(X,Y,2);
flin=zeros(X,Y,2);
klog=4;
for i=1:X
    for j=1:Y
        flin(i,j,2)=(mod_ns(i,j)-Min)/(Mout-Min);
        flin(i,j,1)=1-flin(i,j,2);
        f(i,j,2)=1/(1+exp(- klog* (flin(i,j,2)-0.5) )) ;
        f(i,j,1)=1-f(i,j,2);
    end
end
for k=1:2
SplitImg(:,:,k)=f(:,:,k).*Ngat1;
SplitImgLin(:,:,k)=flin(:,:,k).*Ngat1;
end

%plot  figures
maxc=round(max(max(A(:,:,1))));
maxcg=round(max(max(Ngat1)));
maxc1=round(max(max(SplitImg(:,:,1))));

subplot(2,3,1)
colormap(hot)
imagesc(Ngat1, [maxcg-maxcg,maxcg])
title('Sum of stack');
axis image

subplot(2,3,2)
colormap(hot)
imagesc(SplitImg(:,:,1), [maxc1-maxc1,maxc1])
title('SPLIT-PIN');
axis image

subplot(2,3,3)
colormap(hot)
imagesc(SplitImg(:,:,2), [maxc1-maxc1,maxc1])
title('out-of-focus');
axis image

%plot  figures with zomm
xzoom=1+floor(X/2);
yzoom=1+floor(Y/2);
xzoom1=1+floor(X/8);
yzoom1=1+floor(Y/8);

subplot(2,3,4)
colormap(hot)
imagesc(Ngat1(xzoom-xzoom1:xzoom+xzoom1,yzoom-yzoom1:yzoom+yzoom1), [maxcg-maxcg,maxcg])
title('Sum of stack');
axis image

subplot(2,3,5)
colormap(hot)
imagesc(SplitImg(xzoom-xzoom1:xzoom+xzoom1,yzoom-yzoom1:yzoom+yzoom1,1), [maxc1-maxc1,maxc1])
title('SPLIT-PIN');
axis image
subplot(2,3,6)
colormap(hot)
imagesc(SplitImg(xzoom-xzoom1:xzoom+xzoom1,yzoom-yzoom1:yzoom+yzoom1,2), [maxc1-maxc1,maxc1])
title('out-of-focus');
axis image



%end here SPLIT-PIN operations   
prompt3 = {'M-in','M-out','First frame to sum','Smoothing cycles image:' };
dlg_title3 = 'SPLIT parameters'; 
num_lines3 = 1;
def3 = {num2str(Min,2),num2str(Mout,2),num2str(Tg),num2str(smr) };
answer3 = inputdlg(prompt3,dlg_title3,num_lines3,def3);
figcheck3=~isempty(answer3); 
if figcheck3==1
Min=str2num(answer3{1});
Mout=str2num(answer3{2});
Tg=str2num(answer3{3});
smr=str2num(answer3{4});
end
end

answer = questdlg('Save data?');
if strcmp(answer,'Yes')==1
filenameout=filenamefull(1:end) ;

    
% export 16-bit tiff image scaled (0-65535 correspond to 0-max) scaling factor indicated in file (for shared version)
A1=double(SplitImg(:,:,1));
MaxValImg=max(A1,[],'all');
% MinValImg=min(A1,[],'all');
MinValImg=0;
Aout=(A1-MinValImg)/(MaxValImg-MinValImg);
Aoutscaled=uint16(Aout*65535);
Factor=round(65535/MaxValImg);
% outputFileName = [filenameout, '_min',num2str(MinValImg,2),'max',num2str(MaxValImg),'.tiff'];
outputFileName = [filenameout, '_SPLITPIN_x',num2str(Factor,3),'.tiff'];
delete outputFileName
% imwrite(Aout, outputFileName);
imwrite(Aoutscaled, outputFileName);

% export 16-bit tiff modulation scaled x1000 (for shared version)
modscaled=uint16( B.*mod_ns * 1000 ) ;
outputFileName2 = [filenameout, '_mod_x1000.tiff'];
delete outputFileName2
imwrite(modscaled, outputFileName2);

% % export images SPLIT in, SPLIT out, selected frames sum
% (requires OMEX library)
% fnameout = [filenameout, '.obf'];
% res = {  SplitImg(:,:,1)' SplitImg(:,:,2)'  Ngat1'  };    % save images to imspector format
% hf = omas_bf_open(fnameout,1);
% omas_bf_write(hf,res);
% omas_bf_close(hf);

% % export images SPLIT in, SPLIT out, selected frames sum, last frame
% % (LINEAR, might have neg values)
% filenameout=filenamefull(1:end) ;
% fnameout = [filenameout, '-LIN.obf'];
% res = {  SplitImgLin(:,:,1)' SplitImgLin(:,:,2)'  Ngat1'  Nlast'};    % save images to imspector format
% hf = omas_bf_open(fnameout,1);
% omas_bf_write(hf,res);
% omas_bf_close(hf);

% export   elongation, mod 
% filenameout=filenamefull(1:end) ;
% fnameout = [filenameout, '_2Fmod.obf'];
% res = { B'.*mod_ns'  mod_ns' };    % save images to imspector format
% hf = omas_bf_open(fnameout,1);
% omas_bf_write(hf,res);
% omas_bf_close(hf);

%save the analysis in Matlab
save([filenameout,'.mat']);

%save main analysis parameters in Matlab
save([filenameout,'_par.mat'], 'Tg','sm','smp','smr','Thrg', 'Min','Mout' );
end


end

function A=simpleSPLIT_readfiletif(filename)

fname = [filename, '.tif'];
info = imfinfo(fname);
nslice = numel(info);

A=imread(fname, 1);  % read tif stack data
for k = 2:nslice
    B = imread(fname, k);
    A=cat(3,A,B);
end


end

function y=simpleSPLIT_smooth_simple(M,sm,n)
y=M;
if sm>0
filt = (1/(8+1/sm))*[1 1 1; 1 1/sm 1; 1 1 1]; % sm factor <=1 
    for i=1:n
    y = filter2(filt,y);
    end
end
    
end
function B=simpleSPLIT_Threshold(A,thr)
  
if length(thr)==1
B=A;
B(B<=thr)=0;
B(B>thr)=1;
else
B=A;
B(B>thr(2))=0;
B(B<=thr(1))=0;
B(B>0)=1;
end

% figure
subplot(2,2,1)
imagesc(A,[min(nonzeros(A)),max(nonzeros(A))])
axis image
colorbar;
title('Intensity')
subplot(2,2,2)
imagesc(B)
colormap(hot)
axis image
colorbar;
title('modulation mask')

% [m,n]=size(A);
% m0=1+floor(m/2);
% n0=1+floor(n/2);
% m1=1+floor(m/8);
% n1=1+floor(n/8);
% subplot(2,2,3)
% imagesc(A(m0-m1:m0+m1,n0-n1:n0+n1),[min(nonzeros(A)),max(nonzeros(A))])
% axis image
% subplot(2,2,4)
% imagesc(B(m0-m1:m0+m1,n0-n1:n0+n1))
% colormap(hot)
% axis image

end

function [mod_ns, mod, Mmin, maxmod]=simpleSPLIT_calc_mod(A,Tg,Thrg,sm,smr, imagetoplot)
% outputs: mod_ns is the modulation; modimg and phase1 are the modulation and phase of the represented phasor. 
% T=1;
A=double(A);
dimimg=size(A);
X=dimimg(1);
Y=dimimg(2);
% if length(dimimg)>2
T=dimimg(3);
% end
DC=sum(A,3);
Ntot=simpleSPLIT_smooth_simple(DC,sm,1);
max1=max(max(A(:,:,1)));
Ngat=sum(A(:,:,Tg:T),3);
mod_ns= (A(:,:,1) -  A(:,:,2))./ (A(:,:,1) + A(:,:,2))  ;    % calc mod from 2 frames
mod_ns = simpleSPLIT_smooth_simple(mod_ns,sm,smr);
maxmod=max(mod_ns(Ngat>Thrg));
mod=mod_ns;
mod(isnan(mod) | Ngat<Thrg)=0;
% %mod image for export
% modimg=mod;
% modimg(isnan(modimg) | Ngat<Thrg)=0;
subplot(2,2,4);
hist(mod_ns(Ngat>Thrg),50)
% calc of histo tails from max value
% [elem, centers]=hist(mod_ns(Ngat>Thrg),50);
% [C,I] = max(elem);
% ind = find(elem>0.15*C, 1, 'first') ;
%     if ind==1
%         Mmin=0;
%     elseif ind==I
%         ind=I-1;
%     else
% %         Mmin=mean([centers(ind),centers(ind-1)]);
%           Mmin=centers(ind);
%     end
% ind2 = find(elem>0.15*C, 1, 'last') ;    
% maxmod=centers(ind2);
% if maxmod<Mmin
%     Mmin=-1;
%     maxmod=1;
% end

% calc of histo tails from percentiles
modlist=sort(mod_ns(Ngat>Thrg & ~isnan(mod_ns)));
Mmin=modlist(round(length(modlist)*0.1));
maxmod=modlist(round(length(modlist)*0.9));

title(['min=',num2str(Mmin,2), ' max=',num2str(maxmod,2)])

subplot(2,2,3);
switch imagetoplot
    case 'modulation'  
    if isfinite(maxmod)
    imagesc(mod,[maxmod-maxmod,maxmod]);
    else
        imagesc(mod);
    end
    title(['modulation'])

    case 'elongation'  
    imagesc(mod,[Mmin,maxmod]);
    title(['modulation'])
    
    case 'phase'  
    imagesc(phase1);
    title(['phase'])
end
axis image
colormap('jet');
colorbar;

end

function [modin, modout] = simpleSPLIT_2PH_Get_const(PH1size,PH2size)

PHval=[ 2 1.5 1 0.5 0.2 ] ;
[PH1, PH2] = meshgrid((PHval));

modavin = [         0    0.1026    0.2659    0.6444    0.9011
   -0.1026         0    0.1676    0.5814    0.8810
   -0.2659   -0.1676         0    0.4600    0.8385
   -0.6444   -0.5814   -0.4600         0    0.6217
   -0.9011   -0.8810   -0.8385   -0.6217         0] ;

modavout = [          0    0.1695    0.3724    0.7824    0.9607
   -0.1695         0    0.2170    0.7085    0.9454
   -0.3724   -0.2170         0    0.5806    0.9179
   -0.7824   -0.7085   -0.5806         0    0.7628
   -0.9610   -0.9458   -0.9188   -0.7628         0 ] ;

[Xq,Yq] = meshgrid( flip(0.2:0.1:2) , flip(0.2:0.1:2) ); 
modincal = interp2(PH1,PH2,modavin',Xq,Yq);
modoutcal = interp2(PH1,PH2,modavout',Xq,Yq);


test=abs(complex((Xq-PH1size),(Yq-PH2size)));
[~, pos1]=min(test(:));
[I_row, I_col] = ind2sub(size(test),pos1)
    
modin=modincal(I_row, I_col);
modout=modoutcal(I_row, I_col);
end


function [mod_gen, Mmin, maxmod ]=simpleSPLIT_phasor_plot_only(A,Tg,gfraw,B,Thrg,h1,sm,smp, imagetoplot)   % INCOMPLETE % mod_ns  (A,Tg,Thrg,sm,smr, imagetoplot)
A=double(A);
dimimg=size(A);
X=dimimg(1);
Y=dimimg(2);
T=dimimg(3);
Ngat=sum(A(:,:,Tg:T),3);
% DC=sum(A,3);
gf=gfraw;
test1=isfinite(gf);
test2=prod(test1,3);
for vk=1:T
    gf(:,:,vk)=gf(:,:,vk).*test2;
end
if T==2
gfsin(:,:,:)=0;
end
  for k=1:T          % smooth g and s
        gf(:,:,k) = simpleSPLIT_smooth_simple(gf(:,:,k),sm,smp);
  end
g_ns=real(gf(:,:,2)) ;
for s=2:T
  g_ns=cat(3, g_ns, real(gf(:,:,ceil(s/2)+1))*0.5*(1+ (-1)^(s+1) ) - imag(gf(:,:,ceil(s/2)+1))*0.5*(1+ (-1)^(s) ));
end
g_ns(isnan(g_ns))=0;
%   
% g_s=real(gf_s(:,:,2)) ;
% for s=2:T
%   g_s=cat(3, g_s, real(gf_s(:,:,ceil(s/2)+1))*0.5*(1+ (-1)^(s+1) ) - imag(gf_s(:,:,ceil(s/2)+1))*0.5*(1+ (-1)^(s) ));
% end
% g_s(isnan(g_s))=0;
harm1=1+2*(h1-1);
mod= (g_ns(:,:,harm1).^2 +  g_ns(:,:,harm1+1).^2 ).^0.5;    % calc mod of stack
maxmod=max(max(mod));
mod(isnan(mod) | B==0)=0;
mod_ns=mod;
phase1=angle(complex( conj(gf(:,:,1+h1)) ));  % calc phase of stack
% phase1=phase1 + pi ; % define from 0 to 2pi
phase1(isnan(phase1) | B==0)=0;
% tphase1=(Tns/(2*pi*h1)) * ( g_ns(:,:,harm1+1)./g_ns(:,:,harm1) );  % calc tau-phase of stack
% tphase1(isnan(tphase1) | B>0)=0;
%mod image for export
modimg=mod;
modimg(isnan(modimg) | B==0)=0;

subplot(2,2,4);
hist(mod_ns(Ngat>Thrg),50)
modlist=sort(mod_ns(Ngat>Thrg));
Mmin=modlist(round(length(modlist)*0.1));
maxmod=modlist(round(length(modlist)*0.9));

% hist(phase1(Ngat>Thrg),50)
% phaselist=sort(phase1(Ngat>Thrg));
% Phasemin=phaselist(round(length(phaselist)*0.1));
% maxphase=phaselist(round(length(phaselist)*0.9));

%V0=[

mod_gen=mod_ns;

title(['min=',num2str(Mmin,2), ' max=',num2str(maxmod,2)])

subplot(2,2,3);
switch imagetoplot
    case 'modulation'  
    if isfinite(maxmod)
    imagesc(mod,[maxmod-maxmod,maxmod]);
    else
        imagesc(mod);
    end
    title(['modulation'])

    case 'elongation'  
    imagesc(mod,[Mmin,maxmod]);
    title(['modulation'])
    
    case 'phase'  
    imagesc(phase1);
    title(['phase'])
end
axis image
colormap('jet');
colorbar;


end

function [gfraw ]=simpleSPLIT_phasor_raw(A)    
A=double(A);
dimimg=size(A);
X=dimimg(1);
Y=dimimg(2);
T=dimimg(3);
DC=sum(A,3);

for k = 2:T
    DC=cat(3,DC,sum(A,3));
end
if T>4
gfraw=fft(A, [], 3)./DC;   %temporal fft of the stack
% gfexc=fft(Iexc,[],2)/sum(Iexc);
else
     for ix=1:X
        for j=1:Y
            for vk=1:T
   gfcos(ix,j,vk)=( sum( squeeze(A(ix,j,:))'.*cos(2*pi*(vk-1)*(0:T-1)/T)) )/DC(ix,j,1);  
   gfsin(ix,j,vk)=( sum( squeeze(A(ix,j,:))'.*sin(2*pi*(vk-1)*(0:T-1)/T)) )/DC(ix,j,1);           
   gfraw(ix,j,vk)=gfcos(ix,j,vk)-1i*gfsin(ix,j,vk);      
            end
        end
     end
end

end