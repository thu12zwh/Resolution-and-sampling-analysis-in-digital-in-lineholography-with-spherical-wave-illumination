clc
close all
clear all

% hologram recording
%--------------------------------------------------------------
% illumination wavelength
lam = 500e-6;
k = 2*pi/lam;

% input object
o = imread('test.tif');
o = im2double(o);
o = imresize(o,2);
o = imbinarize(o,0.7);
o = im2double(o);
o = 1-o;

figure,imshow(o,[],'border','tight')

% sampling number 
[m,n] = size(o);

%sampling pitch at the object plane
pitch = 0.001;

% object size
lx = n*pitch;
ly = m*pitch;

% coordination at the object plane in spatial and frequency domain
[x,y] = meshgrid(linspace(-lx/2,lx/2,n),linspace(-ly/2,ly/2,m));
[fx,fy] = meshgrid(linspace(-1/2/pitch,1/2/pitch,n),linspace(-1/2/pitch,1/2/pitch,m));

% source-object distance
z = 3;

% illumiantion spherical wave
r = sqrt(x.^2+y.^2+z^2);
illu = exp(1i*k*r);
illu1 = exp(1i*k/2/z*(x.^2+y.^2));

% complex amplitude behind the object t(x,y)
o_illu = o.*illu;

% zero-padded t(x,y)
o_illu = padarray(o_illu,[2000,2000]);

% sampling number of zero-padded t(x,y)
[m1,n1] = size(o_illu);

% updated parameters for wave propagation
lx1 = n1*pitch;
ly1 = m1*pitch;

[x1,y1] = meshgrid(linspace(-lx1/2,lx1/2,n1),linspace(-ly1/2,ly1/2,m1));
[fx1,fy1] = meshgrid(linspace(-1/2/pitch,1/2/pitch,n1),linspace(-1/2/pitch,1/2/pitch,m1));

% object0-sensor distance
Z = 12;
% transfer function of the angular spectrum method
H = exp(1i*k*Z*sqrt(1-(fx1*lam).^2-(fy1*lam).^2));
% kernel of the convolution-based Fresnel approximation
kernel = exp(1i*k/2/Z*(x1.^2+y1.^2));
kernel_FT = fftshift(fft2(fftshift(kernel)));

o_illu_FT = fftshift(fft2(fftshift(o_illu)));
%complex amplitude of the field at the sensor plane
u = ifftshift(ifft2(ifftshift(o_illu_FT.*kernel_FT)));

%  hologram
h = u.*conj(u);
% downsample the hologram to hologram_1 as in the real situation
pp = 5;
h1 = imresize(h,[m1/pp,n1/pp]);
%--------------------------------------------------------------

% hologram reconstruction with plane wave
[m2,n2] = size(h1);
pitch2 = pitch*pp;
figure,imshow(h1,[],'border','tight')
lx2 = n2*pitch2;
ly2 = m2*pitch2;

[x2,y2] = meshgrid(linspace(-lx2/2,lx2/2,n2),linspace(-ly2/2,ly2/2,m2));
[fx2,fy2] = meshgrid(linspace(-1/2/pitch2,1/2/pitch2,n2),linspace(-1/2/pitch2,1/2/pitch2,m2));

h1_FT = fftshift(fft2(fftshift(h1)));
Z_r1 = 60;
H_r1 = exp(1i*k*Z_r1*sqrt(1-(fx2*lam).^2-(fy2*lam).^2));
kernel_r1 = exp(1i*k/2/Z_r1*(x2.^2+y2.^2));
kernel_r1_FT = fftshift(fft2(fftshift(kernel_r1)));
recons1 = ifftshift(ifft2(ifftshift(h1_FT.*kernel_r1_FT)));
figure,imshow(abs(recons1),[],'border','tight')
a1 = abs(recons1);

%-------------------------------------------------------------------
% hologram reconstruction with original spherical  wave

% oversample hologran_1 to a oversampled hologram_2
qq = 5;
h2 = imresize(h1,[m2*qq, n2*qq]);
pitch3 = pitch2/qq;
[m3,n3] = size(h2);
figure,imshow(h2,[],'border','tight')
lx3 = n3*pitch3;
ly3 = m3*pitch3;

[x3,y3] = meshgrid(linspace(-lx3/2,lx3/2,n3),linspace(-ly3/2,ly3/2,m3));
[fx3,fy3] = meshgrid(linspace(-1/2/pitch3,1/2/pitch3,n3),linspace(-1/2/pitch3,1/2/pitch3,m3));

% model the spherical wave
refer = exp(1i*k*sqrt(x3.^2+y3.^2+(Z+z)^2));


h2_FT = fftshift(fft2(fftshift(h2.*refer)));
Z_r2 = -Z;
H_r2 = exp(1i*k*Z_r2*sqrt(1-(fx3*lam).^2-(fy3*lam).^2));
kernel_r2 = exp(1i*k/2/Z_r2*(x3.^2+y3.^2));
kernel_r2_FT = fftshift(fft2(fftshift(kernel_r2)));
recons2 = ifftshift(ifft2(ifftshift(h2_FT.*H_r2)));
figure,imshow(abs(recons2),[],'border','tight')




