%Código Para hacer la transformada de Hermite a las imagenes
close all; clear all; clc;

%Cargando imágenes
ds=imageDatastore('C:\Users\ferch\Desktop\Clasificacion_de_imagenes\Nada\*.jpeg');
fname = ds.Files;
[lx,ly] = size(fname);
for k=1:lx
    iim = im2double(rgb2gray(readimage(ds,k))); %Convierto a escala de grises
    [xx,yy] = size(iim);

%Transformada de Hermite
    D = 2; % Maximum order (pyramid width)
    M = 2; % Number of layers (pyramid height)
    T=1; % Sampling distance (1 <= T <= N)
% Analysis
    H_Im = dht2(iim,M,D,T);
    for i=1:6
        H_Im(2:224,2:224,i) = H_Im(2:224,2:224,i)/max(max(H_Im(2:224,2:224,i)));
    end
    iimm=H_Im(2:224,2:224,5);
    im(2:xx,2:yy,k) = 1-iimm;    
end
[xx1,yy1,zz] = size(im);
iiim = zeros(224,224,lx);

for k=1:lx
    iiim(:,:,k) = im(:,:,k);
end

%Guardando las imágenes
for k=1:lx
    %imwrite(iiim(:,:,k),'Image_%d.jpeg',k,'JPEG');
    imwrite(iiim(:,:,k), sprintf('image_%d.jpg',k));
end
