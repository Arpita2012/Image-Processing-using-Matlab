function res = ihp(im,thresh)

% inputs
% im is the fourier transform of the image
% thresh is the cutoff circle radius

%outputs
% res is the filtered image

[r,c]=size(im);
d0=thresh;

d=zeros(r,c);
h=zeros(r,c);

for i=1:r
    for j=1:c
     d(i,j)=  sqrt( (i-(r/2))^2 + (j-(c/2))^2);
    end
end

for i=1:r
    for j=1:c
        if d0>=d(i,j)
          h(i,j)= 0;
        end
         if d0<d(i,j)
          h(i,j)= 1;
        end
    end
end


for i=1:r
    for j=1:c
    res(i,j)=(h(i,j))*im(i,j);

    end
end