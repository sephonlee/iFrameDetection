function [img] = data2column(im)

[a b c] = size(im);
img = reshape(im,a*b,c);
