function img = filterTimmer(img)


img(255:310,1110:1230) = 0;
        img(200:230,1110:1230) = 0;
        img(250:275,1:125)= 0;
        img(725:830,260:335) = 0;
        img(790:825,1100:1230) = 0;

end