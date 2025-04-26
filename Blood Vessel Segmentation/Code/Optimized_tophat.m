function [TH] = Optimized_tophat(Se1_size,Se2_size,I)
    Icomp = imcomplement(I);
    Se = strel('disk',Se1_size);
    Isc = imopen(Icomp,Se);
    Se2 = strel('disk',Se2_size);
    Ift = imclose(Isc,Se2);
    TH = Icomp - Ift;
end