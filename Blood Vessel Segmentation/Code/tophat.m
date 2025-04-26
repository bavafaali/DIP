function [TH] = tophat(Se_size,I)
    Se = strel('disk',Se_size);
    TH = I - imopen(I,Se);
end