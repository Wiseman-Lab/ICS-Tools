
FRETstack16bits = ratio2uint(FRET_Stack, 16);
FRETstack16bits_GF = uint16(zeros(size(FRETstack16bits)));
for i = 1:size(FRETstack16bits, 3)
%     imwrite(FRETstack16bits(:,:,i), 'C:\Users\malvi\OneDrive - McGill University\Data\Hayer\201013\201013-15-ImRatio_16bits.tif', "tif", "WriteMode", "append", "Compression", "none");
    FRETstack16bits_GF(:,:,i) = imgaussfilt(FRETstack16bits(:,:,i), 2);
    imwrite(FRETstack16bits_GF(:,:,i), 'C:\Users\malvi\OneDrive - McGill University\Data\Hayer\201013\201013-15-ImRatio_16bits_GF9.tif', "tif", "WriteMode", "append", "Compression", "none");
end