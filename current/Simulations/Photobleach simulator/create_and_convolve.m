
function [matrix_final,particleimage_out,NUM,matrix_conv,matrix] = create_and_convolve(totalsize,den,agg,radius,extra)
    area = (pi*radius^2);  
    matrix = zeros(totalsize+2*extra*radius);
    matrix1 = zeros(totalsize+2*extra*radius,totalsize+2*extra*radius,5);
    particleimage_out = zeros(totalsize+2*extra*radius,totalsize+2*extra*radius);
    for it = 1:length(den)
        matrixtemp    = zeros(size(matrix,1),size(matrix,2));
        num           = round(den(it)*((totalsize+2*extra*radius)^2)/(pi*radius^2));
        particleimage = createImage(totalsize+2*extra*radius,num);       
        NUM(it)       = sum(sum(particleimage((extra*radius+1):(extra*radius+totalsize),(extra*radius+1):(extra*radius+totalsize))));
        matrixtemp    = area*agg(it)*particleimage;
        sizey = size(matrix,1);
        sizex = size(matrix,2);
        for i = 1:sizey
            for j = 1:sizex
                if matrix(i,j) > 0
                else
                    matrix(i,j) = matrixtemp(i,j);
                    particleimage_out(i,j) = particleimage(i,j)*agg(it);
                end
            end
        end
        %matrix  = matrix + matrixtemp;
        %particleimage_out = particleimage_out + particleimage*agg(it);
%         figure(2)
%         colormap(hot(256));
%         imagesc(particleimage_out(:,:))
    end
    matrix_conv  =  round(convolveGaussian(matrix,radius*6,radius));
    matrix_final = matrix_conv((extra*radius+1):(extra*radius+totalsize),(extra*radius+1):(extra*radius+totalsize));
    %matrix = matrix((extra*radius+1):(extra*radius+totalsize),(extra*radius+1):(extra*radius+totalsize));
    %particleimage_out = particleimage_out((extra*radius+1):(extra*radius+totalsize),(extra*radius+1):(extra*radius+totalsize));