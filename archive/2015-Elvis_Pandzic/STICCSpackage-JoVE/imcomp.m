function O = imcomp(im1,im2,show,save,titre,num)

im1 = double(im1);      % convert to double for "convenience"
im2 = double(im2);

c3 = zeros(size(im1,1),size(im1,2));
% O = zeros(size(im1,1),size(im1,2),3,size(im1,3));

for  i = 1:size(im1,3);  
    im1 = imadjust(normalisation(im1));
    im2 = imadjust(normalisation(im2));
    
    O(:,:,1,i) = im1;
    O(:,:,2,i) = im2;
    O(:,:,3,i) = c3;
end


if nargin < 3
    show = 'y';
end
if show == 'y'
    scrsz = get(0,'ScreenSize');
    fig = figure('Position',[50 scrsz(4)*.08 scrsz(4)*.85 scrsz(4)*.85]);	% creates square figure
    imagesc(O), axis image
elseif show == 'n'
    scrsz = get(0,'ScreenSize');
    fig = figure('Position',[50 scrsz(4)*.08 scrsz(4)*.85 scrsz(4)*.85],'Visible','Off');
    % creates square figure INVISIBLE
    imagesc(O), axis image
end
if nargin < 4
    save = 'n';
end
if nargin < 5
    titre = [];
end
if nargin < 6
    num = 0;
end

if save == 'y'
    title(titre,'Interpreter','none');
    if num == 0
        a = 1;      % check if a file of the same name exist, if so increment file name
        while exist([pwd,'\composite',num2str(a),'.jpg'],'file') == 2
            a = a+1;
        end
        %     imwrite(O, ['composite',num2str(a),'.jpg'])
        saveas(gcf,['composite',num2str(a),'.jpg'])
    else
        saveas(gcf,['composite',num2str(num),'.jpg'])
    end
end
if show == 'n'
    close(fig)
end

% normalisation function
function matricenorm = normalisation(matrice)
maxi = max(max(matrice));
mini = min(min(matrice));
matricenorm = (matrice-mini)/(maxi-mini);