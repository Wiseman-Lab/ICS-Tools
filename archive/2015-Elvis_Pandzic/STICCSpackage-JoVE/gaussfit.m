function [a] = gaussfit(corr,type,pixelsize,whitenoise,radius)

% Usage: a = gaussfit(corr,type,pixelsize,whitenoise);

%set(gcbf,'pointer','watch');

[X,Y] = meshgrid(-((size(corr,2)-1)/2)*pixelsize:pixelsize:((size(corr,2)-1)/2)*pixelsize,-((size(corr,1)-1)/2)*pixelsize:pixelsize:(size(corr,1)-1)/2*pixelsize);
grid = [X Y];

test=(ismember(corr,max(max(corr))));
 for i=1:size(test,3)
[x,y]=find(test(:,:,i));
if length(x)>1 
X0(i,1)=round(mean(x));
Y0(i,1)=round(mean(y));
else
X0(i,1)=x;
Y0(i,1)=y;
end
end


X0 = mod(X0,size(corr,2));

% Find X0 and Y0 are where remainder from mod was zero -- these are set to
%the "max" (ie size) of the corr
X0(ismember(X0,0)) = size(corr,2);
X0(ismember(Y0,0)) = size(corr,1);

% Sets curve fit options, and sets lower bounds for amplitude and beam
% radius to zero
lb = [0 0 -1 min(min(grid)) min(min(grid))];
ub = [];

if nargin < 5
    weights = ones(size(corr));
else
    weights = ones(size(corr));
    for i=1:size(corr,3)
        weights(:,:,i) = circle(size(corr,2),size(corr,1),X0(i), Y0(i),radius);
    end
end
    
% If there's whitenoise, 2 highest values (NB this might be more than
% two points!) in corr func are set to zero, and given no weight in the fit

if strcmp(whitenoise,'y')&&strcmp(type,'2d')
    for j=1:1
            i = find(ismember(corr(:,:,:),max(max(corr(:,:,:)))));
            %ZerChan = i;
            corr(i) = 0;
            weights(i) = 0;
    end
end

y0 = min(min(corr));
y0 = squeeze(y0);
g0 = squeeze(max(max(corr))) - y0;

% wguess = zeros(size(corr,3),1);
% for i=1:size(corr,3)
% [Wy, Wx] = find(ismember(abs((corr(:,:,i)/g0(i) - exp(-1))),min(min(abs(corr(:,:,i)/g0(i) - exp(-1))))));
% Wx = mod(Wx,size(corr,2));
% wguess(i) = mean(( (Wx - X0(i)).^2  + (Wy - Y0(i)).^2   ).^(1/2))*pixelsize;
% end
wguess = 0.4*ones(size(g0));

% Converts from matrix index to LOCATION in pixelsize units
for i=1:size(corr,3)
    X0(i) = X(1,X0(i));
end
for i=1:size(corr,3)
    Y0(i) = Y(Y0(i),1);
end

options = optimset('Display','off');
displayWaitbar = 'n';

if strcmp(displayWaitbar,'y')
    h = waitbar(0,'Fitting correlation functions...');
end

a = zeros(size(corr,3),6);

% Fits each corr func separately
switch lower(type)
    case '2d'
        initguess = [g0 wguess y0 X0 Y0];
        for i=1:size(corr,3)
            if strcmp(displayWaitbar,'y'); waitbar(i/size(corr,3),h); end
            %a0 = initguess(i,:);
            a0xy(1:2) = initguess(i,1:2);
            a0xy(3) = a0xy(2);
            a0xy(4:6) = initguess(i,3:5);
            a(i,:) = lsqcurvefit(@gauss2dwxy,a0xy,grid,corr(:,:,i).*weights(:,:,i),lb,ub,options,weights(:,:,i));
        end
    case 'time'
        initguess = [g0 wguess wguess y0 X0 Y0];
        for i=1:size(corr,3)
            if strcmp(displayWaitbar,'y'); waitbar(i/size(corr,3),h); end
            if i==1
                a0 = initguess(i,:);
            else
                a0 = a(i-1,:);
            end
            % sneak weights into end of grid matrix
            grid(:,(size(X,2)*2+1):(size(X,2)*3)) = weights(:,:,i);
            funlist = {1, @(a,grid) exp(-((grid(:,1:size(grid,2)/3)-a(2)).^2+(grid(:,size(grid,2)/3+1:2*size(grid,2)/3)-a(3)).^2)/(a(1)^2)) .* grid(:,2*size(grid,2)/3+1:end)  };
            NLPstart = [a0(2) a0(5) a0(6)];
            warning('off','MATLAB:rankDeficientMatrix');
            [INLP,ILP] = pleas(funlist,NLPstart,grid,corr(:,:,i),options);
            warning('on','MATLAB:rankDeficientMatrix');
            a(i,1) = ILP(2);
            a(i,2) = INLP(1);
            a(i,3) = INLP(1);
            a(i,4) = ILP(1);
            a(i,5) = INLP(2);
            a(i,6) = INLP(3);
            % "old" way -- all parameters determined in a nonlinear fit
            %[a(i,:),res(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = lsqcurvefit(@gauss2d,a0,grid,corr(:,:,i).*weights(:,:,i),lb,ub,curvefitoptions,weights(:,:,i));
        end
    case 'timeasym'
        initguess = [g0 wguess wguess y0 X0 Y0];
        for i=1:size(corr,3)
            if strcmp(displayWaitbar,'y'); waitbar(i/size(corr,3),h); end
            if i==1
                a0 = initguess(i,:);
            else
                a0 = a(i-1,:);
            end
            funlist = {1, @(a,grid) exp(-(    ((grid(:,1:size(grid,2)/2)-a(3))/a(1)).^2  +  ((grid(:,size(grid,2)/2+1:end)-a(4))/a(2)).^2 ) )}; % instead of 0:  grid(:,1:size(grid,2)/2).*grid(:,size(grid,2)/2+1:end)
            NLPstart = [a0(2) a0(3) a0(5) a0(6)];
            warning('off','MATLAB:rankDeficientMatrix');
           %  options = optimset('disp','iter');

            [INLP,ILP] = pleas(funlist,NLPstart,grid,corr(:,:,i),options);
            warning('on','MATLAB:rankDeficientMatrix');
            a(i,1) = ILP(2);
            a(i,2) = INLP(1);
            a(i,3) = INLP(2);
            a(i,4) = ILP(1);
            a(i,5) = INLP(3);
            a(i,6) = INLP(4);
            %[a(i,:),res(i),RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = lsqcurvefit(@gauss2d,a0,grid,corr(:,:,i).*weights(:,:,i),lb,ub,curvefitoptions,weights(:,:,i));
        end
    otherwise
        error('Fitting mode must be ''2d'', ''time'', or ''timeasym''.');
end

% If the peak moves "past" the edge of the correlation function, it will
% appear on the other side; this unwraps the positions so there are not
% discontinuities in the Gaussian position.  Does it separately for the x-
% and y-coordinates.
if strcmpi(type,'time') || strcmpi(type,'timeasym')
    a(:,5) = unwrapCustom(a(:,5),size(corr,2)/2*pixelsize,1);
    a(:,6) = unwrapCustom(a(:,6),size(corr,1)/2*pixelsize,1);
end

try close(h); catch end % if it waitbar was closed or never open to begin with

%set(gcbf,'pointer','arrow');