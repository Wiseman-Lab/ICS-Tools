function [postConv] = simul8tr2(sizeXdesired,sizeYdesired,sizeT,density,bleachType,bleachDecay,qYield,pixelsize,timesize,PSFType,PSFSize,PSFZ,noBits,diffCoeff,flowX,flowY,flowZ,countingNoise,backgroundNoise,FCS);

% Matlab CLSM simulator
% By DK and SC
% April 23, 2005
% Syntax: [postConv] =
% simul8tr(sizeXdesired,sizeYdesired,sizeT,density,bleachType,bleachDecay,qYield,pixelsize,timesize,PSFType,PSFSize,PSFZ,noBits,diffCoeff,flowX,flowY,flowZ,countingNoise,backgroundNoise);
% * bleachType must be 'none', 'mono', or 'qd'
% for 'mono', bleachDecay should be a single number 0 < k < 1
% for 'qd', bleachDecay should be [mOn mOff], both from ~0.5 to ~1.2
% * PSFType must be 'g' for Gaussian, or 'a' for Airy disk
% * For 2D simulations, set PSFZ = 0

if strcmp(FCS,'y')
    FCS=1;
elseif strcmp(FCS,'n')
    FCS=0;
else
    error('FCS experiment must be either ''y'' or ''n''.')
end
% Adds a "border" around simulation -- cropped at end
sizeX = round(PSFSize/pixelsize)*4+sizeXdesired;
sizeY = round(PSFSize/pixelsize)*4+sizeYdesired;
% sizeX = sizeXdesired;
% sizeY = sizeYdesired;

% Generates an odd (ie, not even!) sized filter for convolution
% '6' is arbitrary factor so Gaussian is small by edge of filter
if (strcmp(PSFType,'g') | strcmp(PSFType,'l'))
    if mod(ceil(PSFSize/pixelsize*6),2)==0
        filtersizeXY = ceil(PSFSize/pixelsize*6)+1;
    else
        filtersizeXY = ceil(PSFSize/pixelsize*6);
    end
elseif strcmp(PSFType,'a')
    % PSFSize is given as the first zero of the airy disk    
    % 10.1735 is the 3rd zero (aka minima) of a Bessel function
    % 3.8317 is the 1st zero
    if mod(ceil(PSFSize/pixelsize/3.8317*7.0156),2)==0
        filtersizeXY = ceil(PSFSize/pixelsize/3.8317*7.0156*2)+1;
    else
        filtersizeXY = ceil(PSFSize/pixelsize/3.8317*7.0156*2);
    end    
else 
    error('PSFType must be either ''g'' or ''a''.')
end

if PSFZ ~= 0
    sizeZ = round(PSFZ/pixelsize*12);
else
    sizeZ = 0;
end 

PSFLUT = fspecial('gaussian',filtersizeXY,PSFSize/pixelsize/2);
    
preConv = zeros(sizeY,sizeX,size(density,2));
% Preallocate output arrays
if ~FCS
postConv = zeros(sizeY,sizeX,sizeT);
elseif FCS
    postConv=zeros(1,sizeT);

end
numToBleach = zeros(size(sizeT,size(density,2)));

% Sets up "population" structure
for i=1:size(density,2)
    % Generate random "seed" for initial positions
    if PSFZ == 0
        numParticles(i) = round(pixelsize^2*sizeX*sizeY*density(i));
    else
        numParticles(i) = round(pixelsize^3*sizeX*sizeY*sizeZ*density(i));
        zCoor = sizeZ*rand(numParticles(i), 1);
        population(i).zCoor = zCoor;
    end    
    xCoor=sizeX*rand(numParticles(i), 1);
    yCoor=sizeY*rand(numParticles(i), 1);
    population(i).xCoor = xCoor;
    population(i).yCoor = yCoor;
    population(i).xCoorDisplay = xCoor;
    population(i).yCoorDisplay = yCoor;
    population(i).diffCoeff = diffCoeff(i);
    population(i).flowX = flowX(i);
    population(i).flowY = flowY(i);
    population(i).flowZ = flowZ(i);
    population(i).qYield = ones(numParticles(i), 1).*qYield(i);
    population(i).blink = ones(numParticles(i), sizeT);
    population(i).numToBleach(1) = 0;
    %Evaluate number of fluorophores left due to photobleaching
    switch lower(bleachType)
        case 'mono'
            bleachHowMany = @(bleachDecay,timesize,numFromBefore,t) (poissrnd(numFromBefore - numFromBefore*exp(-bleachDecay(i)*timesize)));
        case 'bi'
         %   bleachHowMany = @(bleachDecay,timesize,numFromBefore,t) (poissrnd(numFromBefore - numFromBefore*( bleachDecay(1,i)*exp(-bleachDecay(2,i)*timesize) + bleachDecay(3,i)*exp(-bleachDecay(4,i)*timesize)  ))));
        case 'none'
            bleachHowMany = @(bleachDecay,timesize,numFromBefore,t) (0);
        case 'qd'
            mOn = bleachDecay(1);
            mOff = bleachDecay(2);
            numPerStep = 100;
            bleachHowMany = @(bleachDecay,timesize,numFromBefore,t) (0);
            for q=1:numParticles(i)
                population(i).blink(q,:) = QDBlink(mOn,mOff,sizeT,numPerStep);
            end
        otherwise
            error('bleachType must be ''none'', ''mono'', or ''qd''.')
    end
    for j=1:sizeT-1;
        population(i).numToBleach(j+1) = bleachHowMany(bleachDecay,timesize,numParticles(i) - sum(population(i).numToBleach(1:j)),j*timesize);
    end
    
end

set(gcbf,'pointer','watch');
h = waitbar(0,'Simulating...');
%xdiff=zeros(size(density,2),numParticles(i),sizeT);
%ydiff=zeros(size(density,2),numParticles(i),sizeT);
 FCSCounter = 0;
for i = 1:sizeT;
    waitbar(i/(sizeT-1),h);
    for j = 1:size(density,2)
        % Adjust particle positions
        population(j) = simul8trMovement(population(j),timesize,pixelsize,sizeX,sizeY,sizeZ);
%         xdiff(j,:,i)=population(j).xCoor;
%         ydiff(j,:,i)=population(j).yCoor;
%         %Bleaches stuff
        if (population(j).numToBleach(i) ~= 0) & (isnan(population(j).numToBleach(i)) ~= 1)
            population(j).qYield(find(population(j).qYield,population(j).numToBleach(i),'first'))=0;
        end
        % Creates object positions
        if PSFZ ==0 
            preConv(:,:,j)=full(sparse(population(j).yCoorDisplay,population(j).xCoorDisplay, population(j).qYield.*population(j).blink(:,i), sizeY, sizeX));
        elseif PSFZ~=0 
            preConv(:,:,j)=full(sparse(population(j).yCoorDisplay,population(j).xCoorDisplay, population(j).qYield.*exp(-(population(j).zCoor-(sizeZ/2)).^2/(2*(PSFZ/pixelsize/2).^2)).*population(j).blink(:,i), sizeY, sizeX));
          end
    end
    % Convolve new positions
    if strcmp(PSFType,'g') && ~FCS
        postConv(:,:,i) = convolveGaussian(sum(preConv,3),filtersizeXY,PSFSize/pixelsize);
    elseif strcmp(PSFType,'l') && ~FCS
        postConv(:,:,i) = convolveGaussianLine(sum(preConv,3),filtersizeXY,PSFSize/pixelsize);    
    elseif strcmp(PSFType,'a') && ~FCS
        postConv(:,:,i) = convolveAiry(sum(preConv,3),filtersizeXY,PSFSize/pixelsize);
    end
    
    %%% do the FCS
    if FCS
        
         YLine = sizeY/2;
         XPixel = sizeX/2;
         % Convolve new positions
        PSFfilter = zeros(sizeY,sizeX);
        PSFfilter(YLine-(filtersizeXY-1)/2:YLine+(filtersizeXY-1)/2,XPixel-(filtersizeXY-1)/2:XPixel+(filtersizeXY-1)/2) = PSFLUT;
        FCSCounter = FCSCounter + 1;
        %convFCS = imgaussfilt3(sum(preConv,3),[PSFSize PSFSize PSFZ]);
        %postConv(FCSCounter) =sum(convFCS(:));
        postConv(FCSCounter) = sum(sum(PSFfilter.*sum(preConv,3)))*500;
    end   
               
    
    
end

if ~FCS
% Crops the "border" from the simulation
postConv = postConv(round(PSFSize/pixelsize)*2+1:sizeY-round(PSFSize/pixelsize)*2,round(PSFSize/pixelsize)*2+1:sizeX-round(PSFSize/pixelsize)*2,:);
% Adds background noise
if backgroundNoise > 0
    maxIntensity = mean(max(max(postConv)));
    postConv=addBackgroundNoise(postConv,backgroundNoise);
    StoB = maxIntensity/(backgroundNoise*0.60272);
else
    StoB = Inf;
end
postConv = formatSeriesLikeMicroscope(postConv, noBits);
% Adds counting noise
if countingNoise > 0
    postConv=addCountingNoise(postConv,countingNoise);
end

% Adjusts the bit depth of the image
postConv = uint16(formatSeriesLikeMicroscope(postConv, noBits));
end
close(h);