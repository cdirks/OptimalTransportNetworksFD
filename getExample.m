function [v_init,range,x,y,bdy] = getExample(n,example)
h = 1/(n-1); 
% Get example 
%-------------------------------------------------------------------------
if ( strcmp(example,'1To2Points') ) % One source and two sinks 
    % Grid and range
    x = 0:h:1; y = 0:h:1;
    range = 0:1/4:1; 
    [X,~,~] = meshgrid(x,y,range);
    % Boundaries 
    bounds = 'quadratic';
    % Initial image 2D
    vInit2D = 0.75*ones(length(y),length(x));
    vInit2D(:,floor(size(vInit2D,2)/2)+1:end) = 0.25;
    vInit2D(end,2:end-1) = 0.5;
%-------------------------------------------------------------------------
elseif ( strcmp(example,'4To4Points') ) % Four sources and four sinks
    % Grid and range
    x = 0:h:1; y = 0:h:1;
    range = 0:1/8:1; 
    [X,~,~] = meshgrid(x,y,range);
    % Boundaries 
    bounds = 'quadratic';
    % Initial image 2D
    vInit2D = 0.25*ones(length(y),length(x));
    sizey = size(vInit2D,2);
    vInit2D(:,1:floor(sizey/5)) = 0.75;
    vInit2D(:,floor(sizey/5)+1:floor(2*sizey/5)) = 0.625;
    vInit2D(:,floor(2*sizey/5)+1:floor(3*sizey/5)) = 0.5;
    vInit2D(:,floor(3*sizey/5)+1:floor(4*sizey/5)) = 0.375;
%-------------------------------------------------------------------------
elseif ( strcmp(example,'16To16Points') ) % 16 sources and 16 sinks
    % Grid and range
    x = 0:h:1; y = 0:h:1;
    range = -1/16:1/16:17/16; 
    [X,~,~] = meshgrid(x,y,range);
    % Boundaries 
    bounds = 'quadratic';
    % Initial image 2D
    vInit2D = zeros(length(y),length(x));
    for i = 1:17
        vInit2D(:,8*(i-1)+1:8*i) = 1-(i-1)/16;
    end
%-------------------------------------------------------------------------
elseif ( strcmp(example,'3To3Circular1') ) % 3 sources and three sinks on circular domain 
    % Grid and range
    x = 0:h:1; y = 0:h:1;
    range = -1/8:1/8:9/8; 
    [X,~,~] = meshgrid(x,y,range);
    % Boundaries 
    bounds = 'circular';
    % Initial image 2D
    vInit2D = zeros(n,n); M = [floor(n/2),floor(n/2)];
    for i = 1:n
        for j = 1:n
            if ( sqrt((i-M(1))^2+(j-M(2))^2) <=floor(n/2)-1 && sqrt((i-M(1))^2+(j-M(2))^2) >= floor(n/2)-5 )
                vInit2D(i,j) = 1;
                xCoord = j-M(2);
                yCoord = M(1)-i;
                if ( xCoord >= 0 )
                    phi = atan(yCoord/xCoord);
                    if ( phi > -pi/2 && phi <= -5*pi/12 )
                        vInit2D(i,j) = 1;
                    elseif ( phi > -5*pi/12 && phi <= -5*pi/12+pi/4+pi/30 )
                        vInit2D(i,j) = 1/4;
                    elseif ( phi > -5*pi/12+pi/4-pi/30 && phi <= 7*pi/12-pi/4+pi/30 )
                        vInit2D(i,j) = 3/8;
                    elseif ( phi > 7*pi/12-pi/4+pi/30 && phi <= pi/2 )
                        vInit2D(i,j) = 1/4;
                    end
                elseif ( xCoord < 0 && yCoord >= 0 )
                    phi = atan(yCoord/xCoord)+pi;
                    if ( phi > pi/2 && phi <= 7*pi/12 )
                        vInit2D(i,j) = 1/4;
                    elseif ( phi > 7*pi/12 && phi <= 7*pi/12+pi/4+pi/30 )
                        vInit2D(i,j) = 1;
                    elseif ( phi > 7*pi/12+pi/4+pi/30 && phi <= pi )
                        vInit2D(i,j) = 7/8;
                    end
                else 
                    phi = atan(yCoord/xCoord)-pi;
                    if ( phi > -pi && phi <= -5*pi/12-pi/4+pi/30 ) 
                        vInit2D(i,j) = 7/8;
                    elseif ( phi > -5*pi/12-pi/4+pi/30 && phi <= -pi/2 )
                        vInit2D(i,j) = 1;
                    end
                end
            end
        end
    end
%-------------------------------------------------------------------------
elseif ( strcmp(example,'3To3Circular2') ) % 3 sources and three sinks on circular domain of different mass 
    % Grid and range
    x = 0:h:1; y = 0:h:1;
    range = -1/8:1/8:9/8; 
    [X,~,~] = meshgrid(x,y,range);
    % Boundaries 
    bounds = 'circular';
    % Initial image 2D
    vInit2D = zeros(n,n); M = [floor(n/2),floor(n/2)];
    for i = 1:100
        for j = 1:100
            if ( sqrt((i-M(1))^2+(j-M(2))^2) <=floor(n/2)-1 && sqrt((i-M(1))^2+(j-M(2))^2) >= floor(n/2)-5 )
                vInit2D(i,j) = 1;
                xCoord = j-M(2);
                yCoord = M(1)-i;
                if ( xCoord >= 0 )
                    phi = atan(yCoord/xCoord);
                    if ( phi >= -pi/2 && phi <= -5*pi/12 )
                        vInit2D(i,j) = 1;
                    elseif ( phi > -5*pi/12 && phi <= -5*pi/12+pi/4+pi/30 )
                        vInit2D(i,j) = 1/2;
                    elseif ( phi > -5*pi/12+pi/4-pi/30 && phi <= 7*pi/12-pi/4+pi/30 )
                        vInit2D(i,j) = 1/4;
                    elseif ( phi > 7*pi/12-pi/4+pi/30 && phi <= pi/2 )
                        vInit2D(i,j) = 1/2;
                    end
                elseif ( xCoord < 0 && yCoord >= 0 )
                    phi = atan(yCoord/xCoord)+pi;
                    if ( phi > pi/2 && phi <= 7*pi/12 )
                        vInit2D(i,j) = 1/2;
                    elseif ( phi > 7*pi/12 && phi <= 7*pi/12+pi/4+pi/30 )
                        vInit2D(i,j) = 1;
                    elseif ( phi > 7*pi/12+pi/4+pi/30 && phi <= pi )
                        vInit2D(i,j) = 7/8;
                    end
                else 
                    phi = atan(yCoord/xCoord)-pi;
                    if ( phi > -pi && phi <= -5*pi/12-pi/4+pi/30 ) 
                        vInit2D(i,j) = 7/8;
                    elseif ( phi > -5*pi/12-pi/4+pi/30 && phi <= -pi/2 )
                        vInit2D(i,j) = 1;
                    end
                end
            end
        end
    end
%-------------------------------------------------------------------------
elseif ( strcmp(example,'1To32Circular') ) % Source in the middle and 32 sinks on circle boundaries 
    % Grid and range
    x = 0:h:1; y = 0:h:1;
    range = -1/32:1/32:33/32; 
    [X,~,~] = meshgrid(x,y,range);
    % Boundaries 
    bounds = 'circularWithMidpoint';
    % Initial image 2D
    vInit2D = zeros(n,n); M = [floor(n/2),floor(n/2)];
    for i = 1:100
        for j = 1:100
            if ( ( sqrt((i-M(1))^2+(j-M(2))^2) <=floor(n/2)-1 && sqrt((i-M(1))^2+(j-M(2))^2) >= floor(n/2)-5 ) || sqrt((i-M(1))^2+(j-M(2))^2) <= 5 )
                vInit2D(i,j) = 1;
                xCoord = j-M(2);
                yCoord = M(1)-i;
                if ( xCoord >= 0 )
                    phi = atan(yCoord/xCoord);
                    if ( phi >= -pi/2 && phi <= -8*pi/18 )
                        vInit2D(i,j) = 8/32;
                    elseif ( phi > -8*pi/18 && phi <= -7*pi/18 )
                        vInit2D(i,j) = 7/32;
                    elseif ( phi > -7*pi/18 && phi <= -6*pi/18 )
                        vInit2D(i,j) = 6/32;
                    elseif ( phi > -6*pi/18 && phi <= -5*pi/18 )
                        vInit2D(i,j) = 5/32;
                    elseif ( phi > -5*pi/18 && phi <= -4*pi/18 )
                        vInit2D(i,j) = 4/32;
                    elseif ( phi > -4*pi/18 && phi <= -3*pi/18 )
                        vInit2D(i,j) = 3/32;  
                    elseif ( phi > -3*pi/18 && phi <= -2*pi/18 )
                        vInit2D(i,j) = 2/32;
                    elseif ( phi > -2*pi/18 && phi <= -1*pi/18 )
                        vInit2D(i,j) = 1/32;
                    elseif ( phi > -1*pi/18 && phi <= 0 )
                        vInit2D(i,j) = 0;
                    elseif ( phi > 0 && phi <= pi/18 )
                        vInit2D(i,j) = 1;
                    elseif ( phi > pi/18 && phi <= 2*pi/18 )
                        vInit2D(i,j) = 31/32;
                    elseif ( phi > 2*pi/18 && phi <= 3*pi/18 )
                        vInit2D(i,j) = 30/32; 
                    elseif ( phi > 3*pi/18 && phi <= 4*pi/18 )
                        vInit2D(i,j) = 29/32;  
                    elseif ( phi > 4*pi/18 && phi <= 5*pi/18 )
                        vInit2D(i,j) = 28/32;  
                    elseif ( phi > 5*pi/18 && phi <= 6*pi/18 )
                        vInit2D(i,j) = 27/32; 
                    elseif ( phi > 6*pi/18 && phi <= 7*pi/18 )
                        vInit2D(i,j) = 26/32; 
                    elseif ( phi > 7*pi/18 && phi <= 8*pi/18 )
                        vInit2D(i,j) = 25/32;  
                    elseif ( phi > 8*pi/18 && phi <= pi/2 )
                        vInit2D(i,j) = 24/32; 
                    end                   
                elseif ( xCoord < 0 && yCoord >= 0 )
                    phi = atan(yCoord/xCoord)+pi;
                    if ( phi > pi/2 && phi <= 10*pi/18 )
                        vInit2D(i,j) = 24/32;
                    elseif ( phi > 10*pi/18 && phi <= 11*pi/18 )
                        vInit2D(i,j) = 23/32;
                    elseif ( phi > 11*pi/18 && phi <= 12*pi/18 )
                        vInit2D(i,j) = 22/32;
                    elseif ( phi > 12*pi/18 && phi <= 13*pi/18 )
                        vInit2D(i,j) = 21/32;
                    elseif ( phi > 13*pi/18 && phi <= 14*pi/18 )
                        vInit2D(i,j) = 20/32;
                    elseif ( phi > 14*pi/18 && phi <= 15*pi/18 )
                        vInit2D(i,j) = 19/32;
                    elseif ( phi > 15*pi/18 && phi <= 16*pi/18 )
                        vInit2D(i,j) = 18/32;
                    elseif ( phi > 16*pi/18 && phi <= 17*pi/18 )
                        vInit2D(i,j) = 17/32;
                    elseif ( phi > 17*pi/18 && phi <= pi )
                        vInit2D(i,j) = 16/32;                    
                    end
                else 
                    phi = atan(yCoord/xCoord)-pi;
                    if ( phi >= -pi && phi < -17*pi/18 )
                        vInit2D(i,j) = 16/32;
                    elseif ( phi >= -17*pi/18 && phi < -16*pi/18 )
                        vInit2D(i,j) = 15/32;
                    elseif ( phi >= -16*pi/18 && phi < -15*pi/18 )
                        vInit2D(i,j) = 14/32;
                    elseif ( phi >= -15*pi/18 && phi < -14*pi/18 )
                        vInit2D(i,j) = 13/32;
                    elseif ( phi >= -14*pi/18 && phi < -13*pi/18 )
                        vInit2D(i,j) = 12/32;
                    elseif ( phi >= -13*pi/18 && phi < -12*pi/18 )
                        vInit2D(i,j) = 11/32;
                    elseif ( phi >= -12*pi/18 && phi < -11*pi/18 )
                        vInit2D(i,j) = 10/32;
                    elseif ( phi >= -11*pi/18 && phi < -10*pi/18 )
                        vInit2D(i,j) = 9/32;
                    elseif ( phi >= -10*pi/18 && phi < -pi/2 )
                        vInit2D(i,j) = 8/32;
                    end
                end
            end
        end
    end
    vInit2D(49,50:end-1) = 1;
    vInit2D(50,50:end-1) = 0;
end
%-------------------------------------------------------------------------
% Boundaries 
if ( strcmp(bounds,'quadratic') )
    bdy = zeros(size(X));
    bdy(1,:,:) = 1; bdy(end,:,:) = 1; bdy(:,1,:) = 1;
    bdy(:,end,:) = 1; bdy(:,:,1) = 1; bdy(:,:,end) = 1;
    bdy = logical(bdy);
elseif ( strcmp(bounds,'circular') )
    bdy = zeros(size(X));
    M = [floor(size(bdy,1)/2),floor(size(bdy,2)/2)];
    for k = 1:size(bdy,3)
        for i = 1:size(bdy,1)
            for j = 1:size(bdy,2)
                if ( sqrt((i-M(1))^2+(j-M(2))^2) >= floor(n/2)-5 )
                    bdy(i,j,k) = 1;
                end
            end
        end
    end
    bdy = logical(bdy);
elseif ( strcmp(bounds,'circularWithMidpoint') )
    bdy = zeros(size(X)); M = [floor(n/2),floor(n/2)];
    for k = 1:size(bdy,3)
        for i = 1:size(bdy,1)
            for j = 1:size(bdy,2)
                if ( sqrt((i-M(1))^2+(j-M(2))^2) >= floor(n/2)-5 || sqrt((i-M(1))^2+(j-M(2))^2) <= 5 )
                    bdy(i,j,k) = 1;
                end
            end
        end
    end
    bdy(floor(n/2)-1:floor(n/2),floor(n/2):end,:) = 1;
    bdy(:,:,1) = 1; bdy(:,:,end) = 1;
    bdy = logical(bdy);
end
% Initial image 3D 
v_init = lift(vInit2D,range);
end