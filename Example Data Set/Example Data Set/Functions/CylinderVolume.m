    function [V, Vmin, Vdiff] = CylinderVolume(Ca,Cyl);
        % This function provides the cylinder volume as function of 
        % Ca : Crankangle [degrees]
        % Cyl :  a struct containing
        % Cyl.S : Stroke
        % Cyl.B                   : Bore
        % Cyl.ConRod              : Connecting Rod length
        % Cyl.CompressionRatio    : Compession Ratio
        % Cyl.TDCangle            : Angle associated with the Top Dead Center
        %
        % Gives output of:
        % V: Cylinder volume in cubic units
        %----------------------------------------------------------------------

        B   = Cyl.Bore;
        S   = Cyl.Stroke;
        cr  = Cyl.CompressionRatio;
        r   = S/2;
        l   = Cyl.ConRod;
        %----------------------------------------------------------------------
 
        % Calculate minimum and maximum cylinder volumes
        Vdiff = pi / 4 * B^2 * S; % Displacement volume
        Vmin = Vdiff / (cr - 1); % Clearance volume
    
        % Calculate instantaneous volume
        x = r + l - (r * cosd(Ca) + sqrt(l^2 - r^2 * (sind(Ca)).^2));
        V = Vmin + x * pi / 4 * B^2;
    end



