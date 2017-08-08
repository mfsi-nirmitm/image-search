

function[gaborSquareEnergy, gaborMeanAmplitude] = CCF(varargin)

       
    [im, nscale, norient, minWaveLength, mult, sigmaOnf,  dThetaOnSigma,k, ...
     polarity] = checkargs(varargin(:));  

    v = version; Octave = v(1) < '5';  % Crude Octave test    
    epsilon         = .0001;         % Used to prevent division by zero.

   
    thetaSigma = pi/norient/dThetaOnSigma;  
    
    [rows,cols] = size(im);
    imagefft = fft2(im);                
    zero = zeros(rows,cols);
    
    totalEnergy = zero;                 
                                        
    totalSumAn  = zero;                
                                        
    orientation = zero;                 
                                        
    estMeanE2n = [];
    EO = cell(nscale, norient);         
    ifftFilterArray = cell(1, nscale);  

    
    
    if mod(cols,2)
	xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
    else
	xrange = [-cols/2:(cols/2-1)]/cols;	
    end
    
    if mod(rows,2)
	yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
    else
	yrange = [-rows/2:(rows/2-1)]/rows;	
    end
    
    [x,y] = meshgrid(xrange, yrange);
                      
    radius = sqrt(x.^2 + y.^2);      
    theta = atan2(-y,x);              
                                      
                                      

    radius = ifftshift(radius);       
    theta  = ifftshift(theta);        
    radius(1,1) = 1;                  
                                      
                                     

    sintheta = sin(theta);
    costheta = cos(theta);
    clear x; clear y; clear theta;    % save a little memory

   
    lp = lowpassfilter([rows,cols], .4, 10);   % Radius .4, 'sharpness' 10

    logGabor = cell(1, nscale);
    
    for s = 1:nscale
        wavelength = minWaveLength*mult^(s-1);
        fo = 1.0/wavelength;                  % Centre frequency of filter.
        logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
        logGabor{s} = logGabor{s}.*lp;        
        logGabor{s}(1, 1) = 0;                 
                                              
    end

    % Then construct the angular filter components...
    spread = cell(1, norient);
    
    for o = 1:norient
        angl = (o-1)*pi/norient;          
        
        
        
        ds = sintheta * cos(angl) - costheta * sin(angl);    
        dc = costheta * cos(angl) + sintheta * sin(angl);    
        dtheta = abs(atan2(ds,dc));                          
        spread{o} = exp((-dtheta.^2) / (2 * thetaSigma^2));  
                                                             
    end

count = 1;
gaborSquareEnergy = [];
gaborMeanAmplitude = [];
    % The main loop...
    for o = 1:norient,                  
        fprintf('Processing orientation %d \r', o); 
        if Octave fflush(1); end
	
        sumAn_ThisOrient  = zero;      
        Energy_ThisOrient = zero;      

        for s = 1:nscale,                 

            filter = logGabor{s} .* spread{o};  
                                                

            ifftFilt = real(ifft2(filter))*sqrt(rows*cols); 
            ifftFilterArray{s} = ifftFilt;                   

            % Convolve image with even and odd filters returning the result in EO
            EO{s, o} = ifft2(imagefft .* filter);
            An  = abs(EO{s,o});                        
            sumAn_ThisOrient = sumAn_ThisOrient + An; 
            


		gaborSquareEnergy(count) = sum(sum( An.^2 ) );
		gaborMeanAmplitude(count) = mean2( An );
 		count = count + 1;

%             end
            
            if s==1
                EM_n = sum(sum(filter.^2)); 
            end                             
        end                                

       

        if polarity == 0    
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                    + abs(real(EO{s,o})) - abs(imag(EO{s,o}));
            end
            
        elseif polarity == 1  % Just look for 'white' spots
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                    + real(EO{s,o}) - abs(imag(EO{s,o}));
            end
            
        elseif polarity == -1  % Just look for 'black' spots
            for s = 1:nscale,                  
                Energy_ThisOrient = Energy_ThisOrient ...
                    - real(EO{s,o}) - abs(imag(EO{s,o}));
            end      
      
        end
        
            
        medianE2n  = median(reshape(abs(EO{1,o}).^2,1,rows*cols));
        meanE2n    = -medianE2n/log(0.5);
        estMeanE2n = [estMeanE2n meanE2n];
        
        noisePower = meanE2n/EM_n; % Estimate of noise power.
        
       
        
        EstSumAn2 = zero;
        for s = 1:nscale
            EstSumAn2 = EstSumAn2+ifftFilterArray{s}.^2;
        end
        
        EstSumAiAj = zero;
        for si = 1:(nscale - 1)
            for sj = (si + 1):nscale
                EstSumAiAj = EstSumAiAj + ifftFilterArray{si} .* ifftFilterArray{sj};
            end
        end
        
        EstNoiseEnergy2 = 2*noisePower*sum(sum(EstSumAn2)) + 4*noisePower*sum(sum(EstSumAiAj));
        
        tau = sqrt(EstNoiseEnergy2/2);                
        EstNoiseEnergy = tau*sqrt(pi/2);              
        EstNoiseEnergySigma = sqrt( (2-pi/2)*tau^2 );
        
        T =  EstNoiseEnergy + k*EstNoiseEnergySigma;  % Noise threshold
        
        
        T = T/1.7;    

        % Apply noise threshold 
        Energy_ThisOrient = max(Energy_ThisOrient - T, zero); 
                          
     
        totalSumAn  = totalSumAn + sumAn_ThisOrient;
        totalEnergy = totalEnergy + Energy_ThisOrient;
      
        
        if(o == 1),
            maxEnergy = Energy_ThisOrient;
        else
            change = Energy_ThisOrient > maxEnergy;
            orientation = (o - 1).*change + orientation.*(~change);
            maxEnergy = max(maxEnergy, Energy_ThisOrient);
        end
        
    end  % For each orientation
    fprintf('                                   \r');
    

    
 
    averageDirectionalEnergy = zero;
   
    
    for sc = 1:nscale
    clear XA;
    clear XE;
    display( sprintf( 'Taking  scale %d/%d', sc, nscale ) );
    scale_current = sc;
    % for a fixed scale, iterate thru each orientation
    for ori=1:norient
        XA(:, :, ori) = abs( EO{scale_current, ori} );
        XE(:, :, ori) = abs( real(EO{scale_current, ori}) ) - abs( imag(EO{scale_current, ori}) );
    end
   
   
     appr_r_XA = reshape( XA, [ size(XA,1)*size(XA,2) norient ] );
     appr_r_median_XA = median( appr_r_XA, 2 );
     mA = reshape( appr_r_median_XA, [size(XA,1) size(XA,2) ] );
     
     appr_r_XE = reshape( XE, [ size(XE,1)*size(XE,2) norient ] );
     appr_r_median_XE = median( appr_r_XE, 2 );
     mE = reshape( appr_r_median_XE, [size(XE,1) size(XE,2) ] );
     

    end
    A = sum( mA, 3 );
    E = sum( mE, 3 );
    
    averageDirectionalEnergy = E ./ (A + epsilon);






    phaseSym = totalEnergy ./ (totalSumAn + epsilon);
    
    % Convert orientation matrix values to degrees
    orientation = orientation * (180 / norient);
    
    

function [im, nscale, norient, minWaveLength, mult, sigmaOnf, ...
          dThetaOnSigma,k, polarity] = checkargs(arg); 

    nargs = length(arg);
    
    if nargs < 1
        error('No image supplied as an argument');
    end    
    
   
    im              = [];
    nscale          = 5;       
    norient         = 6;     
    minWaveLength   = 3;       
    mult            = 2.1;      
    sigmaOnf        = 0.55; 
                             
                             
                                
    dThetaOnSigma   = 1.2;      
                            
                            
                            
    k               = 2.0;  
                            
                             

    polarity        = 0;     

       
    allnumeric   = 1;      
    keywordvalue = 2;       
                            
    readstate = allnumeric;
    
    if readstate == allnumeric
        for n = 1:nargs
            if isa(arg{n}, 'char')
                readstate = keywordvalue;
                break;
            else
                if     n == 1, im            = arg{n}; 
                elseif n == 2, nscale        = arg{n};              
                elseif n == 3, norient       = arg{n};
                elseif n == 4, minWaveLength = arg{n};
                elseif n == 5, mult          = arg{n};
                elseif n == 6, sigmaOnf      = arg{n};
                elseif n == 7, dThetaOnSigma = arg{n};
                elseif n == 8, k             = arg{n};              
                elseif n == 9, polarity      = arg{n};                              
                end
            end
        end
    end

    % Code to handle parameter name - value pairs
    if readstate == keywordvalue
        while n < nargs
            
            if ~isa(arg{n},'char') | ~isa(arg{n+1}, 'double')
                error('There should be a parameter name - value pair');
            end
            
            if     strncmpi(arg{n},'im'      ,2), im =        arg{n+1};
            elseif strncmpi(arg{n},'nscale'  ,2), nscale =    arg{n+1};
            elseif strncmpi(arg{n},'norient' ,2), norient =   arg{n+1};
            elseif strncmpi(arg{n},'minWaveLength',2), minWavelength = arg{n+1};
            elseif strncmpi(arg{n},'mult'    ,2), mult =      arg{n+1};
            elseif strncmpi(arg{n},'sigmaOnf',2), sigmaOnf =  arg{n+1};
            elseif strncmpi(arg{n},'dthetaOnSigma',2), dThetaOnSigma =  arg{n+1};
            elseif strncmpi(arg{n},'k'       ,1), k =         arg{n+1};
            elseif strncmpi(arg{n},'polarity',2), polarity =  arg{n+1};
            else   error('Unrecognised parameter name');
            end

            n = n+2;
            if n == nargs
                error('Unmatched parameter name - value pair');
            end
            
        end
    end
    
    if isempty(im)
        error('No image argument supplied');
    end

    if ~isa(im, 'double')
        im = double(im);
    end
    
    if nscale < 1
        error('nscale must be an integer >= 1');
    end
    
    if norient < 1 
        error('norient must be an integer >= 1');
    end    

    if minWaveLength < 2
        error('It makes little sense to have a wavelength < 2');
    end          
        
    if polarity ~= -1 & polarity ~= 0 & polarity ~= 1
        error('Allowed polarity values are -1, 0 and 1')
    end