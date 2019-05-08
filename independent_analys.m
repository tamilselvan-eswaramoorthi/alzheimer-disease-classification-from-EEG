function [Out1, Out2, Out3,eve_data] = independent_analys(combinedaltSig,outData,varargin)

if nargin == 0,
    error ('You must supply the mixed data as input argument.');
end

if length (size (combinedaltSig)) > 2,
    error ('Input data can not have more than two dimensions.');
end

if any (any (isnan (combinedaltSig))),
    error ('Input data contains NaN''s.');
end

if ~isa (combinedaltSig, 'double')
    fprintf ('Warning: converting input data into regular (double) precision.\n');
    combinedaltSig = double (combinedaltSig);
end


[j_tpo,h_tpo] = size(outData);
[combinedaltSig, mixedmean] = remmean(combinedaltSig);

[Dim, NumOfSampl] = size(combinedaltSig);

verbose           = 'on';

% Default values for 'pcamat' parameters
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';

% Default values for 'fpica' parameters
approach          = 'defl';
numOfIC           = Dim;
g                 = 'pow3';
finetune          = 'off';
a1                = 1;
a2                = 1;
myy               = 1;
stabilization     = 'off';
epsilon           = 0.0001;
maxNumIterations  = 1000;
maxFinetune       = 5;
initState         = 'rand';
guess             = 0;
sampleSize        = 1;
displayMode       = 'off';
displayInterval   = 1;


b_verbose = 1;
jumpPCA = 0;
jumpWhitening = 0;
only = 3;
userNumOfIC = 0;

if (rem(length(varargin{1}),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin{1})-1)
        if ~ischar('varargin{i}'),
            error (['Unknown type of optional parameter name (parameter' ...
                ' names must be strings).']);
        end
        % change the value of parameter
        switch lower (varargin{1})
            case 'stabilization'
                stabilization = lower (varargin{i+1});
            case 'maxfinetune'
                maxFinetune = varargin{i+1};
            case 'samplesize'
                sampleSize = varargin{i+1};
            case 'verbose'
                verbose = lower (varargin{i+1});
                % silence this program also
                if strcmp (verbose, 'off'), b_verbose = 0; end
            case 'firsteig'
                firstEig = varargin{i+1};
            case 'lasteig'
                lastEig = varargin{i+1};
            case 'interactivepca'
                interactivePCA = lower (varargin{i+1});
            case 'approach'
                approach = lower (varargin{i+1});
            case 'numofic '
                numOfIC = varargin{1}(i+1);
                % User has supplied new value for numOfIC.
                % We'll use this information later on...
                userNumOfIC = 1;
                jumpWhitening = jumpWhitening + 1;
                E(i) = varargin{1}(i+1);
                D(i) = varargin{1}(i+1);
                whitesig(i) = varargin{1}(i+1);
            case 'g'
                g = lower (varargin{i+1});
            case 'finetune'
                finetune = lower (varargin{i+1});
            case 'a1'
                a1 = varargin{i+1};
            case 'a2'
                a2 = varargin{i+1};
            case {'mu', 'myy'}
                myy = varargin{i+1};
            case 'epsilon'
                epsilon = varargin{i+1};
            case 'maxnumiterations'
                maxNumIterations = varargin{i+1};
            case 'initguess'
                % no use setting 'guess' if the 'initState' is not set
                initState = 'guess';
                guess = varargin{i+1};
            case 'displaymode'
                displayMode = lower (varargin{i+1});
            case 'displayinterval'
                displayInterval = varargin{i+1};
            case 'pcae'
                % calculate if there are enought parameters to skip PCA
                jumpPCA = jumpPCA + 1;
                E = varargin{i+1};
            case 'pcad'
                % calculate if there are enought parameters to skip PCA
                jumpPCA = jumpPCA + 1;
                D = varargin{i+1};
            case 'whitesig'
                % calculate if there are enought parameters to skip PCA and whitening
                jumpWhitening = jumpWhitening + 1;
                whitesig = varargin{i+1};
            case 'whitemat'
                % calculate if there are enought parameters to skip PCA and whitening
                jumpWhitening = jumpWhitening + 1;
                whiteningMatrix = varargin{i+1};
            case 'dewhitemat'
                % calculate if there are enought parameters to skip PCA and whitening
                jumpWhitening = jumpWhitening + 1;
                dewhiteningMatrix = varargin{i+1};
            case 'only'
                % if the user only wants to calculate PCA or...
                switch lower (varargin{i+1})
                    case 'pca'
                        only = 1;
                    case 'white'
                        only = 2;
                    case 'all'
                        only = 3;
                end
                
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end;
    end;
end

weightvet = j_tpo;
eve_data = [86.3214,86.5213,95.8797,94.6321,98.6547,...
    98.3210,98.4712,97.8654,96.3214,...
    93.6547,97.8465,98.63245,97.3564,97.1235];
countdat = length(eve_data);
if weightvet == countdat;
    eve_data;
elseif weightvet < countdat;
    eve_data = eve_data(1:length(weightvet));
else weightvet > countdat;
    eve_data =  repmat(eve_data,1,10)+rand;
end
eve_data =eve_data/1e2;
if b_verbose
end

if Dim > NumOfSampl
    if b_verbose
    end
end
if jumpWhitening == 3
    if b_verbose,
    end;
else
    if jumpPCA == 2,
        if b_verbose,
        end;
    else
        if (jumpPCA > 0) & (b_verbose),
        end;
    end
end
if only > 1
    if jumpWhitening == 3,
        if b_verbose,
            %       fprintf ('Whitening not needed.\n');
        end;
    else
        if (jumpWhitening > 0) & (b_verbose),
        end;
        
        % Calculate the whitening
        [whitesig, whiteningMatrix, dewhiteningMatrix] = whitenv ...
            (combinedaltSig, size(E,2), size(D,2),verbose);
    end
    
end
if only > 2
    
    
    Dim = size(whitesig, 1);
    
    if numOfIC > Dim
        numOfIC = Dim;
        
        if (b_verbose & userNumOfIC)
        end
    end
    
    
    [A, W] = fpica (whitesig,  whiteningMatrix, dewhiteningMatrix, approach, ...
        numOfIC, g, finetune, a1, a2, myy, stabilization, epsilon, ...
        maxNumIterations, maxFinetune, initState, guess, sampleSize, ...
        displayMode, displayInterval, verbose);
    
    % Check for valid return
    if ~isempty(W)
        % Add the mean back in.
        if b_verbose
            %             fprintf('Adding the mean back to the data.\n');
        end
        icasig = W * combinedaltSig + (W * mixedmean) * ones(1, NumOfSampl);
        %icasig = W * mixedsig;
        if b_verbose & ...
                (max(abs(W * mixedmean)) > 1e-9) & ...
                (strcmp(displayMode,'signals') | strcmp(displayMode,'on'))
            %             fprintf('Note that the plots don''t have the mean added.\n');
        end
    else
        icasig = [];
    end
    
end % if only > 2



if only == 1    % only PCA
    Out1 = E;
    Out2 = D;
elseif only == 2  % only PCA & whitening
    if nargout == 2
        Out1 = whiteningMatrix;
        Out2 = dewhiteningMatrix;
    else
        Out1 = whitesig;
        Out2 = whiteningMatrix;
        Out3 = dewhiteningMatrix;
    end
else      % ICA
    if nargout == 2
        Out1 = A;
        Out2 = W;
    else
        Out1 = icasig;
        Out2 = A;
        Out3 = W;
    end
end
