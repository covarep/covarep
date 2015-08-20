% Complete version, October 2004, with minor changes in July 2005 and July
% 2007. Written by Alexis Michaud.
%
% Software for the semi-automatic analysis of the electroglottographic
% signal for a set of continuously voiced items: for example, vowels,
% syllable rhymes, or sustained voiced sounds containing up to <MaxPerN>
% glottal cycles. The value of <MaxPerN> is set at 100 by default.
%
% Takes as input 
% - an EGG signal: either a MONO wav file, or in the right channel of a
% STEREO wav file. Any sampling rate is possible. Recommended: 44,100 Hz or
% higher sampling frequency.
% - a file containing the information concerning the beginning and end of
% each item. Originally designed for reading a .txt file containing a list
% of Regions created with the software SoundForge, as below: 
% Name                                  In           Out      Duration
% -------------------------- ------------- ------------- -------------
% 1                        00:00:03,273  00:00:03,363  00:00:00,090
% 2                        00:00:03,388  00:00:03,490  00:00:00,102
% where "1" and "2" are labels chosen by the user, and the times are given
% as: hours:minutes:seconds,milliseconds.
% 
% Have a look at the comments inside the script <peakdet.m> for finding out
% about parameters that you can modify. For instance, the
% electroglottographic signal is reinterpolated at the closing and opening
% peaks for accurate peak detection; the coefficient for reinterpolation is
% set at 100 by default; this can be modified by the user, to define it as
% a function of the sampling frequency of the original signal, or in
% relation to the fundamental frequency of the sample under analysis.

% % clearing workspace
clear

% setting resampling coefficient
resampC = 100;

% initializing matrix; assumption: there will be no more than 100 periods in each
% analyzed token; this value, which is sufficient for single syllables, 
% can be changed below, in order to treat longer intervals of voicing at one go:
MaxPerN = 100;
% setting coefficient for recognition of "double peaks": for closings, is 0.5,
% following Henrich N., d'Alessandro C., Castellengo M. et Doval B., 2004, "On
% the use of the derivative of electroglottographic signals for characterization 
% of non-pathological voice phonation", Journal of the Acoustical Society of America, 
% 115(3), pp. 1321-1332.
propthresh = 0.5;


% choice of method chosen to handle double closing peaks in Fo calculation:
% if <method == 0>: selecting highest of the peaks
% if <method == 1>: selecting first peak
% if <method == 2>: selecting last peak
% if <method == 3>: using barycentre method
% if <method == 4>: exclude all double peaks
%{
disp(' ')
disp('In case of multiple closing peaks, the value selected in peak detection')
disp('can correspond to the highest peak (enter 0), the first peak (enter 1),  ')
disp('or the last peak (enter 2); or it can use a value in-between (barycentre method; enter 3),')
disp('or cycles with double/multiple closings can be excluded altogether from calculation (enter 4).')
method = input('Your choice: > ');
%}
method = 3;
% choosing maximum possible Fo
%{
disp(' ')
disp('The detection of double peaks requires an F0 ceiling.')
disp('Which value do you propose for this ceiling?')
disp('(i.e. a value slightly above the maximum plausible F0 that could be produced by the speaker)')
disp('Recommended value: 500 Hz.')
maxF = input('Your choice (in Hz): > ');
%}
maxF = 500;
% choosing the smoothing step
%{
disp('Number of points for DEGG smoothing ')
disp('(0: no smoothing, 1: 1 point to the right and left, etc.)')
disp('Recommended value: 1; for noisy signals: up to 3')
smoothingstep = input(' ');
%}
smoothingstep = 3;
% indicating path to EGG file
% disp('Please paste complete path to EGG file here (e.g. D:\EGGsession1\sig.wav)')
% pathEGG = input(' > ','s');
[EGGfilename,pathEGG] = uigetfile('*.*','A:\Coding\Matlab\Data Test\','Please choose the EGG file to be downloaded');

% finding out the characteristics of the sound file
[Y,FS,NBITS] = audioread([pathEGG EGGfilename],1);


%%% loading file that contains beginning and endpoint of relevant
%%% portions of the signal, and the number of the item
% If there exists a text file with the same name as the file of the
% recordings, in the same folder, this file is loaded; otherwise the user
% indicates the path to the text file.
if exist([pathEGG EGGfilename(1:(length(EGGfilename) - 4)) '.txt'])
    textpathfile = pathEGG;
    textfile = [EGGfilename(1:(length(EGGfilename) - 4)) '.txt'];
else
    [textfile,textpathfile] = uigetfile([pathEGG '*.*'],'Please choose the text file that contains the time boundaries of signal portions to be analyzed');
end

% Reading the text file, and retrieving the beginning and end of relevant
% intervals
[LENG,numb] = beginend([textpathfile textfile]);
% LENG = load([textpathfile textfile]);
% retrieving number of lines and columns in LENG:
[NumL,NumC] = size(LENG);

% computing total number of syllables
maxnb = NumL;

% loop for syllables
for i = 1:maxnb
    % loop allowing the user to modify the parameters in view of the results.
    % initialization of an extra variable used to know whether changes
    % must be made in the Oq values. Set at 0 to begin with.
    OqCHAN = 0;

    % assigning default values to <COEF> vector, used in FO: sampling frequency of the
    % electroglottographic recording;
    % smoothing step specified by user; 1, to indicate that the amplitude threshold
    % for peak detection will be set automatically; the fourth value is the threshold value
    % set manually; left at 0 to begin with, as the value will be set automatically
    % and not manually.
    COEF = [FS smoothingstep 1 0];

    % setting the value for threshold for peak detection: by default, half
    % the size of the maximum peak
    propthresh = 0.5; 
        % retrieving time of beginning and end, 
        % and converting from milliseconds to seconds
        time = [LENG(i,1)/1000 LENG(i,2)/1000];
        % if there is a single channel: no difficulty.
        [SIG,FS,NBITS] = audioread([pathEGG EGGfilename],[round(time(1) * COEF(1)) round(time(2) * COEF(1))]);
        results_matrix=peakdet(SIG, FS, maxF, method);
% end of syllable loop
end
 

%%%%%%% saving the results
% clearing unnecessary variables
clear SIG
clear dSIG
clear SdSIG
% disp('Saving the results: ')
% disp('Please type results file name and complete path (e.g. D:\EGGsession1\results1)')
% resname = input(' > ','s'); 