%%

%----------------------------------------------------------------------
%                       project sara saljoughi
%----------------------------------------------------------------------

clc
clear
close all
Screen('Preference', 'SkipSyncTests', 1);
[wPtr, rect] = Screen('OpenWindow',max(Screen('Screens')),[128 128 128] );
trialnumber = 3;              %number of trial
N=100;                          %number of color
trial=1; 
respMat = NaN(trialnumber,10);
%----------------------------------------------------------------------
%                       consideration 
%----------------------------------------------------------------------
    
% Start experiment with instructions
myText = 'select the right color by mouse click.\n If you see a circle, you have to remember the top color.\nIf you see a triangle, you have to remember the bottom color.\n try to remember;)\n ESC to exit.\n Press SPACE to start.';
DrawFormattedText(wPtr,myText,'center',rect(4)/2,3);
Screen('Flip',wPtr);
KbWait();                            % wait till space is pressed

% Check the keyboard.
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');

 while trial<trialnumber+1
       
%----------------------------------------------------------------------
%                       fixation point
%----------------------------------------------------------------------
    [x2,  y2] = meshgrid(-128:127,  128:-1:-127);
    M = 127*(1- ((y2 == 0 & x2 > -15 & x2 < 15)|(x2 == 0 & y2 > -15 & y2 < 15 ))) + 1;

    fixation = Screen('MakeTexture',wPtr,M); 
    Screen('DrawTexture',wPtr,fixation);
    Screen('Flip',wPtr);
    WaitSecs ( .5);
%----------------------------------------------------------------------
%                       stimuli
%----------------------------------------------------------------------
   cmp = hsv(N);
    color_panel=255*cmp;
    mrect = [0 0 100 100];
    mrect = CenterRect(mrect,rect);
    color = randperm(N,2);
    color1 = 255*cmp(color(1),:);
    Screen('FillRect',wPtr,color1,mrect-100); 
    color2 = 255*cmp(color(2),:);
    Screen('FillRect',wPtr,color2,mrect+100);
    Screen('Flip',wPtr); 
    WaitSecs ( .5); 
    sample_color = [color1; color2];

%----------------------------------------------------------------------
%                       fixation point
%----------------------------------------------------------------------
    fixation = Screen('MakeTexture',wPtr,M); 
    Screen('DrawTexture',wPtr,fixation);
    Screen('Flip',wPtr);
    delay1 = randperm(2,1);
    if delay1 == 1
       delay1= WaitSecs ( .5);
    else 
        delay1= WaitSecs ( 1);
    end

%----------------------------------------------------------------------
%                       Spatial cue
%----------------------------------------------------------------------
%{
    pic = imread('question.PNG');
    pict =  Screen('MakeTexture',wPtr,pic);
    loc = randperm(2,1);
    if loc == 1
       loc1= mrect-100;
       que_color = sample_color(1, :);
    else 
        loc1 = mrect+100;
        que_color = sample_color(2, :);
    end
    Screen('DrawTexture',wPtr,pict ,[], loc1 );
    Screen('Flip',wPtr);
    WaitSecs ( .3);
    %}
    
    %or
    fixation = Screen('MakeTexture',wPtr,M); 
    Screen('DrawTexture',wPtr,fixation);
    shape = randperm(2,1);
    if shape == 1
       %circle
        Screen('FrameArc', wPtr, [135 10 200],[rect(3)/2-100 rect(4)/2-100 rect(3)/2+100 rect(4)/2+100 ],0,360,10)
       que_color = sample_color(1, :);
    else 
        %triangle
        Screen('DrawLine', wPtr  , [135 10 200],(rect(3)/2)-100  ,(rect(4)/2)-100 ,(rect(3)/2)+100, (rect(4)/2)-100, 10);
       Screen('DrawLine', wPtr  , [135 10 200],(rect(3)/2)-100  ,(rect(4)/2)-100 ,(rect(3)/2), (rect(4)/2)+100, 10);
       Screen('DrawLine', wPtr  , [135 10 200],(rect(3)/2)  ,(rect(4)/2)+100 ,(rect(3)/2)+100, (rect(4)/2)-100, 10);
       
        que_color = sample_color(2, :);
    end
                                       
    Screen('Flip',wPtr);
    WaitSecs ( .3);
    
%----------------------------------------------------------------------
%                       fixation point
%----------------------------------------------------------------------
    fixation = Screen('MakeTexture',wPtr,M); 
    Screen('DrawTexture',wPtr,fixation);
    Screen('Flip',wPtr);
    delay2 = 500+randperm(200,1);
    delay2 = delay2/1000;
    WaitSecs (delay2);


%----------------------------------------------------------------------
%                       target
%----------------------------------------------------------------------

%creating ring
mrect = [rect(3)/2-200 rect(4)/2-200 rect(3)/2+200 rect(4)/2+200 ];

for i=1:N
    color1 = 255*cmp(i,:);
    
    degree = 360/N;           %degree of each sector

    Screen('FrameArc', wPtr, color1,mrect,i*degree,degree ,80)
end
 Screen('Flip',wPtr);
    
%----------------------------------------------------------------------
%                      %subject response 
%---------------------------------------------------------------------- 
    t1 = GetSecs();   
%mouse position
    ShowCursor()
    Xcenter    = rect(3)/2;
    Ycenter    = rect(4)/2;
    buttons = 0;
    while ~buttons
        [keyIsDown,secs, keyCode] = KbCheck(-1);
        if keyCode(escapeKey)
            ShowCursor()
            sca;
            return 
         end
         if GetSecs()-t1>20                   %time limit for response is 20 secs
            ShowCursor()
            sca;
            return 
         end
         [xM,yM,buttons] = GetMouse();
         Object = complex(Xcenter-xM,Ycenter-yM);
         R_Object = abs(Object);
         theta_Object = angle(Object)*(180/pi);
         t2 = GetSecs();
         if theta_Object < 0
             theta_Object = theta_Object + 360;
         end
         theta_Object=theta_Object-90;
         
         
    end    
     
        if theta_Object < 0
             theta_Object = theta_Object + 360;
        end
        if R_Object>216 || R_Object<110
            ShowCursor;
            sca;
            return 
         end
           
    RT = t2-t1;
    
%----------------------------------------------------------------------
%                      %analysis level 1
%---------------------------------------------------------------------- 
 
    
flag = 0;
dif_angle = 0;
Q = theta_Object/(360/N);   
Q = floor(Q);               %selected by user
if Q==0
 
    selectedcolor = color_panel(end,:);
else
    selectedcolor = color_panel(Q,:);
end
W=0;
for j=1:size(color_panel,1)
   if color_panel(j,:)== que_color
       W = j;           %selected by computer
   end
end

if selectedcolor == que_color
    flag = 1;
else
    dif_angle = (Q-W)*(360/N);
    if dif_angle>180
       dif_angle = -360+ dif_angle;
        
    end
    if dif_angle<-180
       dif_angle = 360+ dif_angle;
        
    end
end


%----------------------------------------------------------------------
%                      % Record the response
%----------------------------------------------------------------------

respMat(trial, 1) = flag;          %correct response or not
respMat(trial, 2) = dif_angle ;    %Angle difference
respMat(trial, 3) = shape ;        %shape=1 is for circle, shape=2 is for triangle
respMat(trial, 4) = RT ;           %reaction time
respMat(trial, 5) = Q ;
respMat(trial, 6) = W ;
respMat(1, 7) = N ;
respMat(trial, 8) = Q - W; 

trial = trial+1;

 end
 
    x=[];
    y=[];
    Xcenter    = rect(3)/2+10;
    Ycenter    = rect(4)/2-10;
    R1 = 0.5 ;                  % inner radius 
    R2 = 1 ;                    % outer radius
    nT = linspace(0,2*pi,N) ;
    nR = linspace(R1,R2,2) ;
    [R, T] = meshgrid(nR,nT) ;
    R = 250;
    alfa = 2*pi / N;
    X = R.*cos(T); 
    Y = R.*sin(T);
 for j=1:2
        for i=1:N
            x(i) = Xcenter+ R*sin (nT(i));
            y(i) = Ycenter+ R*cos (nT(i));
            color4 = 255*cmp(i,:);
            x1(i) = Xcenter+ R*sin(nT(i)+pi/2+alfa);
            y1(i) = Ycenter+ R*cos (nT(i)+pi/2+alfa);
            Screen('FillOval', wPtr, color4,[x1(i)-30 y1(i)-30 x1(i)+30 y1(i)+30 ]);
            Screen('TextFont',wPtr,'Helvetica');
            Screen('TextSize',wPtr,20);
            Screen('DrawText',wPtr,'Experiment Finished,thanks:)',rect(3)/2-150,rect(4)/2);
             Screen('Flip',wPtr);
        end
 end
 
[soundData, freq] = audioread('Applause.wav');
%Prepare sound data (make it two rows for stereo playback)
soundData = [soundData, soundData];
numChannels=2;
%Open the audio driver
pahandle = PsychPortAudio('Open',[], [], 0, freq,numChannels);
%Fill the buffer
PsychPortAudio('FillBuffer',pahandle,soundData');
%Play the sound
playTime = PsychPortAudio('Start',pahandle);

clear Screen
save VisualWorkingMemory_test.mat respMat; 
%%

%----------------------------------------------------------------------
%                      % load data and analysis level 2
%---------------------------------------------------------------------- 
clear
close all
load VisualWorkingMemory_samira.mat
angular_error = respMat(:,2)';
trials =size(respMat,1);
N= respMat(1,7);
degree = 360/N;
kk=1;
x=-180:degree:(180-degree);
degree=360/N;
for f=0:N-1
    density(kk)= sum(round(angular_error,1) == round(-180 + f*(360/N),1));
    kk=kk+1;
end

[xData, yData] = prepareCurveData( x, density );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [4 -3.59999999999999 7.18350240095378];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.

h = plot(fitresult,xData, yData);

set(gca,'color',[183 219 255]/255)
% legend( h, name, 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'density', 'Interpreter', 'none' );
grid on
hold on
title('Fit model to data')


angular_error = respMat(:,2)';
trials =size(respMat,1);
N= respMat(1,7);
k=1;

degree=360/N;
for f=0:N-1
    density(k)= sum(round(angular_error,1) == round(-180 + f*(360/N),1));
    k=k+1;
end


reaction_time = respMat(:,4)';

%performance
p=zeros(1,5);
for i=1:trials
    if abs(angular_error(i))<10                                        %excellent performance < 20 degree
       p(1) = p(1)+1 ;
    end
    if abs(angular_error(i))>10     && abs(angular_error(i))<20        %good performance <40 degree
       p(2) = p(2)+1 ;
    end
     if abs(angular_error(i))>20     && abs(angular_error(i))<30        %fine performance <60 degree 
       p(3) = p(3)+1 ;
     end
     if abs(angular_error(i))>30    && abs(angular_error(i))<50         %medium performance <80 degree 
       p(4) = p(4)+1 ;
     end
     if abs(angular_error(i))>50           
       p(5) = p(5)+1 ;
     end
                                                      
end
overall_performance = p(1)/trials*100;
labels = {'excellent performance','good performance','fine performance','medium performance','weak performance'};

%----------------------------------------------------------------------
%                      %plot results
%---------------------------------------------------------------------- 
figure

%bar(x,density)
scatter(x,density,20,hsv(N),'filled');
set(gca,'color',[198 198 198]/255)
grid on
xlim([-200 200]);
ylim([0 max(density)+2]);
xlabel('angular error ')
ylabel('density')


figure
hold on
bar(reaction_time,'facecolor',[193,0,97]/255,'EdgeColor','k')
set(gca,'color',[198 198 198]/255)
grid on
xlabel('trial')
ylabel('reaction time(s)')

figure
hold on
axis off
axis equal
ax = gca(); 
pie(p)
ax.Colormap = summer(5);
set(gca,'color',[198 198 198]/255);
legend(labels,'Location', 'southeastoutside', 'Interpreter', 'none')

figure
hold on
fb = bar(angular_error);
set(gca,'color',[198 198 198]/255)
grid on
xlabel('trial')
ylabel('angular error')
%----------------------------------------------------------------------
%                      %analysis level 3
%----------------------------------------------------------------------
cool = [];
warm = [];
c = 1;
w = 1;
cmp = hsv(N);

for i=1:N
    if cmp(i,2)+cmp(i,3)> cmp(i,1)
       cool(c,:) = 255*cmp(i,:);
       c=c+1;
    elseif cmp(i,2)+cmp(i,3) <= cmp(i,1)
        warm(w,:) = 255*cmp(i,:);
        w=w+1;
    end
end

CWF = zeros(trials,3);

for i = 1:trials
    CWF(i,1) = i;
    if respMat(i,5) == 0
        place = any(cool == 255*cmp(end,:));
        if place ~= 0
            CWF(i,2) = respMat(i,2);    %Showing angular error in cool colors
        else
            CWF(i,3) = respMat(i,2);    %Showing angular error in warm colors
        end
    else
        place = any(cool == 255*cmp(respMat(i,5),:));
        if place ~= 0
            CWF(i,2) = respMat(i,2);    %Showing angular error in cool colors
        else
            CWF(i,3) = respMat(i,2);    %Showing angular error in warm colors
        end
    end
end

Cerr = 0;
Werr = 0;
CC = 0;
WC = 0;

for i=1:trials
    if CWF(i,2) == 0
        WC = WC + 1;
    else
        CC = CC + 1;
    end
    Cerr = Cerr + abs(CWF(i,2));
    Werr = Werr + abs(CWF(i,3));
end

Cerr = Cerr/CC;
Werr = Werr/WC;
figure
hb = bar(reordercats(categorical({'Cool','Warm'}),{'Cool','Warm'}),[Cerr,Werr],'facecolor', 'flat');
hb.CData =[134 156 211 ;223 145 16]/255 ;
set(gca,'color',[183 219 255]/255)
grid on

ylabel('mean angular error')
%----------------------------------------------------------------------
%                      %analysis level 4
%----------------------------------------------------------------------
shape = respMat(:,3);
circle=0;
triangle=0;
angular_error_circle = 0;
angular_error_triangle = 0 ;
for i=1:trials
   if  shape(i)==1
      circle = circle+1; 
      angular_error_circle =abs(angular_error(i))+ angular_error_circle;
   end
   if  shape(i)==2
      triangle = triangle+1; 
      angular_error_triangle =abs( angular_error(i))+ angular_error_triangle;
   end
   
end
angular_error_circle = angular_error_circle/circle;
angular_error_triangle = angular_error_triangle/triangle;

figure
hd = bar(reordercats(categorical({'Circle','Triangle'}),{'Circle','Triangle'}),[angular_error_circle,angular_error_triangle],'facecolor', 'flat')
hd.CData =[178 198 65; 223 145 16]/255 ;
set(gca,'color',[183 219 255]/255)
grid on

ylabel('mean angular error')
%%
%----------------------------------------------------------------------
%                      %analysis level 5
%---------------------------------------------------------------------- 
clear

myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.mat')); %gets all mat files in struct
ResData = cell(length(myFiles),2);
colorfit = hsv(length(myFiles));
for k = 1:length(myFiles)
  ResData{k,1} = myFiles(k).name;
  name = char(ResData{k,1});
  name = erase(name,'VisualWorkingMemory_');
  name = erase(name,'.mat');
  ResData{k,1} = name;
  load(myFiles(k).name);
  ResData{k,2} = respMat;
 



angular_error = respMat(:,2)';
trials =size(respMat,1);
N= respMat(1,7);
kk=1;
degree = 360/N;
x=-180:degree:(180-degree);
degree=360/N;
for f=0:N-1
    density(kk)= sum(round(angular_error,1) == round(-180 + f*(360/N),1));
    kk=kk+1;
end

[xData, yData] = prepareCurveData( x, density );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [4 -3.59999999999999 7.18350240095378];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.

h = plot(fitresult);
h.Color = colorfit(k,:);
h.LineWidth = 1;
set(gca,'color',[183 219 255]/255)
% legend( h, name, 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'density', 'Interpreter', 'none' );
grid on
hold on
end

legend(ResData(:,1)');
