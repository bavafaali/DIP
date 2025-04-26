clearvars, clc, close all, clear all;

% reading the image and extracting green channel
I_orig = imread('02_test.tif');
I = I_orig(:,:,2);
I = imgaussfilt(I,0.463);

%Thin vessel enhancement
sigma = 20;
Se1_size = 4;
Se2_size = 24;
Thin_TH = Optimized_tophat(Se1_size,Se2_size,I);
Thin_Hf = HomoFilt(sigma,Thin_TH);
Thin_2mfilt = twoMatchedFilter(Thin_Hf);
Thin_final = mcet_hho(Thin_2mfilt);


%Thick vessel enhancement
sigma = 2;
Se1_size = 8;
Se2_size = 16;
Thick_OTH = Optimized_tophat(Se1_size,Se2_size,I);
Thick_TH = tophat(Se1_size,I);
Thick_Hf = HomoFilt(sigma,Thick_OTH);
Thick_median = medfilt2(Thick_Hf,[2 2]);
Thick_final = Optimized_tophat(Se1_size,Se2_size,Thick_median);

Final = Thick_final + Thin_final;


figure(1)
subplot(1,3,1)
imshow(I_orig,[])
title('Original Image', 'FontSize', 13)
axis off
subplot(1,3,2)
imshow(I,[])
title('Green channel', 'FontSize', 13)
axis off
subplot(1,3,3)
imshow(Final,[])
title('Enhancement', 'FontSize', 13)
axis off

figure(2)
subplot(1,2,1)
imshow(Thin_Hf,[])
title('Thin vessel enhancement', 'FontSize', 13)
axis off
subplot(1,2,2)
imshow(Thick_Hf,[])
title('Thick vessel enhancement', 'FontSize', 13)
axis off

figure(3)
subplot(1,2,1)
imshow(Thick_OTH,[])
title('Optimized top hat', 'FontSize', 13)
axis off
subplot(1,2,2)
imshow(Thick_TH,[])
title('Top hat', 'FontSize', 13)
axis off

canny = edge(I,'Canny');
sobel = edge(I,'sobel');
prewitt = edge(I,'prewitt');
lev = multithresh(I, 4);
otsu = imquantize(I,lev);
T = adaptthresh(I, 0.5); % The sensitivity factor is set to 0.5
adaptiveTH = imbinarize(I, T);

figure(4)
subplot(2,5,3)
imshow(Final,[])
title('Paper', 'FontSize', 13)
axis off
subplot(2,5,6)
imshow(canny,[])
title('Canny', 'FontSize', 13)
axis off
subplot(2,5,7)
imshow(sobel,[])
title('Sobel', 'FontSize', 13)
axis off
subplot(2,5,8)
imshow(prewitt,[])
title('Prewitt', 'FontSize', 13)
axis off
subplot(2,5,9)
imshow(adaptiveTH,[])
title('adaptive threshold', 'FontSize', 13)
axis off
subplot(2,5,10)
imshow(otsu,[])
title('Otsu', 'FontSize', 13)
axis off


function out = twoMatchedFilter(image)
    img = double(image);
   
    % Parameters
    s = 0.8;
    L = 7;
    theta = 0:7:182;
    
    out = zeros(size(img));
    
    m = max(ceil(3*s),(L-1)/2);
    [x,y] = meshgrid(-m:m,-m:m); % non-rotated coordinate system, contains (0,0)
    for t = theta
       t = t / 180 * pi;        % angle in radian
       u = cos(t)*x - sin(t)*y; % rotated coordinate system
       v = sin(t)*x + cos(t)*y; % rotated coordinate system
       N = (abs(u) <= 3*s) & (abs(v) <= L/2); % domain
       k = exp(-u.^2/(2*s.^2)); % kernel
       k = k - mean(k(N));
       k(~N) = 0;               % set kernel outside of domain to 0
    
       res = conv2(img,k,'same');
       out = max(out,res);
    end
    
    out = uint8(out); % force output to be in [0,1] interval that MATLAB likes
end

function result = mcet_hho(I)
    [h,nh]=imhist(I);           % Get Histogram
    [m,n]=size(I);              % Image size
    L=length(h);                % Lmax levels to segment 0 - 256
    Nt=size(I,1) * size(I,2);   % Total pixels in the image
    
    % Frequency distribution of each intensity level of the histogram 0 - 256
    for i=1:L 
        probI(i)=h(i)/Nt;
    end

    %% Initial data of the HHO algorithm
    nVar=4;                                 % Number of thresholds (Th)     
    VarSize=[1 nVar];                       % Decision Variables in Matrix
    VarMin=1;                               % Minimum value of Th
    VarMax=255;                             % Maximum value of Th
    %% Harris Hawks Algorithm Parameters
    N=30;                                   % Maximum Number of Hawks
    T=250;                                  % Maximum Number of Iterations
    
    tic
    Rabbit_Location = zeros(1,nVar);          % Initialization of the rabbit's location
    Rabbit_Energy=inf;                      % Initialization of the energy of the rabbit
    
    %% Initialization of the position of the hawks
    
    X=initialization(N,nVar,VarMax,VarMin);
    
    %% Harris Hawks Algorithm Main
    CNVG=zeros(1,T);
    t=0;                                    % Counter
    
    while t<T
        for i=1:size(X,1)
            % Check bounds
            FU=X(i,:)>=VarMax;
            FL=X(i,:)<=VarMin;
            X(i,:)=sort(round((X(i,:).*(~(FU+FL)))+VarMax.*FU+VarMin.*FL));
            % Calculate the fitness for each hawk
            fitness=M_CEM(I,X(i,:),h);
            % Update rabbit location with the best fitness
            if fitness<Rabbit_Energy
                Rabbit_Energy=fitness;
                Rabbit_Location=X(i,:);
            end
        end
        
        E1=2*(1-(t/T)); % factor to show the decreaing energy of rabbit
        % Update the location of Harris' hawks
        for i=1:size(X,1)
            E0=2*rand()-1; %-1<E0<1
            Escaping_Energy=E1*(E0);  % escaping energy of rabbit
            
            if abs(Escaping_Energy)>=1
                %% Exploration:
                
                q=rand();
                rand_Hawk_index = floor(N*rand()+1);
                X_rand = X(rand_Hawk_index, :);
                if q<0.5
                    % perch based on other family members
                    X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
                elseif q>=0.5
                    % perch on a random tall tree (random site inside group's home range)
                    X(i,:)=round((Rabbit_Location(1,:)-mean(X))-rand()*((VarMax-VarMin)*rand+VarMin));
                end
                
            elseif abs(Escaping_Energy)<1
                                                   
                %% Exploitation:
                
                %% phase 1: surprise pounce (seven kills)
                % surprise pounce (seven kills): multiple, short rapid dives by different hawks
                
                r=rand(); % probablity of each event
                
                if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                    X(i,:)=abs((Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:)));
                end
                
                if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                    Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                    X(i,:)=abs((Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:)));
                end
                
                %% phase 2: performing team rapid dives (leapfrog movements)
                if r<0.5 && abs(Escaping_Energy)>=0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                    
                    %Jump_strength=2*(1-rand());
                    Jump_strength=(1-rand());
                    X1=abs(round(Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))));
                    for k=1:length(X1)
                        if X1(:,k)>VarMax
                            X1(:,k)=255;
                        end
                        if X1(:,k)<VarMin
                            X1(:,k)=1;
                        end 
                    end
                    if M_CEM(I,X1,h)<M_CEM(I,X(i,:),h) % improved move?
                        X(i,:)=X1;
                    else % hawks perform levy-based short rapid dives around the rabbit
                        X2=abs(round(Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,nVar).*Levy(nVar)));
                        for k=1:length(X2)
                            if X2(:,k)>VarMax
                                X2(:,k)=255;
                            end
                            if X2(:,k)<VarMin
                                X2(:,k)=1;
                            end 
                        end
                        if M_CEM(I,X2,h)<M_CEM(I,X(i,:),h) % improved move?
                            X(i,:)=X2;
                        end
                    end
                end
                
                if r<0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                    % hawks try to decrease their average location with the rabbit
                    Jump_strength=2*(1-rand());
                    X1=abs(round(Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))));
                    for k=1:length(X1)
                        if X1(:,k)>VarMax
                            X1(:,k)=255;
                        end  
                        if X1(:,k)<VarMin
                            X1(:,k)=1;
                        end 
                    end
                    if M_CEM(I,X1,h)<M_CEM(I,X(i,:),h) % improved move?
                        X(i,:)=X1;
                    else % Perform levy-based short rapid dives around the rabbit
                        X2=abs(round(Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,nVar).*Levy(nVar)));
                        for k=1:length(X2)
                            if X2(:,k)>VarMax
                            X2(:,k)=255;
                            end 
                            if X2(:,k)<VarMin
                            X2(:,k)=1;
                            end   
                        end
                        if M_CEM(I,X2,h)<M_CEM(I,X(i,:),h) % improved move?
                            X(i,:)=X2;
                        end
                    end
                end
            end
        end
        t=t+1;
        CNVG(t)=Rabbit_Energy;
    end

    result = MultiTresh(I,Rabbit_Location);
end
