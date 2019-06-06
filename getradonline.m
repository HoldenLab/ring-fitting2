function [R,rho] = getradonline(I,th,minsep)

    figure;
    subplot(2,2,1) ;
    imagesc(I) ;
    title('Original Image') ;

    % Convert image to black white 
    %BW = edge(I,'Sobel');
    BW=im2bw(I,0.25) ;
    subplot(2,2,2) ;
    imshow(BW); 
    title('BW Image') ;

    % Radon transform 
    % Angle projections  
    theta = [-90:90]' ;
    [R, rho] = radon(I, theta) ;
    subplot(2,2,3) ;
    imshow(R, [], 'XData', theta, 'YData', rho, 'InitialMagnification', 'fit');
    xlabel('\theta'), ylabel('\rho');
    axis on, axis normal, hold on;

    % Detect the peaks of transform output  
    % Threshold value for peak detection

    threshold_val = ceil(th*max(R(:))) ;
    % Maximum nof peaks to identify in the image
    max_nofpeaks=inf;
    ispeak = findpeaksn(R,[1 1],round(th*max(R(:))),minsep);
    max_indexes = find(ispeak);
    max_values = R(max_indexes) ;
    [sorted_max, temp_indexes] = sort(max_values, 'descend') ;
    sorted_indexes = max_indexes(temp_indexes) ;

    % Get the first highest peaks for the sorted array
    if (length(sorted_max) <= max_nofpeaks)
        peak_values = sorted_max(1:end) ; 
        peak_indexes = sorted_indexes(1:end) ;
    else
        peak_values = sorted_max(1:max_nofpeaks) ;
        peak_indexes = sorted_indexes(1:max_nofpeaks) ;
    end
    [y, x]  = ind2sub(size(R), peak_indexes ) ;
    peaks = [rho(y) theta(x)] ;
    plot(peaks(:,2), peaks(:,1), 's', 'color','r');
    title('Radon Transform & Peaks') ;

    % Detected lines on the image
    figure
     imagesc(I), title('Detected lines'), hold on

    x_center = floor(size(I, 2)/2) ;
    y_center = floor(size(I, 1)/2) ;
    for p=1:length(peaks)

        x_1 = [2*(-x_center), 2*x_center] ;
        %y_1 = [0, 0] ;
        y_1 = [0, 0] ;

        % Shift at first
        x_1_shifted = x_1 ;
        y_1_shifted = [y_1(1)-peaks(p,1), y_1(2)-peaks(p,1)] ;

        % Rotate 
        peaks(p,2) = 90 - peaks(p,2) ;
        t=peaks(p,2)*pi/180;
        rotation_mat = [ cos(t) -sin(t) ; sin(t) cos(t) ] ;
        x_y_rotated = rotation_mat*[x_1_shifted; y_1_shifted] ;
        x_rotated = x_y_rotated(1,:) ;
        y_rotated = x_y_rotated(2,:) ;
        plot( x_rotated+x_center, y_rotated+y_center, 'k', 'linewidth', 2 );
   end
   hold off;
