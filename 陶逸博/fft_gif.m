fs = 100;
T = 1 / fs;
N = 1024;
M = 30;
t = linspace(-5, 5, N);
f = fs*(0:(N/2))/N;
a = linspace(5, 50, M); % a = V / X0
E = zeros(M, N);
dE = zeros(M, N);
ddE = zeros(M, N);

for i = 1:M
    [E(i,:), dE(i,:), ddE(i,:)] = reference(t, a(i));
end

figure(1)
for i = [linspace(1,M,M) linspace(M-1,2,M-2)]
    subplot(3,2,1)
    plot(t,E(i,:))
    xlim([-5,5])
    ylim([0,1.5])
    title('E(t)')
    
    subplot(3,2,2)
    [mag, ~] = fft_1D(E(i,:));
    plot(f, mag);
    xlim([-2, fs/2+2])
    ylim([0,0.6])
    title('|E(f)|')
    
    subplot(3,2,3)
    plot(t,dE(i,:))
    xlim([-5,5])
    ylim([-8,8])
    title('dE(t)')
    
    subplot(3,2,4)
    [mag, ~] = fft_1D(dE(i,:));
    plot(f, mag)
    xlim([-2, fs/2+2])
    ylim([0,0.3])
    title('|dE(f)|')
    
    subplot(3,2,5)
    plot(t,ddE(i,:))
    xlim([-5,5])
    ylim([-20,20])
    title('ddE(t)')
    
    subplot(3,2,6)
    [mag, ~] = fft_1D(ddE(i,:));
    plot(f, mag)
    xlim([-2, fs/2+2])
    ylim([0,4])
    title('|ddE(f)|')
    
    frame=getframe(gcf);  
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if i == 1
        imwrite(imind,cm,'suidong.gif','GIF', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,'suidong.gif','GIF','WriteMode','append','DelayTime',0.1);
    end  
end

function [E, dE, ddE] = reference(t, a)
    A = atan(a * t);
    E = atan(5 * cos(A));
    
    dA = a * cos(A).^2;
    dE = -5 * sin(A) .* (cos(E).^2) .* dA;
    
    ddE = [diff(dE) ./ diff(t) 0];
end

function [mag, phase] = fft_1D(X)     
    N = length(X);
    assert(mod(N, 2) == 0)
    Y = fft(X);
    P = abs(Y/N);
    mag = P(1:N/2+1);
    mag(2:end-1) = 2 * mag(2:end-1);
    phase = angle(Y);
    phase = phase(1:N/2+1);
end
