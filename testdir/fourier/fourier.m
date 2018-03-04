tt = linspace(-10,10,1e4);
freq1 = 3; % =omega/2*pi
f = 1+cos(2*pi*freq1*tt);

freq2 = 1;
th = 2*pi*freq2*tt;
r = f;
x = r.*cos(th); y = r.*sin(th);

% there are two frequencies: the frequency of the signal (omega =
% 2*pi/T) and the frequency with which function f is wrapped around in a
% circle (freqw, winding frequency).

figure('name','Cartesian Curve');
hold on; grid on; xlim([0,1]);
xlabel('Time','FontSize',12); ylabel('y=f(t)'); pbaspect([1,0.3,1]);
plot(tt,f, 'r-','LineWidth',1.75,'DisplayName','some wave');
%legend('show');

figure('name','Complex Cuurve');
hold on; grid on;
xlabel('Real'); ylabel('Imaginary');
plot(f.*exp(-2*pi*i*freq2*tt),'r-','LineWidth',1.75); pbaspect([1,1,1]);