% plots
% x asis: # of iterations
% y axis: Infinity norm

Npoints = [50*50,100*100, 500*500, 1000*1000,5000*5000];
Niter = [10,10,10,11,11];

hold on;
loglog(Npoints, Niter,'b-o','LineWidth',2)
%plot([3,5],'x','DisplayName','Garbage','Linewidth',2)
xlabel('Number of Points','Fontsize',12)
ylabel('Number of iterations','Fontsize',12)
grid on
%legend('show')
ylim([0,20])

% % Styling: o = circles, - = solid line, -- = dashed line
% plot(x,y,'bo',var,P,'r--')  % Plot...add as many x and y combinations as needed 
% legend('Data Points','Langrange Approximation') % names in order of plotting above
% 
% loglog , semilogx , semilogy % logarithmic, semi-logarithmic plots
% hold on % between plot calls would create different figures.
% 
% % Plotting with loops.
% hold on
% for c = 1:d
%     len = 5/h(c)+1;
%     plot(T(1:len,c),er_F(1:len,c),'-o','Linewidth',1.75,'DisplayName',['h=' num2str(h(c))])
% end
% title('Absolute Error in Forward Euler with varied grid spacing (h)','FontSize',14)
% xlabel('Time','FontSize',14)
% ylabel('Absolute Error','FontSize',14)
% axis sqaure
% legend('show')