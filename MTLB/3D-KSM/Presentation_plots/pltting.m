
Npoints = 0:10;
res = [1.612997867405e+07,2.609001333090e+04,4.645965947371e+03,1.036037761775e+02,4.640066512186e+00,4.297555704406e-01,6.644053321832e-03,2.738406251049e-04,1.721200196660e-05,6.475705996449e-07,1.141988294757e-07];

hold on;
semilogy(Npoints, res,'LineWidth',2)
xlabel('Iterations','Fontsize',15)
ylabel('L2 norm of residual','Fontsize',15)
grid on

%legend('show')
%ylim([0,20])

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