  %===matlab script to obtain relaxation sequence===
% 1. N is the mesh size.
% 2. Input is W, the vector that contains distinct relaxation factors and 
%    Q, the iteration counts of each relaxation factor in one cycle.
% 3. Output is a vector wt of length M, it is the sequence of relaxation
%    factor in one cycle.
% 4. All the scheme except P=2 schemes has been put in. The user can just 
%    uncomment the corresponding line of W and Q. DO NOT uncomment the 
%    line before W. 
% 5. Q vector is printed to show how the cycle is being used. 

close all
clear all
format long

N=512;                      %mesh size: N by N  ===input===
kmin=sin(pi/2/N)^2;         %kappa_min

k=[kmin:kmin:kmin]';
k=[k;[k(end)+kmin:kmin:2]'];


%for rho>100==========================================
%rho=148 N=512;-------------
%W=[91299,25979,3862.1,549.9,80.217,11.992,1.9595,0.59145];
%Q=[1,3,9,27,81,243,729,1337];
%rho=90 N=1024;-------------
%W=[178919, 8024.1, 349.03, 15.9047, 0.799909];
%Q=[1, 7, 49, 343, 3087];
%rho=90 N=2048;-------------
%W=[1664620, 894851,66205.7, 2623.72,99.9385,0.936974];
%Q=[1,8,64,512,4096,200119];

%=====================for P=3 Schemes=================
%P=3 N=16-------------
W=[64.66,6.215,0.7042];
Q=[1 5 21];
%P=3 N=32-------------
%W=[213.8,11.45,0.7616];
%Q=[1 7 45];
%P=3 N=64-------------
%W=[684.3,20.73,0.8149];
%Q=[1 11 106];
%P=3 N=128-------------
% W=[2114,36.78,0.8611];
% Q=[1 17 252];
%P=3 N=256-------------
%W=[6319,63.99,0.8989];
%Q=[1 27 625];
%P=3 N=512-------------
%W=[18278,109.2,0.9282];
%Q=[1 43 1571];
%P=3 N=1024-----------
%W=[ 51769.1,184.31,0.95025];
%Q=[1 68 3955];

%=====================for P=4 Schemes=================
%P=4 N=16-------------
%W=[80.154, 17.217, 2.6201, 0.62230];
%Q=[1 2 8 20];
%P=4 N=32-------------
%W=[289.46 ,40.791, 4.0877, 0.66277];
%Q=[1 3 14 46];
%P=4 N=64-------------
%W=[1029.4, 95.007, 6.3913, 0.70513];
%Q=[1 5 26 114];
%P=4 N=128-------------
%W=[3596.4, 217.80, 9.9666, 0.74755];
%Q=[1  7  50  285];
%P=4 N=256-------------
%W=[12329, 492.05, 15.444, 0.78831];
%Q=[1,9,86,664];
%P=4 N=512-------------
%W=[41459, 1096.3, 23.730, 0.82597];
%Q=[1,12,155,1650];

%=====================for P=5 Schemes=================
%P=5 N=16-------------
%W=[88.190,30.122,6.8843,1.6008,0.58003];
%Q=[1,2,5,12,23];
%P=5 N=32-------------
%W=[330.57,82.172,13.441,2.2402,0.60810];
%Q=[1, 2, 7, 20, 46];
%P=5 N=64-------------
%W=[1228.8,220.14,26.168,3.1668,0.63890];
%Q=[1, 3, 10, 38, 106];
%P=5 N=128-------------
%W=[4522.0,580.86,50.729,4.5018,0.67161];
%Q=[1,3,16,73,250];
%P=5 N=256-------------
%W=[16459,1513.4,97.832,6.4111,0.70531];
%Q=[1,4,26,142,605];
%P=5 N=512-------------
%W=[59226,3900.56,187.53,9.1194,0.73905];
%Q=[1,6,40,277,1500];


M=sum(Q);
wt=zeros(1,sum(Q));
G=ones(size(k));

wt(1)=W(1);
Q(1)=Q(1)-1;
G=G.*abs(1-k*wt(1));
counter=2;
totPerc=100;
while(sum(Q~=0))
    if(sum(Q==0)~=0)
        index=find(Q==0);
        if(index==1)
            W=W(2:end);
            Q=Q(2:end);
        elseif(index==length(Q))
            W=W(1:end-1);
            Q=Q(1:end-1);
        else
            W(index:end-1)=W(index+1:end);
            W=W(1:end-1);
            Q(index:end-1)=Q(index+1:end);
            Q=Q(1:end-1);
        end
    end
    Q
    ww=1/k(find(G==max(G)));
    dis=abs(W-ww);
    index=find(dis==min(dis));
    wt(counter)=W(index);
    G=G.*abs(1-k*wt(counter));
    Q(index)=Q(index)-1;
    counter=counter+1;
    perc=sum(Q)/length(wt);
end
save wt.mat wt M

h=figure();
set(h,'units','points');
set(h,'Color',[1 1 1],'InvertHardcopy','off','position',[50 100 347*0.8 247*0.8]);
set(h,'DefaultAxesFontName', 'Times New Roman')
set(h,'DefaultTextFontname', 'Times New Roman')
set(h,'DefaultAxesFontSize', 11)
set(h,'DefaultTextFontSize', 11)
plot(wt);
xlabel('Iteration');
ylabel('\omega');
set(gca,'yscale','log')
set(gca,'ytick',[0.1 1 10 100 1000 10000 100000])