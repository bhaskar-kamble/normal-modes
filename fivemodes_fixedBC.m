%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     %
% Program to visualize the normal modes of a          %
% spring-mass system consisting of 5 masses connected %
% with springs, with fixed boundary conditions.       %
%                                                     %
% The program plots and saves a series of jpeg files  %
% showing the positions of the masses in each normal  %
% mode. The files can then be strung together to create
% an animation.                                       %
%                                                     %
% keywords - waves, oscillations, normal modes        %
%                                                     %
% Author:- Bhaskar Kamble                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;
N=5;                      % Number of masses
L_tab=6.0;                % Length of the table
Dn=zeros(N,N);            % The matrix governing the equations of motion
dy=0.06;
nspr=20;                  % Number of kinks in each spring

M=800;                    % Number of image files that will be produced
T=linspace(0.0,40*pi,M);  % Start Time and Stop Time, with M divisions
T=T(1:M-1);
theta1=zeros(1,M-1);      % Position of 1st mass
theta2=zeros(1,M-1);      % Position of 2nd mass
theta3=zeros(1,M-1);      % Position of 3rd mass
theta4=zeros(1,M-1);      % Position of 4th mass
theta5=zeros(1,M-1);      % Position of 5th mass

% Dn is a tri-diagonal matrix with 2 on the diagonals...
% ...and -1 on the upper and lower diagonals:
for ii=1:N
    Dn(ii,ii)=2;
end;
for ii=1:N-1
    Dn(ii+1,ii)=-1;
    Dn(ii,ii+1)=-1;
end;
    
% Diagonalize Dn:    
[Vect,Val]=eig(Dn);
Eval=diag(Val);
% Eigenvectors (Vect) give the relative amplitudes in the normal modes
% Eigenvalues (Eval) give the squares of the frequencies

% Normalize the normal mode frequencies with the lowest frequency:
wp1=1.0;
wp2=sqrt(Eval(2)/Eval(1));
wp3=sqrt(Eval(3)/Eval(1));
wp4=sqrt(Eval(4)/Eval(1));
wp5=sqrt(Eval(5)/Eval(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%MODE ONE%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Spr=L_tab/(N+1);%inter-mass separation
Eqb1=-(L_tab*0.5)+Spr;
Eqb2=Eqb1+Spr;
Eqb3=Eqb2+Spr;
Eqb4=Eqb3+Spr;
Eqb5=Eqb4+Spr;

for iter=1:length(T)
    figure(1);
    clf;
    hold on;
    A=Vect(:,1);
    theta1(iter)=A(1)*cos(wp1*T(iter))+Eqb1;
    theta2(iter)=A(2)*cos(wp1*T(iter))+Eqb2;
    theta3(iter)=A(3)*cos(wp1*T(iter))+Eqb3;
    theta4(iter)=A(4)*cos(wp1*T(iter))+Eqb4;
    theta5(iter)=A(5)*cos(wp1*T(iter))+Eqb5;
    axis off;
    %the spring:
    dspr01=(theta1(iter)+(0.5*L_tab))/nspr;
    dspr12=(theta2(iter)-theta1(iter))/nspr;
    dspr23=(theta3(iter)-theta2(iter))/nspr;
    dspr34=(theta4(iter)-theta3(iter))/nspr;
    dspr45=(theta5(iter)-theta4(iter))/nspr;
    dspr56=((0.5*L_tab)-theta5(iter))/nspr;
    
    for ispr=1:nspr
        bias=2.0;
        plot([-(0.5*L_tab)+ispr*dspr01-0.1,-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1,-(0.5*L_tab)+(ispr+1)*dspr01-0.1],[dy+bias,-dy+bias],'k-');
        hold on;

        plot([theta1(iter)+dspr12*ispr,theta1(iter)+(ispr+0.5)*dspr12],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta1(iter)+(ispr+0.5)*dspr12,theta1(iter)+(ispr+1)*dspr12],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta2(iter)+dspr23*ispr,theta2(iter)+(ispr+0.5)*dspr23],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta2(iter)+(ispr+0.5)*dspr23,theta2(iter)+(ispr+1)*dspr23],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta3(iter)+dspr34*ispr,theta3(iter)+(ispr+0.5)*dspr34],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta3(iter)+(ispr+0.5)*dspr34,theta3(iter)+(ispr+1)*dspr34],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta4(iter)+dspr45*ispr,theta4(iter)+(ispr+0.5)*dspr45],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta4(iter)+(ispr+0.5)*dspr45,theta4(iter)+(ispr+1)*dspr45],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta5(iter)+dspr56*ispr,theta5(iter)+(ispr+0.5)*dspr56],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta5(iter)+(ispr+0.5)*dspr56,theta5(iter)+(ispr+1)*dspr56],[dy+bias,-dy+bias],'k-');
        hold on;
    end;
    plot([Eqb1,Eqb1],[-5,5],'k-.');
    hold on;
    plot([Eqb2,Eqb2],[-5,5],'k-.');
    hold on;
    plot([Eqb3,Eqb3],[-5,5],'k-.');
    hold on;
    plot([Eqb4,Eqb4],[-5,5],'k-.');
    hold on;
    plot([Eqb5,Eqb5],[-5,5],'k-.');
    hold on;
    
    xtrig=-(0.5*L_tab):0.1:(0.5*L_tab);
    thetatrig=xtrig;
    xtrig=xtrig+(0.5*L_tab);
    xtrig=xtrig*(pi)/L_tab;
    trig=0.5*sin(xtrig)+bias;
    plot(thetatrig,trig,'k-')
    hold on;
    
    plot([theta1(iter),theta1(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta2(iter),theta2(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta3(iter),theta3(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;    
    plot([theta4(iter),theta4(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta5(iter),theta5(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2 %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A=Vect(:,2);
    theta1(iter)=A(1)*cos(wp2*T(iter))+Eqb1;
    theta2(iter)=A(2)*cos(wp2*T(iter))+Eqb2;
    theta3(iter)=A(3)*cos(wp2*T(iter))+Eqb3;
    theta4(iter)=A(4)*cos(wp2*T(iter))+Eqb4;
    theta5(iter)=A(5)*cos(wp2*T(iter))+Eqb5;
    %the spring:
    dspr01=(theta1(iter)+(0.5*L_tab))/nspr;
    dspr12=(theta2(iter)-theta1(iter))/nspr;
    dspr23=(theta3(iter)-theta2(iter))/nspr;
    dspr34=(theta4(iter)-theta3(iter))/nspr;
    dspr45=(theta5(iter)-theta4(iter))/nspr;
    dspr56=((0.5*L_tab)-theta5(iter))/nspr;
    
    for ispr=1:nspr
        bias=1;        
        plot([-(0.5*L_tab)+ispr*dspr01-0.1,-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1,-(0.5*L_tab)+(ispr+1)*dspr01-0.1],[dy+bias,-dy+bias],'k-');
        hold on;

        plot([theta1(iter)+dspr12*ispr,theta1(iter)+(ispr+0.5)*dspr12],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta1(iter)+(ispr+0.5)*dspr12,theta1(iter)+(ispr+1)*dspr12],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta2(iter)+dspr23*ispr,theta2(iter)+(ispr+0.5)*dspr23],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta2(iter)+(ispr+0.5)*dspr23,theta2(iter)+(ispr+1)*dspr23],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta3(iter)+dspr34*ispr,theta3(iter)+(ispr+0.5)*dspr34],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta3(iter)+(ispr+0.5)*dspr34,theta3(iter)+(ispr+1)*dspr34],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta4(iter)+dspr45*ispr,theta4(iter)+(ispr+0.5)*dspr45],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta4(iter)+(ispr+0.5)*dspr45,theta4(iter)+(ispr+1)*dspr45],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta5(iter)+dspr56*ispr,theta5(iter)+(ispr+0.5)*dspr56],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta5(iter)+(ispr+0.5)*dspr56,theta5(iter)+(ispr+1)*dspr56],[dy+bias,-dy+bias],'k-');
        hold on;
    end;
    
    xtrig=-(0.5*L_tab):0.1:(0.5*L_tab);
    thetatrig=xtrig;
    xtrig=xtrig+(0.5*L_tab);
    xtrig=xtrig*(2*pi)/L_tab;
    trig=0.5*sin(xtrig)+bias;
    plot(thetatrig,trig,'k-')
    hold on;
    
    plot([theta1(iter),theta1(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta2(iter),theta2(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta3(iter),theta3(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;    
    plot([theta4(iter),theta4(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta5(iter),theta5(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3 %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bias=0;
    A=0.7*Vect(:,3);
    theta1(iter)=A(1)*cos(wp3*T(iter))+Eqb1;
    theta2(iter)=A(2)*cos(wp3*T(iter))+Eqb2;
    theta3(iter)=A(3)*cos(wp3*T(iter))+Eqb3;
    theta4(iter)=A(4)*cos(wp3*T(iter))+Eqb4;
    theta5(iter)=A(5)*cos(wp3*T(iter))+Eqb5;
    %the spring:
    dspr01=(theta1(iter)+(0.5*L_tab))/nspr;
    dspr12=(theta2(iter)-theta1(iter))/nspr;
    dspr23=(theta3(iter)-theta2(iter))/nspr;
    dspr34=(theta4(iter)-theta3(iter))/nspr;
    dspr45=(theta5(iter)-theta4(iter))/nspr;
    dspr56=((0.5*L_tab)-theta5(iter))/nspr;
    for ispr=1:nspr
        
        plot([-(0.5*L_tab)+ispr*dspr01-0.1,-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1,-(0.5*L_tab)+(ispr+1)*dspr01-0.1],[dy+bias,-dy+bias],'k-');
        hold on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot([theta1(iter)+dspr12*ispr,theta1(iter)+(ispr+0.5)*dspr12],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta1(iter)+(ispr+0.5)*dspr12,theta1(iter)+(ispr+1)*dspr12],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta2(iter)+dspr23*ispr,theta2(iter)+(ispr+0.5)*dspr23],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta2(iter)+(ispr+0.5)*dspr23,theta2(iter)+(ispr+1)*dspr23],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta3(iter)+dspr34*ispr,theta3(iter)+(ispr+0.5)*dspr34],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta3(iter)+(ispr+0.5)*dspr34,theta3(iter)+(ispr+1)*dspr34],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta4(iter)+dspr45*ispr,theta4(iter)+(ispr+0.5)*dspr45],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta4(iter)+(ispr+0.5)*dspr45,theta4(iter)+(ispr+1)*dspr45],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta5(iter)+dspr56*ispr,theta5(iter)+(ispr+0.5)*dspr56],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta5(iter)+(ispr+0.5)*dspr56,theta5(iter)+(ispr+1)*dspr56],[dy+bias,-dy+bias],'k-');
        hold on;
    end;
    
    xtrig=-(0.5*L_tab):0.1:(0.5*L_tab);
    thetatrig=xtrig;
    xtrig=xtrig+(0.5*L_tab);
    xtrig=xtrig*(3*pi)/L_tab;
    trig=0.5*sin(xtrig)+bias;
    plot(thetatrig,trig,'k-')
    hold on;
    
    plot([theta1(iter),theta1(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta2(iter),theta2(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta3(iter),theta3(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;    
    plot([theta4(iter),theta4(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta5(iter),theta5(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4 %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bias=-1;
    A=0.5*Vect(:,4);
    theta1(iter)=A(1)*cos(wp4*T(iter))+Eqb1;
    theta2(iter)=A(2)*cos(wp4*T(iter))+Eqb2;
    theta3(iter)=A(3)*cos(wp4*T(iter))+Eqb3;
    theta4(iter)=A(4)*cos(wp4*T(iter))+Eqb4;
    theta5(iter)=A(5)*cos(wp4*T(iter))+Eqb5;
    %the spring:
    dspr01=(theta1(iter)+(0.5*L_tab))/nspr;
    dspr12=(theta2(iter)-theta1(iter))/nspr;
    dspr23=(theta3(iter)-theta2(iter))/nspr;
    dspr34=(theta4(iter)-theta3(iter))/nspr;
    dspr45=(theta5(iter)-theta4(iter))/nspr;
    dspr56=((0.5*L_tab)-theta5(iter))/nspr;
    
    for ispr=1:nspr
        
        plot([-(0.5*L_tab)+ispr*dspr01-0.1,-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1,-(0.5*L_tab)+(ispr+1)*dspr01-0.1],[dy+bias,-dy+bias],'k-');
        hold on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot([theta1(iter)+dspr12*ispr,theta1(iter)+(ispr+0.5)*dspr12],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta1(iter)+(ispr+0.5)*dspr12,theta1(iter)+(ispr+1)*dspr12],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta2(iter)+dspr23*ispr,theta2(iter)+(ispr+0.5)*dspr23],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta2(iter)+(ispr+0.5)*dspr23,theta2(iter)+(ispr+1)*dspr23],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta3(iter)+dspr34*ispr,theta3(iter)+(ispr+0.5)*dspr34],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta3(iter)+(ispr+0.5)*dspr34,theta3(iter)+(ispr+1)*dspr34],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta4(iter)+dspr45*ispr,theta4(iter)+(ispr+0.5)*dspr45],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta4(iter)+(ispr+0.5)*dspr45,theta4(iter)+(ispr+1)*dspr45],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta5(iter)+dspr56*ispr,theta5(iter)+(ispr+0.5)*dspr56],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta5(iter)+(ispr+0.5)*dspr56,theta5(iter)+(ispr+1)*dspr56],[dy+bias,-dy+bias],'k-');
        hold on;
    end;
    
    xtrig=-(0.5*L_tab):0.1:(0.5*L_tab);
    thetatrig=xtrig;
    xtrig=xtrig+(0.5*L_tab);
    xtrig=xtrig*(4*pi)/L_tab;
    trig=0.5*sin(xtrig)+bias;
    plot(thetatrig,trig,'k-')
    hold on;
    
    plot([theta1(iter),theta1(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta2(iter),theta2(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta3(iter),theta3(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;    
    plot([theta4(iter),theta4(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta5(iter),theta5(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5 %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bias=-2;
    A=0.4*Vect(:,5);
    theta1(iter)=A(1)*cos(wp5*T(iter))+Eqb1;
    theta2(iter)=A(2)*cos(wp5*T(iter))+Eqb2;
    theta3(iter)=A(3)*cos(wp5*T(iter))+Eqb3;
    theta4(iter)=A(4)*cos(wp5*T(iter))+Eqb4;
    theta5(iter)=A(5)*cos(wp5*T(iter))+Eqb5;
    %the spring:
    dspr01=(theta1(iter)+(0.5*L_tab))/nspr;
    dspr12=(theta2(iter)-theta1(iter))/nspr;
    dspr23=(theta3(iter)-theta2(iter))/nspr;
    dspr34=(theta4(iter)-theta3(iter))/nspr;
    dspr45=(theta5(iter)-theta4(iter))/nspr;
    dspr56=((0.5*L_tab)-theta5(iter))/nspr;
    
    for ispr=1:nspr
        
        plot([-(0.5*L_tab)+ispr*dspr01-0.1,-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([-(0.5*L_tab)+(ispr+0.5)*dspr01-0.1,-(0.5*L_tab)+(ispr+1)*dspr01-0.1],[dy+bias,-dy+bias],'k-');
        hold on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot([theta1(iter)+dspr12*ispr,theta1(iter)+(ispr+0.5)*dspr12],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta1(iter)+(ispr+0.5)*dspr12,theta1(iter)+(ispr+1)*dspr12],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta2(iter)+dspr23*ispr,theta2(iter)+(ispr+0.5)*dspr23],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta2(iter)+(ispr+0.5)*dspr23,theta2(iter)+(ispr+1)*dspr23],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta3(iter)+dspr34*ispr,theta3(iter)+(ispr+0.5)*dspr34],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta3(iter)+(ispr+0.5)*dspr34,theta3(iter)+(ispr+1)*dspr34],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta4(iter)+dspr45*ispr,theta4(iter)+(ispr+0.5)*dspr45],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta4(iter)+(ispr+0.5)*dspr45,theta4(iter)+(ispr+1)*dspr45],[dy+bias,-dy+bias],'k-');
        hold on;
        
        plot([theta5(iter)+dspr56*ispr,theta5(iter)+(ispr+0.5)*dspr56],[-dy+bias,dy+bias],'k-');
        hold on;
        plot([theta5(iter)+(ispr+0.5)*dspr56,theta5(iter)+(ispr+1)*dspr56],[dy+bias,-dy+bias],'k-');
        hold on;
    end;
    
    xtrig=-(0.5*L_tab):0.1:(0.5*L_tab);
    thetatrig=xtrig;
    xtrig=xtrig+(0.5*L_tab);
    xtrig=xtrig*(5*pi)/L_tab;
    trig=0.5*sin(xtrig)+bias;
    plot(thetatrig,trig,'k-')
    hold on;
    
    plot([theta1(iter),theta1(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta2(iter),theta2(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta3(iter),theta3(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;    
    plot([theta4(iter),theta4(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    plot([theta5(iter),theta5(iter)],[bias,bias],'rs','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','r');
    hold on;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    axis([-3,3,-3,3]); %xmin xmax ymin ymax
    file=sprintf('Frame%d.jpg', 1000+iter);
    disp(sprintf('Saving to %s', file));
    print('-dpng',  '-zbuffer',  '-r100', file);
    %pause(0.1);
end;
