%%
%%%%%%%%%%%%%%%%%%%%%%%%a!=1 not alpha fair with barrier method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all

%%
  

n=1000; %number of iterations
T1=100; %barrier t
eta=0; %small number for ill to make 1/0 feasible in f, gradient and hessian
stop=10e-4; % stop if variables not change in the last following iterations
beta=0.5;
zeta=0.01;
%lambda=0.1; %line search
demand=0.9;
no=1; % number of times experiemtn is computed

%a=[0 0.5 2 10 100 1000 10000];  %alpha=[0-inf]
a=[10];  %alpha=[0-inf]

%N=[10,100,1000,10000,100000,1000000];
N=[1000];
%N=[10,20,30,40,50];
%N=[10];
%s=[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1]; % supplya=[0 0.5 2 3 5 10 20 50 100 200 500 1000 5000 10000];  %alpha=[0-inf]

time = containers.Map; %keyvalue pair data structure
x_values=containers.Map;

for k=1:length(N)
    
    for z=1:length(a)
        disp('%%%%%%%%%%%%%%%%% k %%%%%%%%%%%%%%%%%')
        k
        alpha=a(z); %alpha fairness
        m=N(k);  %number of equations/ variables
        
        A=ones(1,m);
        
        xx_newton=zeros(m,no);
        xx_cvx=zeros(m,no);
        
        for y=1:no
            f =@fun7;%function
            Df=@fun8; %gradient
            Hf=@fun9; %hessian
            inv=@fun6; %inverse
            met='newton'; % method name
            %l=[zeros(m,1);-100000];
            l=zeros(m,1);
            u=10+((100-10)*rand(m,1));
            u=round(u,1);
            %u=10+(100-10)*rand(m+1,1); %random between 10-100

            B=demand*sum(u); % 90% supply

            x=demand*u;
            
            x(m,1)=B-sum(x(1:m-1,1));

             if sum(x)==B
                 disp('initial point is feasible')
                 %break
             else
                 
                 disp('initial point not feasible')
                 %break
             end
            %break;
            
            x_old=zeros(n,m); % to save all solutions
            
            %break
            start_time = tic;

            for i = 1: n

              x_old(i,:)=x;

              disp(i)
              val=f(x,m,B,alpha,u,T1,eta);
              disp('fun done i=')
              z1=Df(x,m,alpha,u,T1,eta);
              disp('Df done i=')
              z2=Hf(x,m,alpha,u,T1,eta);
              disp('Hf done i=')

              %sol = -[diag(z2) A';A 0]\[Df(x,m,alpha,u,T1,eta) ;0];
              sol=zeros(m+1,1);
              [z11,c11,out]=fun5(z2,m);
              
              for p=1:m+1
                disp(p)
                sol(p,1) = -inv(z2,m,p,z11,c11,out)*[Df(x,m,alpha,u,T1,eta) ;0];
              end
              
              disp('sol done i=')
              
              %sol_inv=[diag(z2) A';A 0]*z3;
              disp('solinv done i=')
              v=sol(1:m);
              fprime=z1'*v;

              if (abs(fprime)<stop)
                  break;
              end
              
              lambda=1;
              while (min(x+lambda*v) <=0)
                  lambda=(beta*lambda)+eta
              end
              %xnew=x+lambda*v;
              
              va=f(x+lambda*v,m,B,alpha,u,T1,eta);
              while (va >= val+lambda*zeta*fprime)
                  lambda=(beta*lambda)+eta
              end

              x=x+lambda*v;

              fprintf('value of x: %d\n', x);
              fprintf('value of function: %d\n', val);

            end
            

             endtime=toc(start_time);
             key=strcat('time:: met:',num2str(met),' no equations:',num2str(m),' alpha:',num2str(alpha),' y:',num2str(y));
             time(key)=endtime;
             
             key=strcat('funtion:: met:',num2str(met),' no equations:',num2str(m),' alpha:',num2str(alpha),' y:',num2str(y));
             x_values(key)=val;
             x_newton=x;
             xx_newton(:,y)=x;   
                
             %%%%%%%%%%%%%%%%%%%% cvx method%%%%%%%%%%%%%%%%%%
            
            start_time = tic; 
            met='cvx'; 

            z=ones(1,m)/(1-alpha); %weights + denominator

             cvx_begin
             variable x(m);
             maximize( z*(pow_p(x,1-alpha)) );
             subject to
               x <= u(1:m,1);
               x >= l(1:m,1);
               ones(1,m)*x == B;
             cvx_end

             endtime=toc(start_time);
             key=strcat('time:: met:',num2str(met),' no equations:',num2str(m),' alpha:',num2str(alpha),' y:',num2str(y));
             time(key)=endtime;
             x_cvx=x;
             xx_cvx(:,y)=x;
             
             key=strcat('norm:: no equations:',num2str(m),' alpha:',num2str(alpha),' upper bound',' y:',num2str(y));
             x_values(key)=(norm(x_cvx-x_newton))/m;
             
             key=strcat('funtion:: met:',num2str(met),' no equations:',num2str(m),' alpha:',num2str(alpha),' y:',num2str(y));
             x_values(key)=cvx_optval;
         
        end
        key=strcat('x cvx:: no equations:',num2str(m),' alpha:',num2str(alpha),' upper bound');
        x_values(key)=xx_cvx;
        
        key=strcat('x newton:: no equations:',num2str(m),' alpha:',num2str(alpha),' upper bound');
        x_values(key)=xx_newton;
        
    end
end


%% plots
%%%%%%%%%%%%%%%%%%%%% ploting results %%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc 
close all

%%
load('G:\My Drive\fair electricity distribution\fair allocation of electricity\matlab code\time.mat');
load('G:\My Drive\fair electricity distribution\fair allocation of electricity\matlab code\x_values.mat');

t_newton=zeros(no,length(N));
t_cvx=zeros(no,length(N));
f_diff=zeros(no,length(N));
f_newton=zeros(no,length(N));
f_cvx=zeros(no,length(N));
x_diff=zeros(no,length(N));
for k=1:length(N)
    for j=1:length(a)
        
        alpha=a(j); %alpha fairness
        m=N(k);  %number of equations/ variables
        
        for y=1:no
            
            met='newton'; % method name
            key=strcat('time:: met:',num2str(met),' no equations:',num2str(m),' alpha:',num2str(alpha),' y:',num2str(y));
            %key=strcat('time:: met:',num2str(met),' no equations:',num2str(m),' alpha:',num2str(alpha),' y:',num2str(y));

            t_newton(y,k)=time(key);
            
            met='cvx'; % method name
            key=strcat('time:: met:',num2str(met),' no equations:',num2str(m),' alpha:',num2str(alpha),' y:',num2str(y));
            t_cvx(y,k)=time(key);
            
            key=strcat('norm:: no equations:',num2str(m),' alpha:',num2str(alpha),' upper bound',' y:',num2str(y));
            x_diff(y,k)=x_values(key);
             
            key1=strcat('funtion:: met:',num2str('newton'),' no equations:',num2str(m),' alpha:',num2str(alpha),' y:',num2str(y));
            key2=strcat('funtion:: met:',num2str('cvx'),' no equations:',num2str(m),' alpha:',num2str(alpha),' y:',num2str(y));

            f_diff(y,k)=abs(x_values(key1)-x_values(key2));
            f_newton(y,k)=(x_values(key1));
            f_cvx(y,k)=(x_values(key2));
           
        end
    end
end

%% error bar plot time

x = 1:1:length(N);
m1 = mean(t_newton);
m2 = mean(t_cvx);

std1=std(t_newton);
std2=std(t_cvx);

%e=errorbar(x,m1,std1,"-s","MarkerSize",5,...
%    "MarkerEdgeColor","red","MarkerFaceColor",[0.65 0.85 0.90])

h1=errorbar(x,m1,std1)
hold on 
h2=errorbar(x,m2,std2)

ylim([0 50])
xlim([0 8])
set(gca,'XTickLabel',{'10','50','100','500','1000','5000','10000'});
%set(gca,'XTick',a,'XTickLabel',a2);
%xtickangle(90)
%set(gca,'FontSize',15);
b = get(gca,'XTickLabel');
set(gca,'XTickLabel',b,'FontName','Times','fontsize',15)
c = get(gca,'YTickLabel');
set(gca,'YTickLabel',c,'FontName','Times','fontsize',15)
box on
legend({'Proposed Method','CVX'}','Orientation','horizontal')
%set(h, {'DisplayName'}, {'Proposed Method','CVX'}')
legend('Orientation','vertical')
%set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('N','FontSize', 15)
ylabel('time(s)','FontSize', 15)
%set(get(gca,'YLabel'),'Rotation',90)

h1.Marker = '*';
h1.MarkerSize = 10;
h1.Color = 'red';
h1.CapSize = 15;

h2.Marker = '+';
h2.MarkerSize = 10;
h2.Color = 'blue';
h2.CapSize = 15;


%% function value comparison

h1=bar(mean(x_diff))
%ylim([0 30])
%xlim([0 6.5])
grid on
set(gca,'XTickLabel',{'10','50','100','500','1000','5000','10000'});
b = get(gca,'XTickLabel');
set(gca,'XTickLabel',b,'FontName','Times','fontsize',15)
c = get(gca,'YTickLabel');
set(gca,'YTickLabel',c,'FontName','Times','fontsize',15)

%legend('Orientation','horizontal')
xlabel('N','FontSize', 15)

%% mean error in x

h1=bar(-mean(f_cvx)+mean(f_newton))
%ylim([0 30])
%xlim([0 6.5])
grid on
set(gca,'XTickLabel',{'10','50','100','500','1000','5000','10000'});
b = get(gca,'XTickLabel');
set(gca,'XTickLabel',b,'FontName','Times','fontsize',15)
c = get(gca,'YTickLabel');
set(gca,'YTickLabel',c,'FontName','Times','fontsize',15)

%legend('Orientation','horizontal')
xlabel('N','FontSize', 15)

%%

function y = fun7(x,m,B,alpha,u,T1,eta)

a=-((x.^(1-alpha))/(1-alpha))-((log((-x+u+eta)))/T1) -((log((x+eta)))/T1);
                    
y=sum(a);
end

      
function y = fun8(x,m,alpha,u,T1,eta)

y=-(x.^(-alpha))+(((-x+u+eta).^-1)/T1)-(((x+eta).^-1)/T1);
           
end

function y = fun9(x,m,alpha,u,T1,eta)

y=alpha*(x.^(-alpha-1))+((-x+u+eta).^-2)/T1+((x+eta).^-2)/T1;
 
end


function y = fun6(b,m,a1,z,c,out)

% inverse 


a= zeros(1,m+1);
i =a1;
    
    for j=1:m+1
         
        if ((i==j) && (i~=m+1))
           
            a(1,j)=c(i,1)/z;
            
        elseif ((i==m+1) && (j==m+1)) 
            a(1,j)=-prod( b , 'all' )/(z*out*out);
        
        elseif ((i~=j) && (j~=m+1) && (i~=m+1))
            
            g=b(setdiff(1:end,[i,j]));
            %g=g(setdiff(1:end,j));
            
            a(1,j)=-prod( g , 'all' )/z;
            
         elseif ((i~=j) && (i==m+1))
            g=b(setdiff(1:end,j));
            a(1,j)=prod( g , 'all' )/(z*out);
            
          elseif ((i~=j) && (j==m+1))
              g=b(setdiff(1:end,i));
            a(1,j)=-prod( g , 'all' )/(z*out);
        end

end

 y=a*out;
end

function [z,c,out]=fun5(z2,m)

out=min(setdiff(z2,min(z2))); %second lowest besides 0
out=1/out;

b=z2*(out);

z=0;
c=zeros(m,1);

for i=1:m
    %for j=1:m
    z1=b(setdiff(1:end,i));
    d = prod( z1 , 'all' );
    z=z+ d;
    e=0;
   
   for j=1:length(z1)
       f=z1(setdiff(1:end,j));
    
        e =e+ prod( f , 'all' );
   end
   c(i,1)=e;
end

end
      
