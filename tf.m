function tf


clear all

colors = {'-or', '-ob', '-oc', '-ok','-*r', '-*b', '-*c', '-*k','-xr', '-xb',... 
'-r', '-b', '-c', '-k','-g'};

[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
 cd(filepath1)
 fid= fopen(filename1)
%%  ==== Indicent angle
line=getNewDataLine(fid);

numbs = str2num(line);
theta= numbs(1);
theta2=0;
ndivth=0;
if(length(numbs)>2)
theta2=numbs(2);
ndivth=numbs(3);
end

%%=========  stack befor defect
amp1=0;
omgfact1=0;
dd10=[0 0 0];
eps1=[1 1 1];
gama1=[0 0 0];
line=getNewDataLine(fid);
numbs = str2num(line);
dd10(1)=numbs(1);

if(length(numbs)>1)
dd10(2)=numbs(2);
end
if(length(numbs)>2)
dd10(3)=numbs(3);
end

if(length(numbs)>3)
amp1=numbs(4);
end

if(length(numbs)>4)
omgfact1=numbs(5);
end

line=getNewDataLine(fid);
numbs = str2num(line);
eps1(1)=numbs(1);

if(length(numbs)>1)
eps1(2)=numbs(2);
end
if(length(numbs)>2)
eps1(3)=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
gama1(1)=numbs(1);

if(length(numbs)>1)
gama1(2)=numbs(2);
end
if(length(numbs)>2)
gama1(3)=numbs(3);
end
%%=============

%%=========   defect stack
ddef=[0 0 0];
epsdef=[1 1 1];
gamadef=[0 0 0];
line=getNewDataLine(fid);
numbs = str2num(line);
ddef(1)=numbs(1);

if(length(numbs)>1)
ddef(2)=numbs(2);
end
if(length(numbs)>2)
ddef(3)=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
epsdef(1)=numbs(1);

if(length(numbs)>1)
epsdef(2)=numbs(2);
end
if(length(numbs)>2)
epsdef(3)=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
gamadef(1)=numbs(1);

if(length(numbs)>1)
gamadef(2)=numbs(2);
end
if(length(numbs)>2)
gamadef(3)=numbs(3);
end
%%=============

%%=========  stack befor defect
dd20=[0 0 0];
eps2=[1 1 1];
gama2=[0 0 0];
amp2=0;
omgfact2=0;
line=getNewDataLine(fid);
numbs = str2num(line);
dd20(1)=numbs(1);

if(length(numbs)>1)
dd20(2)=numbs(2);
end
if(length(numbs)>2)
dd20(3)=numbs(3);
end

if(length(numbs)>3)
amp2=numbs(4);
end

if(length(numbs)>4)
omgfact2=numbs(5);
end

line=getNewDataLine(fid);
numbs = str2num(line);
eps2(1)=numbs(1);

if(length(numbs)>1)
eps2(2)=numbs(2);
end
if(length(numbs)>2)
eps2(3)=numbs(3);
end
line=getNewDataLine(fid);
numbs = str2num(line);
gama2(1)=numbs(1);

if(length(numbs)>1)
gama2(2)=numbs(2);
end
if(length(numbs)>2)
gama2(3)=numbs(3);
end
%%==========================


line=getNewDataLine(fid);
numbs = str2num(line);
N1=numbs(1);
Ndef=0;
N2=0;
if(length(numbs)>1)
Ndef=numbs(2);
end
if(length(numbs)>2)
N2=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
transmit=numbs(1);
rotation=numbs(2);
show_crystal=0;
if(length(numbs)>2)
show_crystal=numbs(3);
end

if(transmit>0)
colT=colors{transmit};
else
colT=colors{1};
end

if(rotation>0)
colR=colors{rotation+1};
else
colR=colors{2};
end

line=getNewDataLine(fid);
numbs = str2num(line);

wn1=numbs(1);
wn2=numbs(2);

global ndiv;
ndiv=numbs(3);

global epsin;
global epsout;
epsin=1;
epsout=1;
%line=getNewDataLine(fid);
%numbs = str2num(line);
%if(length(numbs)>0)
%epsin=numbs(1);
%epsout=numbs(2);
%end

wn0=wn1;

if(ndivth>0)
wn1=theta;
wn2=theta2;
ndiv=ndivth;
end

%====================================================
%== Variable thinkness begins
% the following lines till [Variable thinkness ends]..
%define variable thinkness. 
%If amp1=amp2=0. thinkness dont vary from layer to layer
%amp1=.25; % bteween 0 and 1.0. 
%amp2=.25; % bteween 0 and 1.0.
omg1=0;
omg2=0;
if(N1>1)
omg1=omgfact1*pi;%%(N1-1);% you can try another number instead of 1.0
end
if(N2>1)
omg2=omgfact2*pi;%%/(N2-1);% you can try another number instead of 1.0
end
dd1=zeros(N1,3);
dd2=zeros(N2,3);

tk1=dd10(1)+dd10(2)+dd10(3);

L1=N1*tk1;
for i=1:N1
dist=tk1*(i-.5);
sinfact=sin(dist/L1*omg1);
for j=1:3
if(j==1)
dd1(i,j)=dd10(j)+tk1*amp1*sinfact;
elseif(j==2)
dd1(i,j)=dd10(j)-tk1*amp1*sinfact;
else
dd1(i,j)=dd10(j);
end
end
end

%figure(11)
% plot(dd1(:,1),colR);

tk2=dd20(1)+dd20(2)+dd20(3);
L2=N2*tk2;

for i=1:N2
dist=L2-(i-0.5)*tk2;
sinfact=sin(dist/L2*omg2);
for j=1:3
if(j==3)
dd2(i,j)=dd20(j)+tk2*amp2*sinfact;
elseif(j==2)
dd2(i,j)=dd20(j)-tk2*amp2*sinfact;
else 
dd2(i,j)=dd20(j);
end
end
end
%== Variable thinkness ends
%====================================================

TE=0;

t1=cputime;

dwn=(wn2-wn1)/ndiv;

Fr_hom=zeros(ndiv+1,1);
theta0=theta;
for p=1:ndiv+1
    
    %p
    
     Fn(p)=wn1+dwn*(p-1);
     
    if(ndivth>0)
    wvlen=wn0;
    theta=Fn(p);
    else
    wvlen=Fn(p); 
    end

   [Ts Rs,Fr]=calculteFaraday(eps1,epsdef,eps2,dd1,ddef,dd2,N1,Ndef,N2,wvlen,...
   p,theta,TE);
           
     if(real(Ts)>1) 
      %  Ts=1;
      end
      
    Tr(p)=real(Fr);
    Tt(p)=real(Ts);
    
    uu=Ts+Rs;

end


t2=cputime;

comptation_time=t2-t1;
comptation_time;
        
        Tr';
        Tt';
	
result=zeros(ndiv+1,3);


%fidx = fopen('analyt.txt','wt');  % Note the 'wt' for writing in text mode
fid = fopen('results.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'[nGx *  nGy]\n');  
  
fprintf(fid,'**************b\n');  
fprintf(fid,'[wn *  Rotation * Transmitance ]\n');  


  for p=1:ndiv+1
  
  result(p,1)= Fn(p);
  result(p,2)= Tr(p);
  result(p,3)= Tt(p);

  fprintf(fid,'%f, %f, %f\n',result(p,1),result(p,2),result(p,3));
 % fprintf(fidx,'%f\t%f\n',result(p,1),Fr_hom(p));


  end
    fclose(fid);
    % fclose(fidx);
     
  disp('Results:');
  disp('[wn *  Rotation * Transmitance ]');	
  disp(result);
  
  rotation=0;
if(rotation &&length(Tr)>1)
                figure(1)
             plot(Fn,Tr,colR);

             hmax=50;
             
             ttmax=max(Tr);   

            if(ttmax<=0)
             hmax=0; 
            end
            if(ttmax>45) 
             div= round(ttmax/5);
            
             hmax=(div+1)*5;
             
           end
           
              ttmin=min(Tr);
           if(ttmin>=0) hmin=0;
          else 
           hmin=ttmin*1.1;
           end
           if(hmin<-45) 
             div= round(hmin/5);
            
             hmin=-(div+1)*5;
             
             end
hmin=0;
             axis([wn1,wn2,hmin,hmax]);
             hold on
             
        %   plot(Fn,Fr_hom,'+k');
        %     hold on

end
            
if(transmit &&length(Tt)>1)
              figure(2)
               plot(Fn,Tt,colT);
                 axis([wn1,wn2,min(Tt),max(Tt)]);
                axis([wn1,wn2,0,1]);
                 hold on
            
end




if(show_crystal==1)
figure(10)
x0=0;
y0=0;
w0=1;
h=1;
c1=['r','y','c'];

ww=0;
for i=1:N1
for j=1:3
w=dd1(i,j);
if(w>0)
rectangle('Position',[x0+ww y0 w h],'FaceColor',c1(j),'EdgeColor','None');
ww=ww+w;
end
end
end

c2=['g','g','g'];
for i=1:Ndef
for j=1:3
w=ddef(j);
if(w>0)
rectangle('Position',[x0+ww y0 w h],'FaceColor',c2(j),'EdgeColor','None');
ww=ww+w;
end
end
end

c3=['c','y','r'];

for i=1:N2
for j=1:3
w=dd2(i,j);
if(w>0)
rectangle('Position',[x0+ww y0 w h],'FaceColor',c3(j),'EdgeColor','None');
ww=ww+w;
end
end
end

rectangle('Position',[x0 y0 ww h]);

axis([0 ww 0 2]);
end

end



function [Ts Rs Fr]=calculteFaraday(eps1,epsdef,eps2,dd1,ddef,dd2,...
N1,Ndef,N2,wvlen,p,theta, TE)

global ndiv;


 lam=wvlen;%*1e-6;
ev=1.60218e-19;
wp=8.24;
wc=0.048;

lam1=wvlen*1e-6;
c0=299792458.0;
w=2.0*pi()*c0/lam1;
em=1-wp^2/(w*(w+1i*wc));
em;
eps1(1,1)=real(em);
eps2(1,3)=real(em);

A = 4.92719645;
B = 7.27691471;
C =11.5786091e2 ;
D = 42.7173925;
E =100;

esi2=A+B*lam1^2/(lam1^2-C)+D*lam1^2/(lam1^2-E);

eps1(1,3)=esi2;
eps2(1,1)=esi2;

 Ts=0;
 Rs=0;
 Fr=0;%%90*sin(2*pi*p/ndiv);
 
 e1=eps1(1);
 e2=eps1(2);
 e3=eps1(3);
  
 n0=1;

 
 n1=sqrt(e1);
 n2=sqrt(e2);
 n3=sqrt(e3);

 p0=cosd(theta)/n0;

d1=dd1(1,1);
d2=dd1(1,2);
d3=dd1(1,3);


 M1=zeros(2,2);
 M2=zeros(2,2);
 M3=zeros(2,2);
 %%%%%%%%%%% layer 1
kz1=2*pi*n1*cosd(theta)/lam;
  phi1=kz1*d1;
  
   P1=zeros(2,2);
   P1(1,1)=exp(1i*phi1);
   P1(2,2)=exp(-1i*phi1);
  
  D1=zeros(2,2);
   D1(1,1)=1;
     D1(1,2)=1;
   D1(2,1)=-kz1;
    D1(2,2)=kz1;
    
    M1=D1*P1*inv(D1);


%%%%%%%%%%% layer 2
  kz2=2*pi*n2*cosd(theta)/lam;
  phi2=kz2*d2;
   P2=zeros(2,2);
   P2(1,1)=exp(1i*phi2);
   P2(2,2)=exp(-1i*phi2);
  
  D2=zeros(2,2);
   D2(1,1)=1;
     D2(1,2)=1;
   D2(2,1)=-kz2;
    D2(2,2)=kz2;
    
  M2=D2*P2*inv(D2);
  
  %%%%%%%%%%% layer 3
  kz3=2*pi*n3*cosd(theta)/lam;
  phi3=kz3*d3;
   P3=zeros(2,2);
   P3(1,1)=exp(1i*phi3);
   P3(2,2)=exp(-1i*phi3);
  
  D3=zeros(2,2);
   D3(1,1)=1;
     D3(1,2)=1;
   D3(2,1)=-kz3;
    D3(2,2)=kz3;
    
  M3=D3*P3*inv(D3);
%%%%%%%%%%%%%%%%%%%%%%%%%%

  ML=M1*M2*M3;
   MU=M3*M2*M1;
   
   MML=eye(2);

  for n=1:N1
    MML=MML*ML;
  end
  
    MMU=eye(2);
    
  for n=1:N2
    MMU=MMU*MU;
  end

  
  %MM=MML*MMU;
  Mdef=eye(2);
 if(Ndef>0)
  nf=sqrt(epsdef(1));
  df=ddef(1,1);
   %%%%%%%%%%% defect layer
   kzf=2*pi*nf*cosd(theta)/lam;
  phif=kzf*df;

   Pf=zeros(2,2);
   Pf(1,1)=exp(1i*phif);
   Pf(2,2)=exp(-1i*phif);
  
  Df=zeros(2,2);
   Df(1,1)=1;
     D3(1,2)=1;
   Df(2,1)=-kzf;
    Df(2,2)=kzf;
    
    Mdef=Df*Pf*inv(Df);
 end
 
  MM=MML*Mdef*MMU;
   %%%%%%%%%%%

  kz0=2*pi*n0*cosd(theta)/lam;  
    D0=zeros(2,2);
    
   D0(1,1)=1;
   D0(1,2)=1;
   D0(2,1)=-kz0;
   D0(2,2)=kz0;

%   M=D0*MM*inv(D0);
    TT=inv(D0)*MM*D0;

    m11=TT(1,1);
  
  t2=1./(m11*conj(m11));

 Ts=abs(t2);
 % Ts=abs(t2);

end


function line=getNewDataLine(fid)

TF=1;
k=0;
while (k<100 && TF==1)
line=fgets(fid);
TF = strncmpi(line,'/',1);
k=k+1;

end


end

