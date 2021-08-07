function trans_mat_MO

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

t1=cputime;

dwn=(wn2-wn1)/ndiv;

Fr_hom=zeros(ndiv+1,1);

for p=1:ndiv+1
    
    %p
    
     Fn(p)=wn1+dwn*(p-1);
    
    wvlen=Fn(p);


   [Ts Rs,Fr]=TransferMatrixMultiLayer(eps1,epsdef,eps2,gama1,gamadef,...
   gama1,dd1,ddef,dd2,N1,Ndef,N2,wvlen,p,theta,show_crystal);
           
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

             axis([wn1,wn2,hmin,hmax]);
             hold on
             
        %   plot(Fn,Fr_hom,'+k');
        %     hold on

end
            
if(transmit &&length(Tt)>1)
              figure(2)
               plot(Fn,Tt,colT);
                 axis([wn1,wn2,min(Tt),max(Tt)]);
                %axis([wn1,wn2,0,1]);
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



function [Ts Rs Fr]=TransferMatrixMultiLayer(eps1,epsdef,eps2,...,
  gama1,gamadef,gama2,dd1,ddef,dd2,N1,Ndef,N2,wvlen,p,theta,show_crystal)

variable_d1=[0 0 0];
variable_d2=[0 0 0];
maxd=[0 0 0];
mind=[1.0e10 1.0e10 1.0e10];
for i=1:N1
for j=1:3
w=dd1(i,j);
if(w>maxd(j)) maxd(j)=w; end
if(w<mind(j)) mind(j)=w; end
end
end
for j=1:3
  if(maxd(j)~=mind(j))
  variable_d1(j)=1;
 end
end

maxd=[0 0 0];
mind=[1.0e10 1.0e10 1.0e10];
for i=1:N2
for j=1:3
w=dd2(i,j);
if(w>maxd(j)) maxd(j)=w; end
if(w<mind(j)) mind(j)=w; end
end
end

for j=1:3
  if(maxd(j)~=mind(j))
  variable_d2(j)=1;
 end
end
%%%%%%%%%%%%%%%%%%%%
 
lam=wvlen;

omega=2*pi/lam*cosd(theta);

Ts=1;
 Rs=0;
  Fr=0;
 
  %$$$$============= stacks before defect
  T1=TransferMatrix(eps1(1),gama1(1),omega,dd1(1,1));
  T2=TransferMatrix(eps1(2),gama1(2),omega,dd1(1,2));
  T3=TransferMatrix(eps1(3),gama1(3),omega,dd1(1,3));
  %%%%%%%%%%%%%%%
   T123=T1*T2*T3;
   TT1=T123;
     
  for n=2:N1
   if(variable_d1(1)==1)
    T1=TransferMatrix(eps1(1),gama1(1),omega,dd1(n,1));
   end
   if(variable_d1(2)==1)
    T2=TransferMatrix(eps1(2),gama1(2),omega,dd1(n,2));
   end
   if(variable_d1(3)==1)
    T3=TransferMatrix(eps1(3),gama1(3),omega,dd1(n,3));
   end
   if(variable_d1(1)+variable_d1(2)+variable_d1(3)>0)
    T123=T1*T2*T3;
   end
    TT1=T123*TT1;
  end
  
    %$$$$============= defect layer
  if(Ndef>0)
  T1=TransferMatrix(epsdef(1),gamadef(1),omega,ddef(1));
  T2=TransferMatrix(epsdef(2),gamadef(2),omega,ddef(2));
  T3=TransferMatrix(epsdef(3),gamadef(3),omega,ddef(3));
   T123=T1*T2*T3;
   Tdef=T123;
  for n=2:Ndef
    Tdef=T123*Tdef;
  end
  
end


  %$$$$============= stacks after defect
  if(N2>0)

  T1=TransferMatrix(eps2(1),gama2(1),omega,dd2(1,1));
  T2=TransferMatrix(eps2(2),gama2(2),omega,dd2(1,2));
  T3=TransferMatrix(eps2(3),gama2(3),omega,dd2(1,3));
  %%%%%%%%%%%%%%%
   T123=T1*T2*T3;
   TT2=T123;
  for n=2:N2
   if(variable_d2(1)==1)
      T1=TransferMatrix(eps2(1),gama2(1),omega,dd2(n,1));
   end   
   if(variable_d2(2)==1)   
      T2=TransferMatrix(eps2(2),gama2(2),omega,dd2(n,2));
   end   
   if(variable_d2(3)==1)   
      T3=TransferMatrix(eps2(3),gama2(3),omega,dd2(n,3));
   end
     if(variable_d2(1)+variable_d2(2)+variable_d2(3)>0)
      T123=T1*T2*T3;
   end
    TT2=T123*TT2;
  end
end

  T=TT1;
  if(Ndef>0)
  T=Tdef*T;
  end;
  if(N2>0)
  T=TT2*T;
  end;
  
  
  D0=DynamicMatrix(1,0,omega);
  
  TT=inv(D0)*T*D0;
 
  ang=TT(3,3)/TT(1,1);
  
  Fr=angle(ang)*90./pi();
  
  t2=1./(TT(1,1)*conj(TT(1,1)));

 Ts=abs(t2);

end


function [Dn]=DynamicMatrix(e1, Q,omega)

  pv1=(1./sqrt(2))*transpose([1 1i 0]);
  pv2=(1./sqrt(2))*transpose([-1 -1i 0]);
  pv3=(1./sqrt(2))*transpose([1 -1i 0]);  
  pv4=(1./sqrt(2))*transpose([-1 1i 0]);
  

  kz=omega*[0 0 1]';
  
  kv1=sqrt(e1*(1+Q))*kz;
  kv2=-sqrt(e1*(1+Q))*kz;
  kv3=sqrt(e1*(1-Q))*kz; 
  kv4=-sqrt(e1*(1-Q))*kz; 

  
  qv1=1/omega*cross(kv1,pv1);
  qv2=1/omega*cross(kv2,pv2);    
  qv3=1/omega*cross(kv3,pv3);
  qv4=1/omega*cross(kv4,pv4);
  
  Dn=[pv1(1) pv2(1) pv3(1) pv4(1);...
     qv1(2) qv2(2) qv3(2) qv4(2);...
    pv1(2) pv2(2) pv3(2) pv4(2);...
    qv1(1) qv2(1) qv3(1) qv4(1)];
end


function [T]=TransferMatrix(eps,gamma,omega,d)
  
     Q=gamma/eps;

  Dn=DynamicMatrix(eps,Q,omega);
  
  Pn=PropagationMatrix(eps,Q,omega,d);

   T=Dn*Pn*inv(Dn);
    
end

function [Pn]=PropagationMatrix(e1, Q,omega,d)


  kz=omega*[0 0 1]';
  
  kv1=sqrt(e1*(1+Q))*kz;
  kv2=-sqrt(e1*(1+Q))*kz;
  kv3=sqrt(e1*(1-Q))*kz; 
  kv4=-sqrt(e1*(1-Q))*kz; 
  
   phi1=kv1(3)*d;
   phi2=kv2(3)*d;;
   phi3=kv3(3)*d;;  
   phi4=kv4(3)*d;;
   
   Pn=zeros(4,1);
   
   Pn(1,1)=exp(-1i*phi1);
    Pn(2,2)=exp(-1i*phi2);
     Pn(3,3)=exp(-1i*phi3);
      Pn(4,4)=exp(-1i*phi4);

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
