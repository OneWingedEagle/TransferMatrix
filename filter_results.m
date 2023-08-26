function filter_results
  
T1=0.5; %% astaneh obour 
F1=.5; %% astaneh charkhesh  

[filename,filepath]=uigetfile('*.txt', 'Select results file')

file=strcat(filepath,strcat('\',filename));
data=csvread(file,3,0);


ndata=size(data,1);

Fn(:,1)=data(:,1);

wn1=Fn(1);
wn2=Fn(ndata);

rowsize=size(data,2);

fid = fopen('filtered.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'[order * index   * wn     * Fr      * Tz       Tx ]\n');  

p=0;

for k=1:ndata
  
  Fr=data(k,2);
  Tr=data(k,3);

  if(Tr>T1 && abs(Fr)>F1)
 
  p=p+1;
  result(p,1)= p;
  result(p,2)= k;
  result(p,3)= data(k,1);
  result(p,4)= Fr;
  result(p,5)= Tr;
  Tx=0;
  if(rowsize==4)
  Tx=data(k,4);
  end
  result(p,6)= Tx;
  
  fprintf(fid,'%5d, %5d, %8.4f, %8.4f, %8.4f, %8.4f\n',result(p,1),result(p,2),result(p,3),result(p,4),result(p,5),result(p,6));
  end
 
end

    disp('Number of data found = ');
    p
   fclose(fid);
