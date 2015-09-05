clear
clc;
ngp=3;

int_mat=[];
   for i=1:ngp-1;
     for j=i+1:ngp;
      for k=i+1:ngp;
         int_mat=[int_mat; [j k i]];
      end;
     end;
   end;

n=size(int_mat,1);

s=ones(n,1);
numeq=3;
done=0;
while (~done);
[s]=[s backsucc(s(:,end), numeq*ones(n,1))];
done=(sum(s(:,end)-ones(n,1)) == 0)
end
s