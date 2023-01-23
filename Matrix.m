function Zmatrix=Matrix(vector,rw,col);
    Zmatrix=zeros(rw,col);
    I=1;
    for k=1:col;
       for q=1:rw;
          
        Zmatrix(q,k)=vector(I);
       I=I+1;
       end
    end
   
   
end