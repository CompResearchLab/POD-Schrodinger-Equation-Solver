function A=normalize(Z,x,y);
      SqrZ=Z.^2;
        A=sqrt(1/trapz(y,trapz(x,SqrZ,2)));
end

