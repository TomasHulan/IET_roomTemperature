function [f,X,T] = localMin(p,k)
%vyberie pri kmitan� body, kde m� funkcia
%lok�lne minimum
%k - po�et okolit�ch bodov na porovn�vanie
%za lok�lne minimum berie bod pre ktor� plat�,
%�e v�etk�ch k*2 okolit�ch bodov je men��ch ako 
%tento bod
%p - vektor, harmonick� priebeh
%f - hodnoty vystupnych bodov, X - indexy lokalnych maxim
%T - perioda maxim
n = length(p);
j = 1;
pR = p(k+2:2*k+1);%body napravo
pL = p(1:k);%body na�avo
m = 0;%index posledn�ho maxima
T = [];
for i = k+1:n-k
    if all(pL > p(i)) && all(pR > p(i))
        f(j) = p(i);  
        pL = p(i-k+1:i);
        if i<n-k
            pR = p(i+2:i+k+1);
        end
        T(j) = i-m;
        X(j) = i;
        j = j+1;
        m=i;
    else
        if i<n-k
            pL = p(i-k+1:i);
            pR = p(i+2:i+k+1);
        end
    end
end
T = round(mean(T));%priemern� rozostup max�m
end