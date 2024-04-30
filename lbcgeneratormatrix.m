clc;
clear all;
close all;
g=[1 1 0 1 0 0 0;0 1 1 0 1 0 0;1 1 1 0 0 1 0;1 0 1 0 0 0 1];
disp(g);
disp('The order of linear block code for given generator matrix is ');
[n,k]=size(transpose(g));
disp('The cord word length is');
disp(n);
disp('The size of message bits are');
disp(k);
for i=1:2^k
    for j=k:-1:1
        disp(j);
        if rem(i-1,2^(-j+k+1))>=2^(-j+k)
            m(i,j)=1;
        else
            m(i,j)=0;
        end
    end
end
disp('The possible message bits are');
 
disp('c0 c1 c2 c3');
disp (m);
disp('The possible codewords are ');
disp('b0 b1 b2 c0 c1 c2 c3 Hamming weight');%%
c=rem(m*g,2);
d_min=sum((c(1:2^k,:))');
d_min2=d_min';
s=[c d_min2];
disp(s);
disp('The minimum hamming weightfor the given block code is=');
d_min1=min(sum((c(2:2^k,:))'));
disp(d_min1);
G_prime = zeros(size(g));

% Perform row operations to get systematic form
for i = 1:size(g, 1)
    % Find the index of the leading 1 in the current row
    leading_1_index = find(g(i,:) == 1, 1);
    
    % Update G_prime with the current row
    G_prime(i, :) = g(i, :);
    
    % Perform row operations to get zeros below the leading 1
    for j = 1:size(g, 1)
        if j ~= i
            factor = g(j, leading_1_index);
            g(j, :) = mod(g(j, :) - factor * g(i, :), 2);
        end
    end
end

% Display the systematic form of the generator matrix G'
disp('Generator Matrix G'' for Systematic Code C''');
disp(G_prime);
%%DECODING PART
n=input('Enter the Size of codeword N=');
k=input('Enter the K='); 
 
P=zeros(k,(n-k));
 
P=input('Enter the Parity Matrix having size of p=');
d=input('Enter the message ')
I=eye(k);
% generator matrix
G=[I P];
C=d*G;
for i=1:n
 if (rem(C(i),2)==0)
 C(i)=0;
 else
 C(i)=1;
 end
end
%starting for the decoder
p=P';
 
I=eye(n-k);
H=[p I];
 
H1=H';
disp('The H matrix is');
disp(H1);
R=input('Enter the Received Code Word');
disp(R);
 
 
%for syndrome
 S=rem(R*H1,2);
disp('Syndrome of a Given codeword is');
disp(S);
 
for i=1:size(H1)
 if(H1(i,1:3)==S)
 R(i)=1-R(i);
 break;
 end
end
disp('The Error is in bit:');
disp(i);
disp('correct receive code is given by');
disp(R);