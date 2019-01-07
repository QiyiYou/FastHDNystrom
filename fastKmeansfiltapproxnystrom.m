%% Optimized version

function img=fastKmeansfiltapproxnystrom(img,S,h,Centre,spatialkernel,convmethod,fast_flag,imgguide)
if ~exist('imgguide','var')
     % guideimage is same as inputimage
     imgguide=img;
end
if ~exist('fast_flag','var')
     % guideimage is same as inputimage
     fast_flag=1;
end
[m,n,~]=size(img);
num=zeros(size(img));
guided=size(imgguide,3);
Cluster=size(Centre,1);
A=zeros(Cluster,Cluster);
for i=1:Cluster
    A(i,i)=1;
    for j=i+1:Cluster
        A(i,j)=exp(-sum((Centre(i,:)-Centre(j,:)).^2,2)/(2*h*h));
        A(j,i)=A(i,j);
    end
end
[Vhat,Dhat]=eig(A);
Vhat = flip(Vhat, 2);
Dhat = flip(flip(Dhat,2),1);
Dhat=diag(Dhat);
W=zeros(m,n,Cluster);
for i=1:Cluster
    W(:,:,i)=sum((imgguide-reshape(Centre(i,:),1,1,guided)).^2,3);   
end
W=exp(-bsxfun(@rdivide,W,(2*(h^2))));
%% Interpolating numerator and denominator seperately
den=zeros(m,n);

if strcmp(convmethod,'matlab')
% Uncomment for matlab implementation
    if strcmp(spatialkernel,'box')
        filt     = ones(2*S+1,2*S+1);       
    elseif strcmp(spatialkernel,'gaussian')       
        w  = round(6*S); if (mod(w,2) == 0); w  = w+1; end
        filt     = fspecial('gaussian', [w w], S);
    else
    end        
%% Calculating Bilateral image for each cluster centre as index pixel
    for i=1:Cluster
        Eigmat=sum(W.*reshape(Vhat(:,i),1,1,Cluster),3);     
        den=den+bsxfun(@times,(1/Dhat(i,1))*Eigmat,imfilter(Eigmat,filt));
        num=num+bsxfun(@times,(1/Dhat(i,1))*Eigmat,imfilter(bsxfun(@times,Eigmat,img),filt));    
    end
end    

if strcmp(convmethod,'O1')
    for i=1:Cluster
        Eigmat=sum(W.*reshape(Vhat(:,i),1,1,Cluster),3);     
        if strcmp(spatialkernel,'box')
            den=den+bsxfun(@times,(1/Dhat(i,1))*Eigmat,box_filter(Eigmat,S,fast_flag));
            num=num+bsxfun(@times,(1/Dhat(i,1))*Eigmat,box_filter(bsxfun(@times,Eigmat,img),S,fast_flag));                      
        elseif strcmp(spatialkernel,'gaussian')
            den=den+bsxfun(@times,(1/Dhat(i,1))*Eigmat,gauss_filter(Eigmat,S,fast_flag));
            num=num+bsxfun(@times,(1/Dhat(i,1))*Eigmat,gauss_filter(bsxfun(@times,Eigmat,img),S,fast_flag));                      
        else
        end 
    end
end 
img=bsxfun(@rdivide,num,den);
end

