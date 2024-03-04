function WithoutNeg=RemoveNeg(Image)
for i=1:2:size(Image,1)
    for j=1:2:size(Image,2)
        sub=Image(i:i+1,j:j+1);
        Count=length(sub(sub<0));
        if Count ==1 && sub(1,1)<0
            b=sort(sub(:));
            b=b(b>0);
            a=randi([0,floor(b(1))]);
            while -(sub(1,1)-a)/3>sub(1,2) || -(sub(1,1)-a)/3>sub(2,2) || -(sub(1,1)-a)/3>sub(2,1)
                a=a-1;
                if a==0
                    break
                end
            end
            sub(1,2)=sub(1,2)+(sub(1,1)-a)/3;
            sub(2,2)=sub(2,2)+(sub(1,1)-a)/3;
            sub(2,1)=sub(2,1)+(sub(1,1)-a)/3;
            sub(1,1)=a;
        elseif Count==1 && sub(2,2)<0
            b=sort(sub(:));
            b=b(b>0);
            a=randi([0,floor(b(1))]);
            while -(sub(2,2)-a)/3>sub(1,1) || -(sub(2,2)-a)/3>sub(1,2) || -(sub(2,2)-a)/3>sub(2,1)
                a=a-1;
                if a==0
                    break
                end
            end
            sub(1,1)=sub(1,1)+(sub(2,2)-a)/3;
            sub(2,1)=sub(2,1)+(sub(2,2)-a)/3;
            sub(1,2)=sub(1,2)+(sub(2,2)-a)/3;
            sub(2,2)=a;
        elseif Count==1 && sub(1,2)<0
            b=sort(sub(:));
            b=b(b>0);
            a=randi([0,floor(b(1))]);
            while -(sub(1,2)-a)/3>sub(1,1) || -(sub(1,2)-a)/3>sub(2,2) || -(sub(1,2)-a)/3>sub(2,1)
                a=a-1;
                if a==0
                    break
                end
            end
            sub(1,1)=sub(1,1)+(sub(1,2)-a)/3;
            sub(2,1)=sub(2,1)+(sub(1,2)-a)/3;
            sub(2,2)=sub(2,2)+(sub(1,2)-a)/3;
            sub(1,2)=a;
        elseif Count==1 && sub(2,1)<0
            b=sort(sub(:));
            b=b(b>0);
            a=randi([0,floor(b(1))]);
            while -(sub(2,1)-a)/3>sub(1,2) || -(sub(2,1)-a)/3>sub(2,2) || -(sub(2,1)-a)/3>sub(1,1)
                a=a-1;
                if a==0
                    break
                end
            end
            sub(1,1)=sub(1,1)+(sub(2,1)-a)/3;
            sub(1,2)=sub(1,2)+(sub(2,1)-a)/3;
            sub(2,2)=sub(2,2)+(sub(2,1)-a)/3;
            sub(2,1)=a;
        end
        WithoutNeg(i:i+1,j:j+1)=sub;
    end
end
end