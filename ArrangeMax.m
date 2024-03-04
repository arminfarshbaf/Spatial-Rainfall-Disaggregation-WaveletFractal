%% In order to arrange pixels better
function Arrangement=ArrangeMax(Image)
step=input('Step :');
for i=1:step:size(Image,1)
    for j=1:step:size(Image,2)
        sub=Image(i:i+step-1,j:j+step-1);
        b=flip(sort(sub(:)));
        [~,idx]=max(sub(:));    %vaghti be vector tabdil mikonim sutun ha aval mian na row ha
        [MainRow,MainCol] = ind2sub(size(sub),idx);
        if MainRow==1 && MainCol==1
            r=randsample(3,3);
            for k=2:4
                [Row,Column]=find(sub==b(k));
                c=[0,1,1];
                d=[1,0,1];
                sub(Row,Column)=sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)));
                sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)))=b(k);
            end
        elseif MainRow==step && MainCol==1
            r=randsample(3,3);
            for k=2:4
                [Row,Column]=find(sub==b(k));
                c=[-1,-1,0];
                d=[0,1,1];
                sub(Row,Column)=sub(MainRow+c(k-1),MainCol+d(k-1));
                sub(MainRow+c(k-1),MainCol+d(k-1))=b(k);
            end
        elseif MainRow==1 && MainCol==step
            r=randsample(3,3);
            for k=2:4
                [Row,Column]=find(sub==b(k));
                c=[0,1,1];
                d=[-1,-1,0];
                sub(Row,Column)=sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)));
                sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)))=b(k);
            end
        elseif MainRow==step && MainCol==step
            r=randsample(3,3);
            for k=2:4
                [Row,Column]=find(sub==b(k));
                c=[-1,-1,0];
                d=[0,-1,-1];
                sub(Row,Column)=sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)));
                sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)))=b(k);
             end
        elseif MainCol==1 && MainRow<step && MainRow>1
            r=randsample(5,5);
            for k=2:6
                [Row,Column]=find(sub==b(k));
                c=[-1,-1,0,1,1];
                d=[0,1,1,1,0];
                sub(Row,Column)=sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)));
                sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)))=b(k);
            end
        elseif MainRow==1 && MainCol<step && MainCol>1
            r=randsample(5,5);
            for k=2:6
                [Row,Column]=find(sub==b(k));
                c=[0,1,1,1,0];
                d=[-1,-1,0,1,1];
                sub(Row,Column)=sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)));
                sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)))=b(k);
            end
        elseif MainRow==step && MainCol<step && MainCol>1
            r=randsample(5,5);
            for k=2:6
                [Row,Column]=find(sub==b(k));
                c=[0,-1,-1,-1,0];
                d=[-1,-1,0,1,1];
                sub(Row,Column)=sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)));
                sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)))=b(k);
            end
        elseif MainCol==step && MainRow<step && MainRow>1
            r=randsample(5,5);
            for k=2:6
                [Row,Column]=find(sub==b(k));
                c=[-1,-1,0,1,1];
                d=[0,-1,-1,-1,0];
                sub(Row,Column)=sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)));
                sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)))=b(k);
            end
        else
            r=randsample(8,8);
            for k=2:9
                [Row,Column]=find(sub==b(k));
                c=[-1,-1,-1,0,1,1,1,0];
                d=[-1,0,1,1,1,0,-1,-1];
                sub(Row,Column)=sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)));
                sub(MainRow+c(r(k-1)),MainCol+d(r(k-1)))=b(k);
            end
        end
        Arrangement(i:i+step-1,j:j+step-1)=sub;
    end
end
end