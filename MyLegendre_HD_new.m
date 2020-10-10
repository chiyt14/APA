function z = MyLegendre_HD_new(x,B,B_his,LegendreCoef)
%    Author: Yingtian Chi, Yiwei Qiu
%    Date: Oct, 10, 2020
%    Version:1.0
B_size=size(B);
if length(B_size)==2 && B_size(1)==1
    D=1;
else
    D=ndims(B);
end

z=0;
multi_ix=zeros(1,D);

while true
    ix=My_sub2ind(size(B),multi_ix+1);
    coef=B(ix);
    temp=coef;
    for ii=1:D
        temp=temp.*polyval(LegendreCoef{multi_ix(ii)+1},x(ii,:));
    end
    z=z+temp;
    for ii=1:D
        multi_ix(ii)=multi_ix(ii)+1;
        try
            ix=My_sub2ind(size(B),multi_ix+1);
        catch
            multi_ix(ii)=0;
            continue;
        end
        if isempty(B_his{ix})
            multi_ix(ii)=0;
        else
            break;
        end
    end
    if all(multi_ix==0)
        break;
    end
end

end

