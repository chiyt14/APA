function [fun_approx, OrderAndCoef, collocationPoint, funEval, FunEvalTimes, BasisCoef, BasisCoef_history, LegendreCoef] = APA(target_fun, D_in, D_out, err_con, max_order, n_eval, c_point, f_eval, symbols)
% ------Implementing the APA method-----
%    Author: Yingtian Chi, Yiwei Qiu
%    Date: Oct, 10, 2020
%    Version:1.0
% ----------------Input-----------------
%   'target_fun' --- the target function to be approximated. Function handle.
%   'D_in' --- the dimension of inputs (number of inputs). Integer.
%   'D_out' --- the dimension of outputs (number of outputs). Integer.
%   'err_con' --- the error control parameter. 1 * D_in array.
%   'max_order' --- the maximum order. 1 * D_in array.
%   'n_eval' --- the number of evaluated sampling points. Non-negative integer.
%   'c_point' --- the evaluated sampling points (if exists). D_in * n_eval array.
%   'f_eval' --- the value of target_fun on 'c_point'. D_out * n_eval array.
%   'symbols' --- symbolic representations of the inputs. 1 * D_in array.

% ----------------Output-----------------
%   'fun_approx' --- approximating results of target_fun. Symbolic expressions.
%   'collocationPoint' --- the evaluated sampling points (including c_point in the input)
%   'funEval' --- function values on the 'collocationPoint'
%   'OrderAndCoef' --- Coefficients of each basis function on each output
%   'BasisCoef', 'BasisCoef_history', 'LegendreCoef' --- Intermediate variables for constructing the approximated results 
% ----------------------------------------



% -----check the input (type, range and number)----
if rem(D_in,1)~=0 | rem(D_out,1)~=0 | rem(max_order,1)~=0 | D_in<=0 | D_out<=0 | max_order<=0
    error('INPUT ERROR: The maximum order and the dimension of inputs and outputs should be positive integers.');
    return;
end
if size(c_point,2)~=size(f_eval,2) | size(c_point,2)~=n_eval
    error('INPUT ERROR: The dimensions of c_point and f_eval do not match with n_eval.');
    return
end
if n_eval~=0 & (size(c_point,1)~=D_in | size(f_eval,1)~=D_out)
    error('INPUT ERROR: The dimensions of c_point and f_eval do not match with D_in and D_out.');
    return
end
if any(err_con<=0) | size(err_con,2)~=D_out
    error('INPUT ERROR: Err_con should be positive and its dimension should match with D_out.');
    return
end

% --------main body of the APA method--------
collocationPoint=c_point;
funEval=f_eval;
FunEvalTimes=0;
D=D_in;% dimension of input
D_O=D_out; % dimension of output
inctol=err_con;
OrderMax=max_order;

LegendreCoef=cell(1,OrderMax+1+1);
ROOT_LEGENDRE=cell(1,OrderMax+1+1);
syms x;
for i=1:OrderMax+1+1
    LegendreCoef{i}=sym2poly(legendreP(i-1,x));
    ROOT_LEGENDRE{i}=roots(LegendreCoef{i});
end
A_pre=cell(1,OrderMax+1);
for i=1:OrderMax+1
    A_pre{i}=legendreP(ones(i).*[0:i-1],ones(i).*ROOT_LEGENDRE{i+1});
end
orderCandidate=cell(1,D_O);
incCandidate=cell(1,D_O);
BasisCoef=cell(1,D_O);% Results during the adaptive process. Different grids are used for different outputs.
FinalBasisCoef=cell(1,D_O); % Results after the adaptive process. The same grid is used for different outputs.
BasisCoef_history=cell(1,D_O);
FinalGrid=cell(1,D_O+1);% 存储每一维输出最后使用的网格。最后一维是网格的并集（不是简单相加）
CollocationPointForBasis=cell(1,(OrderMax+1)^D);% Stor the requried collocation points (sampling points) for each full basis set.
for i=1:D_O
    orderCandidate{i}=zeros(D,1);
    incCandidate{i}=ones(1,1);
    BasisCoef{i}=zeros(1,(OrderMax+1)^D);
    FinalBasisCoef{i}=zeros(1,(OrderMax+1)^D);
    BasisCoef_history{i}=cell(1,(OrderMax+1)^D);
    FinalGrid{i}=zeros(1,(OrderMax+1)^D);
    FinalGrid{1+D_O}=zeros(1,(OrderMax+1)^D);
    basisnorm=ones(1,(OrderMax+1)^D);
    if D>1
        basisnorm=reshape(basisnorm,ones(1,D).*(OrderMax+1));
        BasisCoef{i}=reshape(BasisCoef{i},ones(1,D).*(OrderMax+1));
        FinalBasisCoef{i}=reshape(FinalBasisCoef{i},ones(1,D).*(OrderMax+1));
        BasisCoef_history{i}=reshape(BasisCoef_history{i},ones(1,D).*(OrderMax+1));
        FinalGrid{i}=reshape(FinalGrid{i},ones(1,D).*(OrderMax+1));
        FinalGrid{1+D_O}=reshape(FinalGrid{1+D_O},ones(1,D).*(OrderMax+1));
        CollocationPointForBasis=reshape(CollocationPointForBasis,ones(1,D).*(OrderMax+1));
    end
end
multi_ix=zeros(1,D);
while true
    for i_temp=1:D
        basisnorm(My_sub2ind(size(basisnorm),multi_ix+1))=basisnorm(My_sub2ind(size(basisnorm),multi_ix+1))*(1/(2*multi_ix(i_temp)+1));
    end
    for i_temp=1:D
        multi_ix(i_temp)=multi_ix(i_temp)+1;
        if multi_ix(i_temp)>OrderMax
            multi_ix(i_temp)=0;
        else
            break;
        end
    end
    if all(multi_ix==0)
        break
    end
end
hasOrderCandidate=true;
while hasOrderCandidate
    for i_OUTPUT=1:D_O
        % choosing new order from candidates
        if isempty(orderCandidate{i_OUTPUT})
            continue;
        end
        [~,ix]=max(incCandidate{i_OUTPUT});
        orderNow=orderCandidate{i_OUTPUT}(:,ix(1))';
        disp(['Calculating order: ',num2str(orderNow),' for output #',num2str(i_OUTPUT),', remained candidates: ',num2str(length(incCandidate{i_OUTPUT})-1)]);
        incCandidate{i_OUTPUT}(ix(1))=[];
        orderCandidate{i_OUTPUT}(:,ix(1))=[];
        multi_ix=ones(1,D);
        Col=[];
        temp=zeros(D,1);
        while true % find all the collocation point (sampling points)
            for ii=1:D
                root=ROOT_LEGENDRE{orderNow(ii)+1+1};
                temp(ii)=root(multi_ix(ii));
            end
            Col=[Col,temp];
            % -----index update-----
            for ii=1:D
                multi_ix(ii)=multi_ix(ii)+1;
                if multi_ix(ii)>orderNow(ii)+1
                    multi_ix(ii)=1;
                else
                    break;
                end
            end
            if all(multi_ix==1)
                break
            end
            % -----index update-----
        end
        CollocationPointForBasis{My_sub2ind(size(CollocationPointForBasis),orderNow+1)}=Col;
        Num=size(Col,2);
        Func=zeros(D_O,Num);
        for k=1:Num % evaluate the collocation points (sampling points)
            if isempty(collocationPoint)
                Func(:,k)=target_fun(Col(:,k));
                collocationPoint=[collocationPoint,Col(:,k)];
                funEval=[funEval,Func(:,k)];
                FunEvalTimes=FunEvalTimes+1;
                disp(['function evaluation #',num2str(FunEvalTimes)]);
                disp(['remaining sampling points to be eval in this loop: ',num2str(Num-k)]);
                continue;
            end
            ix=find(collocationPoint(1,:)==Col(1,k)); % find the sampling points that have been evaluated
            for ii=2:D
                ix=intersect(ix,find(collocationPoint(ii,:)==Col(ii,k)));
            end
            if isempty(ix)
                Func(:,k)=target_fun(Col(:,k));
                collocationPoint=[collocationPoint,Col(:,k)];
                funEval=[funEval,Func(:,k)];
                FunEvalTimes=FunEvalTimes+1;
                disp(['function evaluation #',num2str(FunEvalTimes)]);
                disp(['remaining sampling points to be eval in this loop: ',num2str(Num-k)]);
            else
                Func(:,k)=funEval(:,ix(1));
            end
        end
        rootN=size(Col,2);
        N_A=orderNow+1;
        A=1;
        N_inner=1;
        N_outter=prod(N_A);
        for ii=1:D
            N_outter=N_outter/N_A(ii);
            A=A.*kron(ones(N_outter),kron(A_pre{N_A(ii)},ones(N_inner)));
            N_inner=N_inner*N_A(ii);
        end
        
        b=A\Func(i_OUTPUT,:)'; % calculate the coeffcients
        if D>1
            incNow=reshape(b(:,1),orderNow+1);
        else
            incNow=b(:,1);
        end
        BasisCoef_history{i_OUTPUT}{My_sub2ind(size(BasisCoef_history{i_OUTPUT}),orderNow+1)}=incNow;
        
        if D>1
            incNow=zeros(orderNow+1);
        else
            incNow=zeros(1,orderNow+1);
        end
        multi_ix=zeros(1,D);
        ix=find(orderNow);
        FinalGridFlag=true;
        if FinalGrid{1+D_O}(My_sub2ind(size(FinalGrid{i_OUTPUT}),orderNow+1))==1
            FinalGridFlag=false;
        end
        while true
            N_A=orderNow+1-multi_ix;
            ixx=[];
            N_inner=1;
            N_outter=prod(N_A);
            for ii=1:D
                N_outter=N_outter/N_A(ii);
                ixx=[ixx,kron(ones(N_outter,1),kron((1:N_A(ii))',ones(N_inner,1)))];
                N_inner=N_inner*N_A(ii);
            end
            incNow(My_sub2ind(size(incNow),ixx))=incNow(My_sub2ind(size(incNow),ixx))+(-1)^sum(multi_ix)*...
                reshape(BasisCoef_history{i_OUTPUT}{My_sub2ind(size(BasisCoef_history{i_OUTPUT}),orderNow+1-multi_ix)},...
                size(incNow(My_sub2ind(size(incNow),ixx))));
            FinalGrid{i_OUTPUT}(My_sub2ind(size(FinalGrid{i_OUTPUT}),orderNow+1-multi_ix))=...
                FinalGrid{i_OUTPUT}(My_sub2ind(size(FinalGrid{i_OUTPUT}),orderNow+1-multi_ix))+(-1)^sum(multi_ix);
            if FinalGridFlag
                FinalGrid{1+D_O}(My_sub2ind(size(FinalGrid{1+D_O}),orderNow+1-multi_ix))=...
                    FinalGrid{1+D_O}(My_sub2ind(size(FinalGrid{1+D_O}),orderNow+1-multi_ix))+(-1)^sum(multi_ix);
            end
            for ii=1:length(ix)
                multi_ix(ix(ii))=multi_ix(ix(ii))+1;
                if multi_ix(ix(ii))>1
                    multi_ix(ix(ii))=0;
                else
                    break;
                end
            end
            if ~any(multi_ix(ix))
                break;
            end
        end
        
        inc=norm(reshape(incNow,1,[]),2);
        disp(['The 2-norm of coefficients increment after adding ',num2str(orderNow),' is: ',num2str(inc)]);
        N_A=orderNow+1;
        ixx=[];
        N_inner=1;
        N_outter=prod(N_A);
        for ii=1:D
            N_outter=N_outter/N_A(ii);
            ixx=[ixx,kron(ones(N_outter,1),kron((1:N_A(ii))',ones(N_inner,1)))];
            N_inner=N_inner*N_A(ii);
        end
        BasisCoef{i_OUTPUT}(My_sub2ind(size(BasisCoef{i_OUTPUT}),ixx))=BasisCoef{i_OUTPUT}(My_sub2ind(size(BasisCoef{i_OUTPUT}),ixx))+...
            reshape(incNow,size(BasisCoef{i_OUTPUT}(My_sub2ind(size(BasisCoef{i_OUTPUT}),ixx))));
        if inc>inctol(i_OUTPUT) % adding new order candidate
            for ii=1:D
                if orderNow(ii)==OrderMax
                    continue
                end
                orderNew=orderNow;
                orderNew(ii)=orderNew(ii)+1;
                flag=true; % judge if the new point is admissible
                for jj=1:D
                    temp=orderNew;
                    if temp(jj)>0
                        temp(jj)=temp(jj)-1;
                        if isempty(BasisCoef_history{i_OUTPUT}{My_sub2ind(size(BasisCoef_history{i_OUTPUT}),temp+1)})
                            flag=false;
                            break;
                        end
                    end
                end
                if flag
                    orderCandidate{i_OUTPUT}=[orderCandidate{i_OUTPUT},orderNew'];% add orderNew into candidates
                    incCandidate{i_OUTPUT}=[incCandidate{i_OUTPUT},inc];
                    disp(['adding ',num2str(orderNew),' into candidates']);
                end
            end
        end
    end
    
    hasOrderCandidate=false;
    for i_OUTPUT=1:D_O
        if ~isempty(orderCandidate{i_OUTPUT})
            hasOrderCandidate=true;
            break;
        end
    end
end
% find the final collocation point that is needed
FinalCollocationPoint=cell(1,D_O+2);
FinalCollocationPoint{D_O+2}=zeros(1,D); % D_O+2 is the union of 1~D_O
% D_O+1 is the actual final collocation points when all outputs share the same grid, which generally is not equal to D_O+2
% we can compare the number of D_O+1 and D_O+2 to decide to use either the same
% grid or different grids
for i=1:D_O+1
    ix=find(FinalGrid{i});
    FinalCollocationPoint{i}=zeros(1,D);
    for ii=1:length(ix)
        FinalCollocationPoint{i}=union(FinalCollocationPoint{i},CollocationPointForBasis{ix(ii)}','rows');
    end
    FinalCollocationPoint{D_O+2}=union(FinalCollocationPoint{i},FinalCollocationPoint{D_O+2},'rows');
    FinalCollocationPoint{i}=FinalCollocationPoint{i}';
end
FinalCollocationPoint{D_O+2}=FinalCollocationPoint{D_O+2}';

% --------processing & output--------
fun_approx=cell(1,D_O);
OrderAndCoef=cell(1,D_O);
for i_output=1:D_O
    if D_in>1
        multi_ix=zeros(1,ndims(BasisCoef{i_output}));
    else
        multi_ix=0;
    end
    OrderAndCoef{i_output}=[];
    VarName=symbols;
    ExpressionResult=0;
    while true
        ix=My_sub2ind(size(BasisCoef{i_output}),multi_ix+1);
        OrderAndCoef{i_output}=[OrderAndCoef{i_output},[multi_ix,BasisCoef{i_output}(ix)]'];
        temp=BasisCoef{i_output}(ix);
        for ii=1:D
            temp=temp*legendreP(multi_ix(ii),VarName(ii));
        end
        ExpressionResult=ExpressionResult+temp;
        for ii=1:D
            multi_ix(ii)=multi_ix(ii)+1;
            try
                ix=My_sub2ind(size(BasisCoef{i_output}),multi_ix+1);
            catch
                multi_ix(ii)=0;
                continue;
            end
            if isempty(BasisCoef_history{i_output}{ix})
                multi_ix(ii)=0;
            else
                break;
            end
        end
        if all(multi_ix==0)
            break;
        end
    end
    fun_approx{i_output}=vpa(simplify(ExpressionResult),5);
end

end

