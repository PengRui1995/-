function  Pos =Chan3(BSN,BS,R)
%�����Ĺ��ܣ�����Chan���㷨Nվ��λ�ĺ���,���ǽ�վ�����¶�Դλ������һ���Ĺ���
%������������������̨վ�ĺ������꣬�Ѿ�����������վ���ź�Դվ�ľ����Ϳ��Թ��Ƴ��ź�Դվ��λ��
%������ʹ�ã�Pos=Chan3(BSN,BS,R)
%���룺
%     BSN:���ڶ�λ��̨վ��
%     BS:����̨վ�ĺ�������,ÿ�д���һ������վ����
%     R:����������վ���ź�Դվ�ľ����
%�����
%     Pos:���Ƴ��ź�Դվ��λ�õĺ�������,������������꣬��Ҫ��������֪ʶ����ȡ�ᡣ
%���ӣ���;
%ע��������ú��������÷�Χ,����֪�ź�Դ�����վ��Զ��������ʹ��,��ʵ����㷢�ֽ�վ����������á�
%����˵��:��ʽ�г�����h��h'��G,G'��z��z����ʾ��ͬ�������ƣ��������ha,hb,Ga,Gb,za,zb��Ϊ��ʶ������������ת�þ��������

%k  = xi^2+yi^2
k =zeros(BSN,1);
for i =1:BSN
    k(i) = BS(i,1)^2+BS(i,2)^2;
end
% ha =0.5*(ri1^2+k1) i=2,..m
ha = zeros(BSN-1,1);
for i =1:BSN-1
    ha(i) = 0.5*(R(i)^2-k(i+1)+k(1));
end
%Ga = -[xi1,yi1,ri1]
Ga = zeros(BSN-1,3);
for i=1:BSN-1
    Ga(i,1:2)=BS(i+1,:)-BS(1,:);
    Ga(i,3) = R(i);
end
Ga =-Ga;
%Q ��������Э�������
%Q =cov(R);
Q = 0.5*(eye(BSN-1)+ones(BSN-1,BSN-1));
%% ��һ��WLSԶվ�����
za = (Ga'/Q*Ga)\Ga'/Q*ha;

%% ��һ��WLS��վ�����
%Ba =diag��r2,r3..rm��ri=ri1+r1;
r =zeros(BSN-1,1);
r = R+za(3)*ones(BSN-1,1);
Ba =diag(r');
%
Sigma = Ba*Q*Ba;
W = inv(Sigma);
% za ���¹���
za =(Ga'*W*Ga)\Ga'*W*ha;

%% �ڶ���WLS��վ�����
%hb 
hb =zeros(3,1);
hb(1) =(za(1)-BS(1,1))^2;
hb(2) =(za(2)-BS(1,2))^2;
hb(3) =za(3)^2;
%Gb 
Gb=[1,0;0,1;1,1];
%Bb
Bb =diag(za-[BS(1,:),0]');
%zCov =(Ga'*W*Ga)^-1
zCov =inv(Ga'*W*Ga);
%SigmaB  =4*Bb*zCov*Bb
SigmaB = 4*Bb*zCov*Bb;
%zb
zb=(Gb'/SigmaB*Gb)\Gb'/SigmaB*hb;
%%
res = sqrt(zb);
%% Pos
Pos(1,:) = [res(1) res(2)]+BS(1,:);
Pos(2,:) = [res(1) -res(2)]+BS(1,:);
Pos(3,:) = [-res(1) res(2)]+BS(1,:);
Pos(4,:) = [-res(1) -res(2)]+BS(1,:);

%%  ѡ�����һ��wls�Ľ������ĵ�
dis= zeros(4,1);
for i=1:4
    dis(i) = [Pos(i,:)-za(1:2,1)']* [Pos(i,:)-za(1:2,1)']';
end
[val,ind] = min(dis);
%%
Pos = Pos(ind,:);

end