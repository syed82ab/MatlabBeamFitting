function res = FitAiry2D_rot( P,data,X,Y )
data=imrotate(data,P(end),'bilinear','crop');
res=data-airy2d(P(1:5),X,Y);
end

function fct=airy2d(P,X,Y)
% gaussienne(P,x,y) 2D gaussian
% P parameter 
% P(1)=offset
% P(2)=amplitude
% P(3)=center along x
% P(4)=center along y
% P(5)=sigma along x
% P(6)=sigma along y


% fct=P(1)+P(2)*exp(-(X-P(3)).^2*2/P(5)^2).*exp(-(Y-P(4)).^2*2/P(6)^2);
R=sqrt((X-P(3)).^2+(Y-P(4)).^2);
R(R==0)=eps;

fct=P(1)+P(2)*(2*besselj(1,R./P(5))./(R./P(5))).^2;

end
