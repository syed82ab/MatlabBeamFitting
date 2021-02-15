function res = FitGauss2D_rot( P,data,X,Y )
data=imrotate(data,P(end),'bilinear','crop');
res=data-gauss2d(P(1:6),X,Y);
end

function fct=gauss2d(P,X,Y)
% gaussienne(P,x,y) 2D gaussian
% P parameter 
% P(1)=offset
% P(2)=amplitude
% P(3)=center along x
% P(4)=center along y
% P(5)=sigma along x
% P(6)=sigma along y


fct=P(1)+P(2)*exp(-(X-P(3)).^2*2/P(5)^2).*exp(-(Y-P(4)).^2*2/P(6)^2);

end
