function [P, data_fit, horizontal_cut, vertical_cut] = fit_waist(Data,resize_factor,Pixel_conv,rot_angle,fixed_rot)
%%fit_waist 
conv=ones(Pixel_conv)/(Pixel_conv^2);
data=conv2(Data,conv,'same');

[data]=imresize(data,1/resize_factor);

[aa,a]=max(data);
[xymax ,x0] =max(aa);
y0=a(x0);
[height, width]=size(data);
x=double(1:width);
y=double(1:height);
[X,Y] = meshgrid(x, y);

min_a=double(min(Data(:))); %background is most abundant pixel value
data= data - min_a;
aa=mean(data);
a=max(aa);
xFWHM(1)=find(aa>a/2,1,'first');
xFWHM(2)=find(aa>a/2,1,'last');
Sigmax=diff(xFWHM)/sqrt(2*log(2));
aa=mean(data,2);
a=max(aa);
yFWHM(1)=find(aa>a/2,1,'first');
yFWHM(2)=find(aa>a/2,1,'last');
Sigmay=diff(yFWHM)/sqrt(2*log(2));

P0=[min_a,xymax,x0,y0,Sigmax,Sigmay,rot_angle];

options=optimset( 'Display','final-detailed');
if fixed_rot
    fit_fn= @(param) FitGauss2D_rot([param rot_angle],data,X,Y);
    P=lsqnonlin(fit_fn,P0(1:6),[],[],options);
    P(7)=rot_angle;
else
    fit_fn= @(param) FitGauss2D_rot(param,data,X,Y);
    P=lsqnonlin(fit_fn,P0,[],[],options);
end
data_fit=imrotate(gauss2d(P,X,Y),-P(end),'bilinear','crop');
if ceil(P(3)) > size(data_fit,2)
    P_y_int=size(data_fit,2);
elseif ceil(P(3)) < 1
    P_y_int=1;
else
    P_y_int=ceil(P(3));
end
vertical_cut(:,1)=data_fit(:,P_y_int);

if ceil(P(4)) > size(data_fit,1)
    P_x_int=size(data_fit,1);
elseif ceil(P(4)) < 1
    P_x_int=1;
else
    P_x_int=ceil(P(4));
end

horizontal_cut(:,1)=data_fit(P_x_int,:);
data_rot=imrotate(data,P(end),'bilinear','crop');
horizontal_cut(:,2)=data_rot(P_x_int,:);
vertical_cut(:,2)=data_rot(:,P_y_int);
% disp(['Sigma X is ' num2str(P(5)*pixelsize*resize_factor) ' mm.']);
% disp(['Sigma Y is ' num2str(P(6)*pixelsize*resize_factor) ' mm.']);
% P0(3:4)*pixelsize*resize_factor

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
