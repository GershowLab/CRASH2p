function it = polarinterp1(x, theta, xi, varargin)
%functionit = polarinterp1(x, theta, xi, varargin)
%
%projects all theta onto the unit circle, interpolates projections, finds
%angle
    xx = interp1(x,cos(theta),xi,varargin{:});
    yy = interp1(x,sin(theta),xi,varargin{:});
    it = atan2(yy,xx);
end