function yn = gfun_ssim_v13(xn, un, n, p, ny, nevals) %Returns y(n) from x(n)
d       = p(end);
yn      = nan(ny, nevals);
z       = xn(1,:) + xn(2,:) + xn(3,:) + 0*xn(4,:);
yn(1,:) = d*un(1,:) - z;     %Motor error for current step (~SLA)
end