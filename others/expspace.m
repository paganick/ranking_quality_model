function v = expspace(vmin,vmax,inc) 
    v= vmin * inc.^(0:( log(vmax/vmin)/log(inc)));
    v = cast(v, 'int32');
end