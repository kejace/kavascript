String.prototype.format = function() {
    var formatted = this
    for (var i = 0; i < arguments.length; i++) {
        var regexp = new RegExp('\\{'+i+'\\}', 'gi')
        formatted = formatted.replace(regexp, arguments[i])
    }
    return formatted
}

function randomRange(min, max) { 
	return (min + (Math.random()* (max-min)))
}

function rk4iter(x,y,h)
{
    // compute the four RK coefficients using the RHS of the diff eq.
    k1 = f(x,y)
    k2 = f(x + h/2, y + k1 * h / 2)
    k3 = f(x + h/2, y + k2 * h / 2)
    k4 = f(x + h, y + k3 * h)
    
    // compute the y value at the end of the iteration
    return y + (k1 + 2*k2 + 2*k3 + k4) * h / 6
}

function rk4(n,y0,x1) 
{
    var h = x1/n   // compute the step size along the x-axis
    
    var x = 0      // initial value of x is 0
    var y = (+y0)  // initial value of y is y0
    
    g2 = [[x,y]]   // start with initial conditions
    
    // loop until x reaches its target value
    for (var i = 1; i <= n; i++)
    {
        // Perform iterations of the RK4 algorithm
        y = rk4iter(x,y,h)
        x += h
        
        // Add data to the plot as well.
        g2.push([x,y])
    }
    
    
    return g2
}



function inner_product(psi_r1, psi_i1, psi_r2, psi_i2) {
    var n = psi_r1.length;
    tp = 0;
    tz = 0;
    for (i = 0; i<n; ++i) {
        var z = zScaler(i);
        var p = psi_r1[i]*psi_r2[i]+psi_i1[i]*psi_i2[i];
        tp += p;
        tz += p*z;
    }
    return tz/tp;
}


    function fft(yr, yi) {
        with (Math) {
            n = yr.length;
            g = floor(log(n)/log(2)+0.5);
            p = 2*PI/n;
        
            for (l = 0; l < g; l++) {
                g1 = pow(2, g-l-1);
                var m = 0;
                
                for (i = 0; i < pow(2, l); i++) {
                    k1 = floor(m/g1);
                    k2 = twiddle(k1);
                    y1 = table_c[k2];
                    y2 = -table_s[k2];
                            
                    for (j = 0; j < g1; j++)  {
                        mg1 = m+g1;
                        y3 = yr[mg1]*y1-yi[mg1]*y2;
                        y4 = yr[mg1]*y2+yi[mg1]*y1;
                        yr[mg1] = yr[m]-y3;
                        yi[mg1] = yi[m]-y4;
                        yr[m] = yr[m]+y3;
                        yi[m] = yi[m]+y4;
                        m = m+1;
                    }
                    m = m + g1;
                }
            }
        
        for (i = 0; i < n; i++) {
                k1 = i;
                k2 = twiddle(k1);
                
                if (k2 >= i) continue;
                                
                k3 = yr[i];
                yr[i] = yr[k2];
                yr[k2] = k3;
                k3 = yi[i];
                yi[i] = yi[k2];
                yi[k2] = k3;
            }
        }
    }

    function ifft(yr, yi) {
        fft(yr, yi);
     
        n = yr.length;
        recip = 1/n;
        for (i = 0; i < n; i++) {
            yr[i] = recip*yr[i];
            yi[i] = recip*yi[i];
        }
        for (i = 1; i < Math.floor(n/2); i++) {
            t = yr[n-i];
            yr[n-i] = yr[i];
            yr[i] = t;
            t = yi[n-i];
            yi[n-i] = yi[i];
            yi[i] = t;
        }
    }


function twiddle(k1) {
    k2 = 0;
    
    for (k = 0; k<g; k++) {
        k3 = Math.floor(k1/2);
        k2 = 2*(k2-k3) + k1;
        k1 = k3;
    }
    return k2;
}

function range(n){
	vals = []
	for(var i = 0; i < n; i++){
		vals.push(i)
	}
	return vals
}