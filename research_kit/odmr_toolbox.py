import os
import WrightTools as wt
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
ls = sp.optimize.least_squares
beta = 13.9902  # GHz/T
beta *= 1e-3  # GHZ/mT
g = 2.00232

con = np.concatenate

def Bactual(BsetkG):
    Bset = BsetkG * 1e3
    Bactual =  2e-10*Bset**3 - 8e-6*Bset**2 + 1.0869*Bset-3.6262
    return Bactual * 1e-3

def data_boxcar(d, width):
    Bold, Xold, Yold, Zold = d.Bfield.points, d.X.points, d.Y.points, d.AUX1.points
    num = Bold.size // width
    Bnew, Xnew, Ynew, Znew = np.zeros(num), np.zeros(num), np.zeros(num), np.zeros(num)
    for i in range(num):
        s = slice(width * i, width * (i + 1))
        Bnew[i], Xnew[i], Ynew[i], Znew[i] = np.mean(Bold[s]), np.mean(Xold[s]), np.mean(Yold[s]), np.mean(Zold[s])
    return Bnew, Xnew, Ynew, Znew

def col_boxcar(c, width, newname, newparent, mulfactor=-3.8e9, offset=-0.016, mT=False):
    items = c.item_names
    ns = []
    for item in items:
        if 'Sweep' in item:
            ns.append(item)
    # mulfactor is in A/V
    Bs, Xs, Ys, Zs = [], [], [], []
    for n in ns:
        d = c[n]
        Bs.append(d.B.points)
        Xs.append(d.X.points)
        Ys.append(d.Y.points)
        Zs.append((d.AUX1.points+offset)/mulfactor)
    Bold, Xold, Yold, Zold = con(Bs), con(Xs), con(Ys), con(Zs)   
    Xold /= Zold
    Yold /= Zold
    indx = np.argsort(Bold)
    Bold, Xold, Yold, Zold =  Bold[indx], Xold[indx], Yold[indx], Zold[indx]  
    num = Bold.size // width
    Bnew, Xnew, Ynew, Znew = np.zeros(num), np.zeros(num), np.zeros(num), np.zeros(num)
    for i in range(num):
        s = slice(width * i, width * (i + 1))
        Bnew[i], Xnew[i], Ynew[i], Znew[i] = np.mean(Bold[s]), np.mean(Xold[s]), np.mean(Yold[s]), np.mean(Zold[s])
    d = wt.Data(name=newname, parent=newparent)
    if mT:
        d.create_variable('B', values=Bactual(Bnew)*100, units='mT')
    else:
        d.create_variable('B', values=Bactual(Bnew), units='kG')
    d.transform('B')
    d.create_channel('X', values=Xnew)
    d.create_channel('Y', values=Ynew)
    d.create_channel('Z', values=Znew)
    d.create_channel('R', values=np.sqrt(Xnew**2 + Ynew**2))
    return d

def data_linear_offset(d, numpts, smooth=None, copy=False):
    if copy:
        d = d.copy(verbose=False)
    if smooth != None:
        d.smooth(smooth, verbose=False)
    x = d.B.points 
    X = d.X.points
    Y = d.Y.points
    x1, x2 = np.mean(x[:numpts]), np.mean(x[-1*numpts:])
    y1, y2 = np.mean(X[:numpts]), np.mean(X[-1*numpts:])
    m = (y2-y1)/(x2-x1) 
    Xlin = m*x+ (y1-m*x1)
    y1, y2 = np.mean(Y[:numpts]), np.mean(Y[-1*numpts:])
    m = (y2-y1)/(x2-x1) 
    Ylin = m*x+ (y1-m*x1)
    X -= Xlin
    Y -= Ylin    
    R = np.sqrt(X**2 + Y**2)
    d.X[:] = X
    d.Y[:] = Y
    d.R[:] = R
    return d
    
def mod_single_freq_sim(rs, rt, ds, dt, Gs, Gt, rsnr, rtnr, kISC, a, F, f):
    # implementation of DOI: 10.1103/PhysRevB.86.115204
    # all parameters other than F are in units of 1/s
    # unable to reproduce some figures e.g. 6b, 8d
    # returns the inphase and out-of-phase signal at a single frequency points 
    def IcIs(rs, rt, ds, dt, Gs, Gt, rsnr, rtnr, kISC, a, F, f, ell):
        T = 1/f
        Cs = rs + rsnr + ds
        Ct = rt + rtnr + dt
        w11 = a + kISC*(1-F)
        w21 = a + kISC*F
        w12 = kISC*(1-F)
        w22 = kISC*F
        
        ns1_0 = (w21*Gt + (Ct+w21)*Gs) / ((Cs+w11)*(Ct+w21)-w11*w21)
        ns2_0 = (w22*Gt + (Ct+w22)*Gs) / ((Cs+w12)*(Ct+w22)-w12*w22)
        nt1_0 = (w11*Gs + (Cs+w11)*Gt) / ((Cs+w11)*(Ct+w21)-w11*w21)
        nt2_0 = (w12*Gs + (Cs+w12)*Gt) / ((Cs+w12)*(Ct+w22)-w12*w22)
        
        argj1 = np.sqrt((Cs+w11-Ct-w21)**2 + 4*w11*w21)
        argj2 = np.sqrt((Cs+w12-Ct-w22)**2 + 4*w12*w22)
        m11 = (Cs+w11+Ct+w21-argj1)/2
        m21 = (Cs+w11+Ct+w21+argj1)/2
        m12 = (Cs+w12+Ct+w22-argj2)/2
        m22 = (Cs+w12+Ct+w22+argj2)/2
        
        b11 = (Cs+w11-m11)/w21
        b21 = (Cs+w21-m21)/w11 #b11 different than paper 
        b12 = (Cs+w12-m12)/w21
        b22 = (Cs+w22-m22)/w11 #b12 different than paper 
        g11 = np.exp(-1*m11*T/2)
        g21 = np.exp(-1*m21*T/2)
        g12 = np.exp(-1*m12*T/2)
        g22 = np.exp(-1*m22*T/2)
        
        dns0 = ns2_0 - ns1_0
        dnt0 = nt2_0 - nt1_0
        Bmat = np.array([[1, 1, -1*g12, -1*g22],
                         [b11, b21, -1*b12*g12, -1*b22*g22],
                         [g11, g21, -1, -1],
                         [b11*g11, b21*g21, -1*b12, -1*b22]])
        nvec = np.array([dns0, dnt0, dns0, dnt0])
        
        Avec = np.linalg.solve(Bmat, nvec) # paper writes out full expression 
        A11, A21, A12, A22 = Avec
        B11, B21, B12, B22 = A11*b11, A21*b21, A12*b12, A22*b22
        
        # time dependance not done here
        
        c0 = np.cos(ell*np.pi)
        c1 = 4*ell**2*np.pi**2/T**2
        Ic = 2*f*(m11*(rs*A11+rt*B11)*(1-g11*c0)/(m11**2+c1)
                  + m21*(rs*A21+rt*B21)*(1-g21*c0)/(m21**2+c1)
                  + m12*(rs*A12+rt*B12)*(c0-g12)/(m12**2+c1)
                  + m22*(rs*A22+rt*B22)*(c0-g22)/(m22**2+c1))
        Is = 4*ell*np.pi/T**2 *((rs*A11+rt*B11)*(1-g11*c0)/(m11**2+c1)
                  + (rs*A21+rt*B21)*(1-g21*c0)/(m21**2+c1)
                  + (rs*A12+rt*B12)*(c0-g12)/(m12**2+c1)
                  + (rs*A22+rt*B22)*(c0-g22)/(m22**2+c1)) + (rs*dns0 + rt*dnt0)*(c0-1)/ell/np.pi 
        return Ic, Is
    Ic1, Is1 = IcIs(rs, rt, ds, dt, Gs, Gt, rsnr, rtnr, kISC, a, F, f, 1)
    return Is1/2, Ic1/2    

def mod_freqs_sim(rs, rt, ds, dt, Gs, Gt, rsnr, rtnr, kISC, a, F, fs):
    Vis,Vos = np.zeros(fs.shape), np.zeros(fs.shape)
    for i,f in enumerate(fs):
        Vi,Vo = mod_single_freq_sim(rs, rt, ds, dt, Gs, Gt, rsnr, rtnr, kISC, a, F, f)
        Vis[i], Vos[i] = Vi, Vo    
    return Vis, Vos
    



# --- modeling functions for EPR lineshapes  ------------------------------------------------------

beta = 13.9902  # GHz/T
beta *= 1e-3  # GHZ/mT
g = 2.00232


def LS(x, x0, G, which="L"):
    # area normalized, G is FWHM
    if which == "L":
        out = G / (2 * np.pi) * ((x - x0) ** 2 + (G / 2) ** 2) ** -1
    elif which == "G":
        s = G / (2 * np.sqrt(2 * np.log(2)))
        out = (s * np.sqrt(2 * np.pi)) ** -1 * np.exp(
            -(((x - x0) / (np.sqrt(2) * s)) ** 2)
        )
    else:
        out = np.zeros(x.shape)
    return out


def Bresonance_triplet(D0GHz, E0GHz, thetarad, phirad, RFGHz=0):
    # returns resonance centers in mT
    D = D0GHz/2 * (3 * np.cos(thetarad) ** 2 - 1)
    E = E0GHz/2 * 3 * np.sin(thetarad) ** 2 * np.cos(2 * phirad)
    pre = 1 / (g * beta)
    one = np.sqrt((RFGHz - D - E) ** 2)
    two = np.sqrt((RFGHz + D + E) ** 2)
    return pre * one, pre * two


def ordering(L, thetarad, norm="sine"):
    # This function is normalized.
    x = L / 4
    if norm == "flat":
        c = (np.pi * np.exp(-x) * sp.special.iv(0, x)) ** -1
    elif L == 0:
        c = 1 / 2
    elif norm == "sine":
        imag = sp.special.erf(np.lib.scimath.sqrt(-L / 2)) / np.lib.scimath.sqrt(-L)
        c = (np.sqrt(2 * np.pi) * np.exp(-L / 2) * imag) ** -1
        c = np.abs(c)  # is this what I want?
    else:
        c = 1
    return c * np.exp(L / 2 * (np.cos(thetarad) ** 2 - 1))


def single_triplet(
    BmT, thetarad, phirad, D0GHz, E0GHz, GmT, weight=1, order=0, RFGHz=0, which="L"
):
    # no volume integrals or sine weighting
    c1, c2 = Bresonance_triplet(D0GHz, E0GHz, thetarad, phirad, RFGHz)
    out = weight * LS(BmT, c1, GmT, which) + weight * LS(BmT, c2, GmT, which)
    out *= ordering(order, thetarad)
    return out


def integrate_single_triplet(
    BmT, thetarad, phirad, D0GHz, E0GHz, GmT, weight=1, order=0, RFGHz=0, which="L"
):
    c1, c2 = Bresonance_triplet(D0GHz, E0GHz, thetarad[:, None], phirad[None, :], RFGHz=RFGHz)
    out = LS(BmT[:, None, None], c1[None, ...], GmT, which=which) + LS(
        BmT[:, None, None], c2[None, ...], GmT, which=which
    )
    out *= np.sin(thetarad[None, :, None])
    out = np.trapz(out, x=phirad, axis=-1)
    orderweight = weight * ordering(order, thetarad, norm="sine")
    out *= orderweight[None, :]
    out = np.trapz(out, x=thetarad, axis=-1)
    return out


def integrate_multiple_triplet(
    BmT, thetarad, phirad, D0GHzs, E0GHzs, GmTs, weights, order=0, RFGHz=0, which="L"
):
    B, t, p = (
        BmT[:, None, None, None],
        thetarad[None, :, None, None],
        phirad[None, None, :, None],
    )
    D, E, G, W = (
        D0GHzs[None, None, None, :],
        E0GHzs[None, None, None, :],
        GmTs[None, None, None, :],
        weights[None, None, None, :],
    )
    c1, c2 = Bresonance_triplet(D, E, t, p, RFGHz=RFGHz)
    out = W * LS(B, c1, G, which=which) + W * LS(B, c2, G, which=which)
    out *= np.sin(t) * ordering(order, t, norm="sine")
    out = np.sum(out, axis=-1)
    out = np.trapz(out, x=phirad, axis=-1)
    out = np.trapz(out, x=thetarad, axis=-1)
    return out

def multi_T(p,x, num=3, thetanum=41, phinum=31, which='L'):
    # assumes centered resonances in highfield limit
    theta = np.linspace(0,np.pi,thetanum)
    phi = np.linspace(0,2*np.pi,phinum)    
    RFGHz = 9
    Bcenter = RFGHz /(g*beta)  
    B = x + Bcenter  
    order = p[0]
    Ds = p[1:num+1]
    Gs = p[num+1:2*num+1]
    Ws = p[2*num+1:3*num+1]
    Es = p[3*num+1:4*num+1]
    out = integrate_multiple_triplet(B, theta, phi, Ds, Es, Gs, Ws, order, RFGHz, which)
    return out

def multi_T_fit(pguess, xdata, ydata, pbounds=(-np.inf,np.inf), num=3, thetanum=41, phinum=31, which='L'):
    def cost(p):
        return ydata - multi_T(p, xdata, num, thetanum, phinum, which)
    fit = ls(cost, pguess, bounds=pbounds)
    return fit.x

"""
def multi_T_order_widths_weights(pold, xdata, ydata, num=3, thetanum=41, phinum=31, which='L'):
    def cast(p):
        pnew = pold.copy()
        pnew[0] = p[0]
        pnew[num+1:2*num+1] = p[1:num+1]
        pnew[2*num+1:3*num+1] = p[num+1:2*num+1]
        return pnew        
    def cost(p):
        pnew = cast(p)
        return ydata - multi_T(pnew, xdata, num, thetanum, phinum, which)
    pguess = np.concatenate((np.array([pold[0]]), pold[1:num+1], pold[num+1:2*num+1]),bounds=pbounds)
    fit = ls(cost, pguess)
    return cast(fit.x)
"""