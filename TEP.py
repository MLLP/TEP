import numpy as np; import numpy.matlib as mtlib; import plotly.figure_factory as ff

def BRK(Mx, vA):
    dx = Mx[:,1] - Mx[:,0]; sumEP = 0;
    for i in range(0,len(vA)):
        sumEP = sumEP + vA[i]*dx**i;
    return  Mx[:,0] * Mx[:,1]  * sumEP

def TEP(par, model, Mx, mB):
    #EP123
    if model[0]==0:
        EP123 = 0
    elif model[0]==1:
        EP12 = BRK(Mx[:,[0,1]], mB[0,:]); EP13 = BRK(Mx[:,[0,2]], mB[1,:]); EP23 = BRK(Mx[:,[1,2]], mB[2,:])
        EP123 = EP12 + EP13 + EP23
    elif model[0]==5:
        EP12 = BRK(Mx[:,[0,1]]/mtlib.repmat(Mx[:,[0,1]].sum(1),2,1).T, mB[0,:])
        EP13 = BRK(Mx[:,[0,2]]/mtlib.repmat(Mx[:,[0,2]].sum(1),2,1).T, mB[1,:])
        EP23 = BRK(Mx[:,[1,2]]/mtlib.repmat(Mx[:,[1,2]].sum(1),2,1).T, mB[2,:])
        EP123 = (Mx[:,[0,1]].sum(1)*EP12 + Mx[:,[0,2]].sum(1)*EP13 + Mx[:,[1,2]].sum(1)*EP23)/2.
    
    if model[1] == 0:
        EPcalc = EP123
    elif model[1] == 1:
        EPcalc = EP123 + Mx[:,0]*Mx[:,1]*Mx[:,2]*par
    elif model[1] == 2:
        EPcalc = EP123 + Mx[:,0]*Mx[:,1]*Mx[:,2]*(par[0]+par[1]*Mx[:,0]+par[2]*Mx[:,1])
    elif model[1] == 4:
        dx21 = Mx[:,1] - Mx[:,0]; dx32 = Mx[:,2] - Mx[:,1]; 
        EPcalc = EP123 + Mx[:,0]*Mx[:,1]*Mx[:,2]*(par[0]+par[1]*dx21+par[2]*dx32+par[3]*dx21**2+par[4]*dx32**2)
    elif model[1] == 5:
        dx21 = Mx[:,1] - Mx[:,0]; dx32 = Mx[:,2] - Mx[:,1];
        soma4 = par[0]+par[1]*dx21+par[2]*dx32+par[3]*dx21**2+par[4]*dx32**2 
        EPcalc = EP123 + Mx[:,0]*Mx[:,1]*Mx[:,2]*(soma4 + par[5]*dx21**3 + par[6]*dx32**3)
    elif model[1] == 6:
        dx21 = Mx[:,1] - Mx[:,0]; dx32 = Mx[:,2] - Mx[:,1]; 
        soma5 = par[0]+par[1]*dx21+par[2]*dx32+par[3]*dx21**2+par[4]*dx32**2+par[5]*dx21**3+par[6]*dx32**3
        EPcalc = EP123 + Mx[:,0]*Mx[:,1]*Mx[:,2]*(soma5 + par[7]*dx21**4 + par[8]*dx32**4)
    elif model[1] == 7:
        dx21 = Mx[:,1] - Mx[:,0]; dx32 = Mx[:,2] - Mx[:,1]; 
        soma6 = par[0]+par[1]*dx21+par[2]*dx32+par[3]*dx21**2+par[4]*dx32**2+par[5]*dx21**3+par[6]*dx32**3+par[7]*dx21**4+par[8]*dx32**4
        EPcalc = EP123 + Mx[:,0]*Mx[:,1]*Mx[:,2]*(soma6 + par[9]*dx21**5 + par[10]*dx32**5)    
    return EPcalc

vi = np.linspace(1e-10,1-1e-10,21); vj = np.linspace(1e-10,1-1e-10,21); Mxter = [];
for i in vi:
    for j in vj:
        k = 1-i-j
        if k >= 0.:
            Mxter += [[i, j, k]];
Mxter = np.array(Mxter);

subst1 = '[BMIM][PF6]'; subst2 = 'acetic acid'; subst3 = 'acetophenone'; lmix = [subst1,subst2,subst3];
mixture = str(subst1)+'/'+str(subst2)+'/'+str(subst3)+' at 30°C'
mB = np.array([[475358.71, 228492.63, 290691.21, -27347.978, -350571.44],[138325, -290681, -137488, 57475, -64931],[463197.2, 290965.01, -991881.11, 79472.753, 193296.52]])
EPcalc = np.array(TEP(np.array([-25.084695,   17.113304,  -23.183305,  -72.389953,   10.339557,   169.13194,   22.111202,  0,  0,  0,  0])*1e5, [5, 7], Mxter, mB));

fig = ff.create_ternary_contour(Mxter.T,EPcalc,title=mixture,pole_labels=lmix,interp_mode='ilr',showscale='True',ncontours=40,coloring='lines')
fig.show()
