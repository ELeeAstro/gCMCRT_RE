import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

dire = 'results/'

fname = dire+'estim_pp_grey_1D_pp.txt'
data = np.loadtxt(fname)

n = data[:,0]
tau_l = data[:,1]
Jdot = data[:,2]
Hdot = data[:,3]
Kdot = data[:,4]
Adot = data[:,5]

e1 = Hdot/Jdot
e2 = Kdot/Hdot
e3 = Kdot/Jdot

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(e1,tau_l,label=r'e1 - H/J',ls='dashed',c=col[0])
#plt.plot(e2,tau_l,label=r'e2 - K/H',ls='dashed',c=col[1])
plt.plot(e3,tau_l,label=r'e3 - K/J',ls='dashed',c=col[3])

plt.vlines(1.0/3.0,tau_l[0],tau_l[-1],ls='dotted',colors=col[0])
#plt.vlines(0.41,tau_l[0],tau_l[-1],ls='dotted',colors=col[0])

plt.vlines(0.5,tau_l[0],tau_l[-1],ls='dotted',colors=col[1])
plt.vlines(1.0/np.sqrt(3.0),tau_l[0],tau_l[-1],ls='dotted',colors=col[1])

plt.text(0.31,40,r'1/3',c=col[0],fontsize=12)


plt.text(0.48,40,r'1/2',c=col[1],fontsize=12)
plt.text(0.55,0.1,r'1/$\sqrt{3}$',c=col[1],fontsize=12)

plt.legend(title=r'Ratios',loc='lower right')

plt.xlim(0.3,0.62)

plt.yscale('log')
plt.gca().invert_yaxis()

yticks = [1e2,1e1,1e0,1e-1,1e-2,1e-3,1e-4]
yticks_lab = [r'10$^{2}$',r'10$^{1}$',r'10$^{0}$',r'10$^{-1}$',r'10$^{-2}$',r'10$^{-3}$',r'10$^{-4}$']
plt.yticks(yticks,yticks_lab)

plt.xlabel(r'Ratio',fontsize=16)
plt.ylabel(r'$\tau_{\rm R}$',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)
#plt.savefig('MCRT_RCE_semi_grey_est.pdf',dpi=144,bbox_inches='tight')

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(1.0/e1,tau_l,label=r'D',ls='dashed',c=col[1])

plt.xlim(1,2)

#plt.legend(title=r'Diffusion Coefficents',loc='lower right')

plt.yscale('log')
plt.gca().invert_yaxis()

yticks = [1e2,1e1,1e0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6]
yticks_lab = [r'10$^{2}$',r'10$^{1}$',r'10$^{0}$',r'10$^{-1}$',r'10$^{-2}$',r'10$^{-3}$',r'10$^{-4}$',r'10$^{-5}$',r'10$^{-6}$']
plt.yticks(yticks,yticks_lab)

plt.xlabel(r'Diffusion COefficent',fontsize=16)
plt.ylabel(r'$\tau_{\rm R}$',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)
#plt.savefig('MCRT_RCE_semi_grey_est.pdf',dpi=144,bbox_inches='tight')

plt.show()

