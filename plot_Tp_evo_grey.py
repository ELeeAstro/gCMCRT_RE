import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

dire = 'results/'

fname = dire+'Tp_grey_1D_pp.txt'
f = open(fname,'r')

nit, nlay =  f.readline().split()
nit = int(nit)
nlay = int(nlay)
nlev = nlay + 1

ni = np.zeros(nit)
l = np.zeros((nit,nlay))
pl = np.zeros((nit,nlay))
t = np.zeros((nit,nlay))
Tl = np.zeros((nit,nlay))

for i in range(nit):
  ni[i] = f.readline()
  for n in range(nlay):
   l[i,n], pl[i,n], t[i,n], Tl[i,n] = f.readline().split()

f.close()

Tint = 1000.0
k_IR = 5.5e-3
grav = 1000.0

Tl_M = np.zeros(nlay)
Tl_M[:] = 3.0/4.0*Tint**4 * (t[-1,:] + 2.0/3.0)
Tl_M[:] = Tl_M[:]**(1.0/4.0)

fig = plt.figure()

col = sns.color_palette("husl", nit)

plt.plot(Tl[0,:], pl[0,:],c='darkcyan',label='IC',ls='solid')
for i in range(1,nit-1):
  plt.plot(Tl[i,:], pl[i,:],c=col[i],ls='dotted')
plt.plot(Tl[-1,:], pl[-1,:],c='darkmagenta',label='Last iter.',ls='solid')

plt.plot(Tl_M[:], pl[-1,:],c='black',label='Milne sol.',ls='dashed')

plt.xlabel(r'T$_{\rm gas}$ [K]',fontsize=16)
plt.ylabel(r'p$_{\rm gas}$ [bar]',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)

plt.yscale('log')
plt.legend()

plt.gca().invert_yaxis()

yticks = [10,1,0.1,0.01,1e-3,1e-4]
yticks_lab = [r'10',r'1',r'0.1',r'0.01',r'10$^{-3}$',r'10$^{-4}$']
plt.yticks(yticks,yticks_lab)

xticks = np.arange(750, 1050, step=50)
xticks_lab = ["%d" % xticks for xticks in xticks]
plt.xticks(xticks,xticks_lab)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)
#plt.savefig('MCRT_RCE_semi_grey_Tp.pdf',dpi=144,bbox_inches='tight')

plt.show()
