#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the 2-panel triplet figure from ../data/triplet_data.csv (S9 + S10).
Left: the non-monotone omega-shape (peak at omega=3), model vs measured on S10.
Right: residuals on both shells; the omega=6 point tightens S9->S10 (small-sample),
omega=7 (11 triplets) is small-sample.
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('../data/triplet_data.csv')))
def series(sh):
    om=[int(r['omega']) for r in rows if r['shell']==sh]
    me=[float(r['measured']) for r in rows if r['shell']==sh]
    mo=[float(r['model']) for r in rows if r['shell']==sh]
    tr=[int(r['triplets']) for r in rows if r['shell']==sh]
    return map(np.array,(om,me,mo,tr))
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,5.4))
om,me,mo,tr=series('S10')
ax1.plot(om,me,'o-',color='#c0392b',lw=2.2,ms=9,label='measured (S10)',zorder=4)
ax1.plot(om,mo,'s--',color='#185FA5',lw=1.8,ms=7,label='two-centre CRT model',zorder=3)
ax1.axhline(1,color='gray',ls=':',lw=1)
ax1.annotate('$\\omega=7$: 11 triplets\n(small sample)',xy=(7,me[-1]),xytext=(5.0,0.40),
             fontsize=8.5,color='gray',arrowprops=dict(arrowstyle='->',color='gray',lw=1))
ax1.set_xlabel(r'$\omega_{>3}(N)$',fontsize=11)
ax1.set_ylabel(r'triplet rate, normalised to $\omega=1$',fontsize=11)
ax1.set_title('Prime triplet $(6N{-}1,6N{+}1,6N{+}5)$: non-monotone\n$\\omega$-shape, peak at $\\omega=3$ (S10)',fontsize=12)
ax1.legend(fontsize=9,loc='lower left'); ax1.grid(alpha=.25); ax1.set_xticks(range(1,8))
om9,me9,mo9,tr9=series('S9')
err9=100*(mo9-me9)/me9; err10=100*(mo-me)/me
ax2.axhspan(-3,3,color='#2ca25f',alpha=.10,label='$\\pm3\\%$')
ax2.axhline(0,color='gray',lw=1)
ax2.plot(om9,err9,'o-',color='#888',lw=1.4,ms=6,label='S9',alpha=.7)
ax2.plot(om,err10,'s-',color='#c0392b',lw=1.8,ms=7,label='S10')
ax2.annotate('$\\omega=6$: $6.5\\%\\!\\to\\!0.5\\%$\n(56$\\to$1832 triplets)',xy=(6,-0.5),xytext=(3.2,-14),
             fontsize=8.5,color='#185FA5',arrowprops=dict(arrowstyle='->',color='#185FA5',lw=1))
ax2.set_xlabel(r'$\omega_{>3}(N)$',fontsize=11)
ax2.set_ylabel('model $-$ measured (%)',fontsize=11)
ax2.set_title('Residuals $\\leq2.4\\%$ for $\\omega\\leq6$;\n$\\omega=7$ ($11$ triplets) small-sample',fontsize=12)
ax2.legend(fontsize=9); ax2.grid(alpha=.25); ax2.set_xticks(range(1,8)); ax2.set_ylim(-25,6)
plt.suptitle('The prime triplet as a three-member instance of the two-centre CRT mechanism (Paper X, extended)',fontsize=12.5,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper15_triplet.pdf',bbox_inches='tight')
plt.savefig('fig_paper15_triplet.png',dpi=160,bbox_inches='tight')
print("figure saved")
