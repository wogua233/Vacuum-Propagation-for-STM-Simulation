"""One physical quantity: boundary shape-jet mismatch modulo a common envelope.

F/Fz/Fzz are native DFT source values and normal derivatives. A/Az/Azz are
an explicitly specified radial-collapsed atomic ansatz on that same plane.
No VP, Green-function transfer, FFT vacuum propagation or tip-height image
is evaluated. The scalar spectral visibility is an energy-only barrier
attenuation estimate, not a spatial wavefunction propagator.
"""
import numpy as np
KIN=3.80998212

def normalized_jet(F,Fz,Fzz,w,reg):
 I=np.einsum('n,pn->p',w,abs(F)**2)
 D=2*np.einsum('n,pn->p',w,(F.conj()*Fz).real)
 DD=2*np.einsum('n,pn->p',w,abs(Fz)**2+(F.conj()*Fzz).real)
 mean=float(reg@I)
 if mean<=1e-250:raise ValueError('Source intensity vanishes')
 J=I/mean;Jz=(D-(reg@D)*J)/mean;Jzz=(DD-(reg@DD)*J-2*(reg@D)*Jz)/mean
 return J,Jz,Jzz

def score(F,Fz,Fzz,A,Az,Azz,weights,region,energies,vacuum_reference,distance,length=None,attenuate=True):
 if distance<0:raise ValueError('Observation below the source plane')
 reg=np.asarray(region,float);reg=reg/reg.sum();w=np.asarray(weights,float).copy()
 if attenuate:
  barrier=vacuum_reference-np.asarray(energies)
  if np.any(barrier[w>0]<=0):raise ValueError('Open mode; scalar evanescent visibility is inapplicable')
  k=np.sqrt(np.maximum(barrier,0)/KIN)
  active=w>0;w*=np.exp(-2*distance*(k-k[active].min()))
 J,Jz,Jzz=normalized_jet(F,Fz,Fzz,w,reg);K,Kz,Kzz=normalized_jet(A,Az,Azz,w,reg)
 L=distance/2 if length is None else length
 residual=J-K+L*(Jz-Kz)+.5*L*L*(Jzz-Kzz)
 val=float(np.sqrt(reg@(residual**2)/(reg@(J**2))))
 return val,dict(source=J-K,slope=L*(Jz-Kz),curvature=.5*L*L*(Jzz-Kzz),residual=residual,native_shape=J,reduced_shape=K,weights=w,length_A=float(L))

def harmonics(dr):
 R=np.linalg.norm(dr,axis=-1);ir=1/R;u=dr/R[...,None];x,y,z=[u[...,k] for k in range(3)]
 Y=np.stack([np.ones_like(x),1.7320508*y,1.7320508*z,1.7320508*x,3.8729833*x*y,3.8729833*y*z,1.1180300*(3*z*z-1),3.8729833*x*z,3.8729833/2*(x*x-y*y)],axis=-1)
 Yz=np.stack([np.zeros_like(x),-1.7320508*y*z,1.7320508*(1-z*z),-1.7320508*x*z,-2*3.8729833*x*y*z,3.8729833*y*(1-2*z*z),6*1.1180300*z*(1-z*z),3.8729833*x*(1-2*z*z),-3.8729833*z*(x*x-y*y)],axis=-1)*ir[...,None]
 Yzz=np.stack([np.zeros_like(x),1.7320508*y*(3*z*z-1),-3*1.7320508*z*(1-z*z),1.7320508*x*(3*z*z-1),2*3.8729833*x*y*(4*z*z-1),2*3.8729833*y*z*(4*z*z-3),6*1.1180300*(1-z*z)*(1-4*z*z),2*3.8729833*x*z*(4*z*z-3),3.8729833*(x*x-y*y)*(4*z*z-1)],axis=-1)*ir[...,None]**2
 return R,z,Y,Yz,Yzz

def ansatz_jet(points,centres,C,k):
 dr=points[:,None]-centres[None];R,z,Y,Yz,Yzz=harmonics(dr);decay=np.exp(-k*R)
 angular=np.einsum('pal,aln->pan',Y,C);dz=np.einsum('pal,aln->pan',Yz,C);dzz=np.einsum('pal,aln->pan',Yzz,C)
 A=np.einsum('pa,pan->pn',decay,angular)
 Az=np.einsum('pa,pan->pn',decay,dz-k*z[:,:,None]*angular)
 Azz=np.einsum('pa,pan->pn',decay,dzz-2*k*z[:,:,None]*dz+(k*k*z*z-k*(1-z*z)/R)[:,:,None]*angular)
 return A,Az,Azz


def effective_kappas(energies,vacuum_reference,source_norm,transverse_norm,weights,cluster_tolerance_eV=1e-5):
 """One scalar normal decay per energy subspace, including transverse kinetic energy.

 Degenerate-state trace averages make the rate invariant to basis rotation
 within that energy subspace. No transverse Fourier transfer is performed.
 """
 E=np.asarray(energies);S=np.asarray(source_norm);G=np.asarray(transverse_norm);w=np.asarray(weights);order=np.argsort(E);q=np.zeros_like(E,float);start=0
 while start<len(order):
  stop=start+1
  while stop<len(order) and E[order[stop]]-E[order[start]]<=cluster_tolerance_eV:stop+=1
  idx=order[start:stop];den=float(w[idx]@S[idx]);value=float(w[idx]@G[idx])/max(den,1e-280);q[idx]=max(value,0.);start=stop
 return np.sqrt(np.maximum((vacuum_reference-E)/KIN,0)+q),q
