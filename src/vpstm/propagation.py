"""Forbidden-region Schrodinger boundary solver: periodic xy, Dirichlet z.

Atomic-unit equation: [-laplacian + 2*(U-E)] psi = 0. The prescribed
bottom wavefunction is never inferred from density. The top is zero and
must be moved away until predictions at the requested heights converge.
"""
import numpy as np
from scipy.fft import dst, idst, fft2, ifft2
from scipy.sparse.linalg import LinearOperator, cg, gmres
from .fields import ANG_TO_BOHR, EV_TO_HA


def propagate(boundary, potential_eV, energy_eV, spacing_A, tolerance=1e-9,
              maxiter=500, allow_allowed=False):
    u=np.asarray(potential_eV,dtype=float)
    b=np.asarray(boundary,dtype=np.complex128)
    if u.ndim!=3 or u.shape[0]<3 or u.shape[1:]!=b.shape:
        raise ValueError('Potential must be (nz,ny,nx), boundary (ny,nx)')
    if not np.all(np.isfinite(u)) or not np.all(np.isfinite(b)) or not np.isfinite(energy_eV):
        raise ValueError('Non-finite propagation data')
    dz,dy,dx=np.asarray(spacing_A,dtype=float)*ANG_TO_BOHR
    if min(dx,dy,dz)<=0:
        raise ValueError('Requires positive spacing')
    allowed = bool(np.min(u[1:-1]-energy_eV)<=0)
    if allowed and not allow_allowed:
        raise ValueError('Requires a forbidden interior unless allow_allowed=True')
    scale=float(np.max(abs(b)))
    if scale==0: return np.zeros_like(u),{'iterations':0,'relative_residual':0.0}
    n,ny,nx=u.shape[0]-2,u.shape[1],u.shape[2]
    q=2*EV_TO_HA*(u[1:-1]-energy_eV)
    rhs=np.zeros(q.shape,dtype=np.complex128);rhs[0]=b/scale/dz**2

    def matvec(vector):
        a=vector.reshape(q.shape)
        result=(2/dx**2+2/dy**2+2/dz**2+q)*a
        result-=(np.roll(a,1,axis=2)+np.roll(a,-1,axis=2))/dx**2
        result-=(np.roll(a,1,axis=1)+np.roll(a,-1,axis=1))/dy**2
        result[1:]-=a[:-1]/dz**2;result[:-1]-=a[1:]/dz**2
        return result.ravel()

    # Exact inverse of the constant-potential discrete operator. This retains
    # all transverse Fourier modes, unlike independent vertical WKB columns.
    kz=4*np.sin(np.pi*np.arange(1,n+1)/(2*(n+1)))**2/dz**2
    ky=4*np.sin(np.pi*np.arange(ny)/ny)**2/dy**2
    kx=4*np.sin(np.pi*np.arange(nx)/nx)**2/dx**2
    # A positive reference keeps the preconditioner nonsingular even for
    # locally allowed or indefinite domains; the actual operator is unchanged.
    reference_q=max(float(q.mean()),1e-4) if allowed else float(q.mean())
    denom=kz[:,None,None]+ky[None,:,None]+kx[None,None,:]+reference_q

    def precondition(vector):
        a=vector.reshape(q.shape)
        transformed=fft2(dst(a,type=1,axis=0,norm='ortho'),axes=(1,2))
        value=ifft2(transformed/denom,axes=(1,2))
        return idst(value,type=1,axis=0,norm='ortho').ravel()

    size=q.size
    operator=LinearOperator((size,size),matvec=matvec,dtype=float)
    pre=LinearOperator((size,size),matvec=precondition,dtype=float)
    count=[0]
    def callback(_): count[0]+=1
    if allowed:
        try:
            solution,info=gmres(operator,rhs.ravel(),M=pre,rtol=tolerance,atol=0,
                                restart=40,maxiter=maxiter,callback=callback,callback_type='legacy')
        except TypeError:
            solution,info=gmres(operator,rhs.ravel(),M=pre,tol=tolerance,atol=0,
                                restart=40,maxiter=maxiter,callback=callback,callback_type='legacy')
    else:
        try:
            solution,info=cg(operator,rhs.ravel(),M=pre,rtol=tolerance,atol=0,maxiter=maxiter,callback=callback)
        except TypeError:
            solution,info=cg(operator,rhs.ravel(),M=pre,tol=tolerance,atol=0,maxiter=maxiter,callback=callback)
    residual=float(np.linalg.norm(matvec(solution)-rhs.ravel())/np.linalg.norm(rhs))
    if info!=0 or residual>10*tolerance:
        raise RuntimeError('Propagation did not converge: info={}, residual={}'.format(info,residual))
    result=np.zeros(u.shape,dtype=np.complex128);result[0]=b;result[1:-1]=solution.reshape(q.shape)*scale
    return result,{'iterations':count[0],'relative_residual':residual,
                  'solver':'gmres' if allowed else 'cg',
                  'allowed_fraction':float(np.mean(q<=0)),
                  'minimum_barrier_eV':float(np.min(u[1:-1]-energy_eV))}


def height_slices(values,z,heights):
    out=[]
    for height in heights:
        if height<z[0] or height>z[-1]: raise ValueError('Requested height outside propagated slab')
        lower=min(int(np.searchsorted(z,height,side='right')-1),len(z)-2)
        t=(height-z[lower])/(z[lower+1]-z[lower])
        out.append((1-t)*values[lower]+t*values[lower+1])
    return np.array(out)


def intensity_slices(values, z, heights):
    """Return ``|psi|**2`` at requested heights for STM-style maps."""
    return np.abs(height_slices(values, z, heights))**2


def propagate_region(known, fixed, potential_eV, energy_eV, spacing_A,
                     tolerance=1e-9, maxiter=1000):
    """Dirichlet propagation outside an arbitrary fixed region (periodic xy).

    Interior fixed points supply signed wavefunction boundary values. The upper
    face is zero, the lower face is fixed. No imaginary energy or barrier clip.
    """
    u=np.asarray(potential_eV,dtype=float);known=np.asarray(known,dtype=np.complex128)
    fixed=np.asarray(fixed,dtype=bool)
    if u.shape!=known.shape or fixed.shape!=u.shape or u.ndim!=3 or u.shape[0]<3:
        raise ValueError('Known field, fixed mask and potential must share a 3D grid')
    if not fixed[0].all() or not fixed[-1].all():raise ValueError('Both z faces must be fixed')
    if not np.isfinite(u).all() or not np.isfinite(known).all():raise ValueError('Nonfinite field')
    hz,hy,hx=np.asarray(spacing_A)*ANG_TO_BOHR
    if min(hz,hy,hx)<=0:raise ValueError('Nonpositive spacing')
    active=~fixed[1:-1];shape=active.shape;n,ny,nx=shape
    if not active.any():raise ValueError('No vacuum points')
    q=2*EV_TO_HA*(u[1:-1]-energy_eV);allowed=bool(np.any(q[active]<=0))
    prescribed=np.where(fixed,known,0.);prescribed[-1]=0.
    scale=float(np.max(abs(prescribed)))
    if scale==0:return np.zeros_like(u),{'iterations':0,'relative_residual':0.,'solver':'zero'}
    prescribed/=scale
    rhs=(prescribed[2:]+prescribed[:-2])/hz**2
    rhs+=(np.roll(prescribed[1:-1],1,1)+np.roll(prescribed[1:-1],-1,1))/hy**2
    rhs+=(np.roll(prescribed[1:-1],1,2)+np.roll(prescribed[1:-1],-1,2))/hx**2
    rhs*=active
    def matvec(vector):
        v=vector.reshape(shape);a=v*active
        r=(2/hx**2+2/hy**2+2/hz**2+q)*a
        r-=(np.roll(a,1,1)+np.roll(a,-1,1))/hy**2
        r-=(np.roll(a,1,2)+np.roll(a,-1,2))/hx**2
        r[1:]-=a[:-1]/hz**2;r[:-1]-=a[1:]/hz**2
        return (r*active+v*(~active)).ravel()
    denom=4*np.sin(np.pi*np.arange(1,n+1)/(2*(n+1)))[:,None,None]**2/hz**2
    denom=denom+4*np.sin(np.pi*np.arange(ny)/ny)[None,:,None]**2/hy**2
    denom=denom+4*np.sin(np.pi*np.arange(nx)/nx)[None,None,:]**2/hx**2+max(float(q[active].mean()),1e-4)
    def precondition(vector):
        v=vector.reshape(shape);a=v*active
        transformed=fft2(dst(a,type=1,axis=0,norm='ortho'),axes=(1,2))
        value=idst(ifft2(transformed/denom,axes=(1,2)),type=1,axis=0,norm='ortho')
        return (value*active+v*(~active)).ravel()
    size=int(np.prod(shape));operator=LinearOperator((size,size),matvec=matvec,dtype=float)
    pre=LinearOperator((size,size),matvec=precondition,dtype=float);count=[0]
    def callback(_):count[0]+=1
    if allowed:
        try:
            solution,info=gmres(operator,rhs.ravel(),M=pre,rtol=tolerance,atol=0,restart=40,maxiter=maxiter,callback=callback,callback_type='legacy')
        except TypeError:
            solution,info=gmres(operator,rhs.ravel(),M=pre,tol=tolerance,atol=0,restart=40,maxiter=maxiter,callback=callback,callback_type='legacy')
    else:
        try:
            solution,info=cg(operator,rhs.ravel(),M=pre,rtol=tolerance,atol=0,maxiter=maxiter,callback=callback)
        except TypeError:
            solution,info=cg(operator,rhs.ravel(),M=pre,tol=tolerance,atol=0,maxiter=maxiter,callback=callback)
    residual=float(np.linalg.norm(matvec(solution)-rhs.ravel())/max(np.linalg.norm(rhs),1e-100))
    if info!=0 or residual>10*tolerance:raise RuntimeError('Region solve failed: info={}, residual={}'.format(info,residual))
    result=prescribed.copy();result[1:-1]+=solution.reshape(shape)*active;result*=scale
    return result,{'iterations':count[0],'relative_residual':residual,'solver':'gmres' if allowed else 'cg',
                   'allowed_fraction':float(np.mean(q[active]<=0)),'minimum_barrier_eV':float(q[active].min()/(2*EV_TO_HA)),
                   'vacuum_grid_points':int(active.sum())}
