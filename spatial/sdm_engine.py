import numpy as np, pandas as pd, pyreadstat
from scipy.optimize import minimize_scalar
from scipy import sparse
SCR="/tmp/claude-0/-home-user-claude-workspace/f562582e-c817-5f44-a99c-4b38ba34c0b7/scratchpad"

def load_panel():
    df,_=pyreadstat.read_dta("/home/user/claude-workspace/data_integrated.dta")
    df['lnfin']=np.log(df['fin'])
    df=df.sort_values(['year','city_code']).reset_index(drop=True)
    return df

def twoway_demean(M, cid, tid):
    M=np.asarray(M,float); orig1d=(M.ndim==1)
    if orig1d: M=M.reshape(-1,1)
    for _ in range(2):
        M=M-pd.DataFrame(M).groupby(cid).transform('mean').values
        M=M-pd.DataFrame(M).groupby(tid).transform('mean').values
    return M[:,0] if orig1d else M

def sdm_ml(y, X, W, cid, tid, eig=None):
    # y (NT,), X (NT,k) raw stacked by (year,city); W single-year NxN; block-diag over T
    N=W.shape[0]; NT=len(y); T=NT//N
    if eig is None: eig=np.linalg.eigvals(W).real
    # spatial lags via block-diagonal multiply (reshape T x N)
    def wlag(v):
        v2=v.reshape(T,N).T                # N x T
        return (W@v2).T.reshape(NT, *v.shape[1:]) if v.ndim>1 else (W@v2).T.reshape(NT)
    Wy=wlag(y)
    WX=np.column_stack([wlag(X[:,j]) for j in range(X.shape[1])])
    Z=np.column_stack([X,WX])             # Durbin design
    yt=twoway_demean(y,cid,tid); Wyt=twoway_demean(Wy,cid,tid); Zt=twoway_demean(Z,cid,tid)
    def negll(rho):
        e=(yt-rho*Wyt)
        b,_,_,_=np.linalg.lstsq(Zt,e,rcond=None)
        res=e-Zt@b; s2=(res@res)/NT
        jac=np.sum(np.log(np.abs(1-rho*eig)))*T
        return -(-NT/2*np.log(s2)+jac)
    lo,hi=1/eig.min()+1e-4, 1/eig.max()-1e-4
    r=minimize_scalar(negll,bounds=(lo,hi),method='bounded')
    rho=r.x; ll=-r.fun
    ll0=-negll(0.0)
    LR=2*(ll-ll0)
    # betas at rho
    e=(yt-rho*Wyt); b,_,_,_=np.linalg.lstsq(Zt,e,rcond=None)
    res=e-Zt@b; s2=(res@res)/NT
    XtX_inv=np.linalg.inv(Zt.T@Zt)
    seb=np.sqrt(np.diag(s2*XtX_inv))
    return dict(rho=rho, LR=LR, ll=ll, beta=b, se=seb, eig=eig, s2=s2)

def moran_resid(y,X,W,cid,tid):
    yt=twoway_demean(y,cid,tid); Xt=twoway_demean(X,cid,tid)
    b,_,_,_=np.linalg.lstsq(Xt,yt,rcond=None); r=yt-Xt@b
    N=W.shape[0]; T=len(y)//N; r2=r.reshape(T,N).T
    Wr=W@r2; num=(r2*Wr).sum(); den=(r2*r2).sum()
    return num/den  # row-std W => Moran-like

def oneway_demean(M, cid):
    M=np.asarray(M,float); o=(M.ndim==1)
    if o: M=M.reshape(-1,1)
    M=M-pd.DataFrame(M).groupby(cid).transform('mean').values
    return M[:,0] if o else M

def sdm_ml_fe(y, X, W, cid, tid, fe='twoway', durbin=True, eig=None):
    N=W.shape[0]; NT=len(y); T=NT//N
    if eig is None: eig=np.linalg.eigvals(W).real
    def wlag(v):
        v2=v.reshape(T,N).T
        return (W@v2).T.reshape(v.shape) if v.ndim>1 else (W@v2).T.reshape(NT)
    Wy=wlag(y)
    if durbin:
        WX=np.column_stack([wlag(X[:,j]) for j in range(X.shape[1])]); Z=np.column_stack([X,WX])
    else: Z=X
    dm=(lambda M: twoway_demean(M,cid,tid)) if fe=='twoway' else ((lambda M: oneway_demean(M,tid)) if fe=='time' else (lambda M: oneway_demean(M,cid)))
    yt=dm(y); Wyt=dm(Wy); Zt=dm(Z)
    def negll(rho):
        e=yt-rho*Wyt; b,_,_,_=np.linalg.lstsq(Zt,e,rcond=None); r=e-Zt@b; s2=(r@r)/NT
        return -(-NT/2*np.log(s2)+np.sum(np.log(np.abs(1-rho*eig)))*T)
    from scipy.optimize import minimize_scalar
    lo,hi=1/eig.min()+1e-4,1/eig.max()-1e-4
    r=minimize_scalar(negll,bounds=(lo,hi),method='bounded'); rho=r.x
    LR=2*(-r.fun-(-negll(0.0)))
    e=yt-rho*Wyt; b,_,_,_=np.linalg.lstsq(Zt,e,rcond=None); res=e-Zt@b; s2=(res@res)/NT
    se=np.sqrt(np.diag(s2*np.linalg.inv(Zt.T@Zt)))
    return dict(rho=rho,LR=LR,beta=b,se=se,k=X.shape[1],durbin=durbin)
