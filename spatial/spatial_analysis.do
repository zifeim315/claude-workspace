*==============================================================================*
*  国家5A级景区建设、减污降碳与空间溢出效应 —— 完整可复现脚本（修正版）
*  一个 do + 一个 dta（data_spatial.dta 内含 lat lon mean_pgdp，Mata 现场构造权重矩阵）。
*  所有结果输出到 results 文件夹。
*  【修正】去掉所有 Mata 自定义函数（改为内联），避免 “'}' found where nothing expected”；
*          版本敏感命令(xthreg / re-Hausman)用 capture noisily 包裹，保证脚本不中断。
*
*  需安装： ssc install xsmle ;  ssc install xthreg ;  ssc install estout
*  运行前： cd "…/spatial"     （确保 data_spatial.dta 在当前目录）
*==============================================================================*
clear all
set more off
set matsize 2000
cap mkdir results

use "data_spatial.dta", clear
cap gen double lnfin = ln(fin)
xtset city_code year
global CTRL lnpgdp lndensity urban struc2 gov tech human lnfin

*------------------------------------------------------------------------------*
* 0. 在 Mata 用经纬度与人均GDP均值构造 5 类空间权重矩阵（全部内联，无自定义函数）
*------------------------------------------------------------------------------*
preserve
    bysort city_code (year): keep if _n==1
    sort city_code
    mkmat lat,       matrix(LAT)
    mkmat lon,       matrix(LON)
    mkmat mean_pgdp, matrix(YBAR)
restore

mata:
    LAT = st_matrix("LAT"):*(pi()/180)
    LON = st_matrix("LON"):*(pi()/180)
    Y   = st_matrix("YBAR")
    n   = rows(LAT)
    // 大圆距离矩阵 D (km)
    D = J(n,n,0)
    for (i=1;i<=n;i++) {
        for (j=1;j<=n;j++) {
            dla = LAT[j]-LAT[i]
            dlo = LON[j]-LON[i]
            a   = sin(dla/2)^2 + cos(LAT[i])*cos(LAT[j])*sin(dlo/2)^2
            D[i,j] = 2*6371*asin(sqrt(a))
        }
    }
    // 地理 1/d^2（对角=0），行标准化
    G = J(n,n,0)
    for (i=1;i<=n;i++) for (j=1;j<=n;j++) if (i!=j) G[i,j] = 1/(D[i,j]^2)
    rs = rowsum(G); rs = rs + (rs:==0); Wgeo = G :/ rs
    // 经济距离 1/|ΔȲ|
    E = J(n,n,0)
    for (i=1;i<=n;i++) for (j=1;j<=n;j++) if (i!=j & Y[i]!=Y[j]) E[i,j] = 1/abs(Y[i]-Y[j])
    rs = rowsum(E); rs = rs + (rs:==0); Wecon = E :/ rs
    // 经济地理嵌套 (1/d^2)*(Ȳ_j/Ȳ̄)
    Yr = Y :/ mean(Y)
    EG = J(n,n,0)
    for (i=1;i<=n;i++) for (j=1;j<=n;j++) if (i!=j) EG[i,j] = G[i,j]*Yr[j]
    rs = rowsum(EG); rs = rs + (rs:==0); Wegn = EG :/ rs
    // 经济地理权重 0.5*Wgeo + 0.5*Wecon，再行标准化
    EW = 0.5:*Wgeo + 0.5:*Wecon
    rs = rowsum(EW); rs = rs + (rs:==0); Wegw = EW :/ rs
    // 邻接：距离邻接 ≤163km；孤立城市取最近邻（全向量化，无循环/嵌套花括号）
    A     = (D:<=163) :* (D:>0)
    Dbig  = D + I(n):*1e12                       // 对角置为极大，排除自身
    rmins = rowmin(Dbig)                         // 每行最近邻距离
    NN    = (Dbig :== (rmins*J(1,n,1)))          // 最近邻位置指示矩阵
    iso   = (rowsum(A):==0)                      // 孤立行指示
    A     = A + NN:*(iso*J(1,n,1))               // 仅孤立行补最近邻
    A     = A + A'                               // 对称化
    A     = (A:>0)                               // 回到 0/1
    rs = rowsum(A); rs = rs + (rs:==0); Wadj = A :/ rs
    // 输出到 Stata 矩阵
    st_matrix("Wadj",  Wadj)
    st_matrix("Wgeo",  Wgeo)
    st_matrix("Wecon", Wecon)
    st_matrix("Wegw",  Wegw)
    st_matrix("Wegn",  Wegn)
    st_matrix("Dmat",  D)
end
di as result "== 5 类空间权重矩阵已构造：Wadj Wgeo Wecon Wegw Wegn =="

*------------------------------------------------------------------------------*
* 1. 空间诊断检验：全局 Moran's I（5 矩阵 × 双向FE残差）+ 分年度 Moran's I 表
*------------------------------------------------------------------------------*
sort year city_code
cap drop _res
qui reghdfe lnpoco2 DID $CTRL, a(city_code year) residuals(_res)
di as txt _n "== 全局 Moran's I（双向固定效应残差）=="
foreach W in Wadj Wgeo Wecon Wegw Wegn {
    mata:
        Wm = st_matrix("`W'"); e = st_data(.,"_res")
        N = 289; T = 21; num = 0; den = 0
        for (t=1;t<=T;t++) {
            et = e[((t-1)*N :+ (1::N))]; et = et :- mean(et)
            Wet = Wm*et; num = num + (et'Wet); den = den + (et'et)
        }
        st_numscalar("mI", num/den)
    end
    di as txt "  `W': Moran's I(残差) = " as res %6.3f mI
}
drop _res

* 分年度全局 Moran's I（主推矩阵 Wegw）+ 正态近似 Z 值 → results/
cap postclose MP
postfile MP int year double MoranI double Zscore double Pvalue using "results/moran_by_year.dta", replace
forvalues y = 2003/2023 {
    preserve
        keep if year==`y'
        sort city_code
        mata:
            W = st_matrix("Wegw"); x = st_data(.,"lnpoco2"); nn = rows(x); x = x :- mean(x)
            Wx = W*x; I = (x'Wx)/(x'x)
            S0 = sum(W); S1 = 0.5*sum((W+W'):^2); S2 = sum((rowsum(W)+colsum(W)'):^2)
            EI = -1/(nn-1); VI = (nn^2*S1 - nn*S2 + 3*S0^2)/(S0^2*(nn^2-1)) - EI^2
            z = (I-EI)/sqrt(VI)
            st_numscalar("mI", I); st_numscalar("mZ", z)
        end
        local pv = 2*(1-normal(abs(mZ)))
        post MP (`y') (mI) (mZ) (`pv')
    restore
}
postclose MP
preserve
    use "results/moran_by_year.dta", clear
    format MoranI Zscore Pvalue %9.3f
    list, sep(0) noobs
    cap export excel using "results/moran_by_year.xlsx", replace first(var)
    twoway (connected MoranI year, msymbol(O) lcolor(black) mcolor(black)), ///
        scheme(s1mono) ytitle("全局 Moran's I") xtitle("年份") ///
        title("减污降碳全局空间自相关的时间演变")
    graph export "results/fig_moran_trend.png", replace width(2000) height(1300)
restore

*------------------------------------------------------------------------------*
* 2. 空间杜宾模型 SDM：三列 (1)Time FE (2)Individual FE (3)Two-way FE（完整系数）
*    注：xsmle 对含空间滞后的模型默认报告 直接/间接/总 效应；nsim() 控制效应标准误模拟次数。
*------------------------------------------------------------------------------*
eststo clear
eststo sdm_time: xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) fe type(time) nsim(200)
estadd scalar rho = e(rho)
eststo sdm_ind:  xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) fe type(ind)  nsim(200)
estadd scalar rho = e(rho)
eststo sdm_both: xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) fe type(both) nsim(200)
estadd scalar rho = e(rho)
esttab sdm_time sdm_ind sdm_both using "results/table_SDM_3FE.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("(1) Time FE" "(2) Individual FE" "(3) Two-way FE") ///
    scalars("rho 空间自回归系数rho") stats(N, labels("观测值N")) ///
    nogaps compress title("空间杜宾模型完整估计(W_egw)")

* Hausman 检验（版本敏感，用 capture 包裹；失败时回退到非空间面板 Hausman）
cap noisily {
    qui xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) fe type(ind)
    est store fe_sdm
    qui xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wegw) re
    est store re_sdm
    hausman fe_sdm re_sdm, sigmamore
}
if _rc {
    di as txt "xsmle re 不可用，改用非空间面板 Hausman："
    qui xtreg lnpoco2 DID $CTRL, fe
    est store fe0
    qui xtreg lnpoco2 DID $CTRL, re
    est store re0
    hausman fe0 re0
}

*------------------------------------------------------------------------------*
* 3. 替换五类权重矩阵（完整系数，Time FE）
*------------------------------------------------------------------------------*
eststo clear
foreach W in Wadj Wgeo Wecon Wegw Wegn {
    eststo m_`W': xsmle lnpoco2 DID $CTRL, model(sdm) wmat(`W') fe type(time) nsim(200)
    estadd scalar rho = e(rho)
}
esttab m_Wadj m_Wgeo m_Wecon m_Wegw m_Wegn using "results/table_SDM_5matrices.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("邻接" "地理" "经济" "经济地理权重" "经济地理嵌套") ///
    scalars("rho rho") stats(N, labels("N")) nogaps compress ///
    title("五类权重矩阵 SDM 完整估计(Time FE)")

*------------------------------------------------------------------------------*
* 4. 溢出范围：不同 km 距离阈值下的 SDM（Time FE）→ 表 + 图
*------------------------------------------------------------------------------*
cap postclose PM
postfile PM int km double rho using "results/spillover_range.dta", replace
foreach km in 150 200 250 300 350 400 450 500 {
    mata:
        Dm = st_matrix("Dmat"); n = rows(Dm)
        Bk    = (Dm:<=`km') :* (Dm:>0)
        Dbig  = Dm + I(n):*1e12
        rmins = rowmin(Dbig)
        NN    = (Dbig :== (rmins*J(1,n,1)))
        iso   = (rowsum(Bk):==0)
        Bk    = Bk + NN:*(iso*J(1,n,1))
        Bk    = Bk + Bk'
        Bk    = (Bk:>0)
        rs = rowsum(Bk); rs = rs + (rs:==0); st_matrix("Wk", Bk:/rs)
    end
    qui xsmle lnpoco2 DID $CTRL, model(sdm) wmat(Wk) fe type(time) nsim(1)
    post PM (`km') (e(rho))
}
postclose PM
preserve
    use "results/spillover_range.dta", clear
    twoway (connected rho km, msymbol(O) lcolor(black) mcolor(black)), ///
        yline(0,lpattern(dash) lcolor(gs9)) scheme(s1mono) ///
        xtitle("地理距离阈值 (km)") ytitle("空间自回归系数 rho") ///
        title("减污降碳空间溢出的地理衰减")
    graph export "results/fig_distance_decay.png", replace width(2000) height(1300)
restore

*------------------------------------------------------------------------------*
* 5. 门槛/分段异质性溢出
*------------------------------------------------------------------------------*
* (A) Hansen 面板门槛（版本敏感，capture 包裹）
cap noisily xthreg lnpoco2 DID $CTRL, rx(DID) qx(lnpgdp) thnum(2) trim(0.01 0.01) grid(300) bs(300 300)

* (B) 分区制 SDM（各维度分高/低，Time FE）
eststo clear
foreach v of newlist human lnpgdp ter_gdp er {
    capture confirm variable `v'
    if _rc continue
    cap drop hi_`v' DIDlo_`v' DIDhi_`v' m_`v'
    bysort city_code: egen m_`v' = mean(`v')
    qui sum m_`v', detail
    gen byte hi_`v' = m_`v' > r(p50)
    gen double DIDlo_`v' = DID*(1-hi_`v')
    gen double DIDhi_`v' = DID*hi_`v'
    eststo reg_`v': xsmle lnpoco2 DIDlo_`v' DIDhi_`v' $CTRL, model(sdm) wmat(Wegw) fe type(time) nsim(100)
}
cap drop east DIDlo_east DIDhi_east
gen byte east = region=="东部"
gen double DIDlo_east = DID*(1-east)
gen double DIDhi_east = DID*east
eststo reg_east: xsmle lnpoco2 DIDlo_east DIDhi_east $CTRL, model(sdm) wmat(Wegw) fe type(time) nsim(100)
esttab reg_* using "results/table_regime_spillover.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) nogaps compress ///
    title("分区制异质性溢出(Time FE, W_egw)")

* (C) 控制组溢出污染检验（queen 邻接）
sort year city_code
cap drop spill
mata:
    Wb = (st_matrix("Wadj"):>0); did = st_data(.,"DID")
    N = rows(Wb); T = rows(did)/N; sp = J(rows(did),1,0)
    for (t=1;t<=T;t++) {
        idx = (t-1)*N :+ (1::N); dt = did[idx]
        nb = (Wb*dt) :> 0; sp[idx] = (dt:==0) :* nb
    }
    st_store(., st_addvar("byte","spill"), sp)
end
eststo clear
eststo base_did:  reghdfe lnpoco2 DID $CTRL,       a(city_code year) vce(cluster city_code)
eststo spill_did: reghdfe lnpoco2 DID spill $CTRL, a(city_code year) vce(cluster city_code)
esttab base_did spill_did using "results/table_spillover_contamination.rtf", replace ///
    b(%9.3f) t(%9.3f) star(* 0.1 ** 0.05 *** 0.01) keep(DID spill) ///
    mtitles("基准" "加溢出虚拟变量") nogaps compress title("控制组溢出污染检验")

*------------------------------------------------------------------------------*
* 6. Moran 散点图：2003 与 2023 年（主推矩阵 Wegw）→ results/
*------------------------------------------------------------------------------*
foreach yr in 2003 2023 {
    preserve
        keep if year==`yr'
        sort city_code
        qui sum lnpoco2
        cap drop z wz1
        gen double z = (lnpoco2 - r(mean))/r(sd)
        mkmat z, matrix(zz)
        mata:
            W = st_matrix("Wegw"); zv = st_matrix("zz")
            st_matrix("Wz", W*zv)
            st_numscalar("MI", (zv'*(W*zv))/(zv'zv))
        end
        svmat Wz, names(wz)
        local mi = MI
        twoway (scatter wz1 z, mcolor(black) msize(small) msymbol(Oh)) ///
               (lfit wz1 z, lcolor(black)), ///
            yline(0,lpattern(dash) lcolor(gs10)) xline(0,lpattern(dash) lcolor(gs10)) ///
            scheme(s1mono) legend(off) xtitle("标准化 lnpoco2") ytitle("空间滞后 W·z") ///
            title("Moran 散点图 `yr'  (Moran's I=" + string(`mi',"%5.3f") + ")")
        graph export "results/fig_moran_`yr'.png", replace width(1600) height(1400)
    restore
}

di as result _n "================ 完成：全部结果见 results 文件夹 ================"
